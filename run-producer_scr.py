#!/usr/bin/env python3
"""Send one long-term MONICA ensemble run per SCR point and parameter group."""

from __future__ import annotations

import copy
import gzip
import json
import os
import posixpath
import re
import sqlite3
import sys
import time
from collections import defaultdict
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np
import zmq
from pyproj import CRS

import monica_io3
import monica_run_lib as Mrunlib
import soil_io3
from scr_simulation_lib import (
    PIPELINE_ID,
    as_bool,
    file_sha256,
    is_complete_output,
    limited,
    parse_key_value_args,
    read_parameter_groups,
    read_points,
    scenario_id,
    scenario_output_path,
)


PATHS = {
    "localProducer-localMonica": {
        "path-to-climate-dir": "./data/",
        "monica-path-to-climate-dir": "./data/",
        "path-to-data-dir": "./data/",
    },
    "re-local-remote": {
        "path-to-climate-dir": "data/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        "path-to-data-dir": "./data/",
    },
    "mbm-local-remote": {
        "path-to-climate-dir": "/run/user/1000/gvfs/sftp:host=login01.cluster.zalf.de,user=rpm/beegfs/common/data/climate/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        "path-to-data-dir": "./data/",
    },
    "remoteProducer-remoteMonica": {
        "path-to-climate-dir": "/data/",
        "monica-path-to-climate-dir": "/monica_data/climate-data/",
        "path-to-data-dir": "./data/",
    },
}

DATA_SOIL_DB = "germany/buek200.sqlite"
DATA_GRID_SOIL = "germany/buek200_1000_25832_etrs89-utm32n.asc"
DATA_GRID_HEIGHT = "germany/dem_1000_25832_etrs89-utm32n.asc"
DATA_GRID_SLOPE = "germany/slope_1000_25832_etrs89-utm32n.asc"
TEMPLATE_PATH_HARVEST = "projects/monica-germany/ILR_SEED_HARVEST_doys_{crop_id}.csv"

DEFAULTS = {
    "mode": "remoteProducer-remoteMonica",
    "server": "login01.cluster.zalf.de",
    "server-port": "6666",
    "shared_id": "scr-lhs",
    "setups-file": "sim_setups_scr.csv",
    "setup-id": "1",
    "points-file": "scr_points.csv",
    "params-file": "scr_parameter_groups.csv",
    "out": "scr-csv-out",
    "start-date": "1980-01-01",
    "end-date": "2020-12-31",
    "expected-seasons": "40",
    "max-points": "-1",
    "max-params": "-1",
    "resume": "true",
    "check-climate-files": "false",
    "flat-climate-dir": "",
    "preflight-only": "false",
    "progress-every": "25",
}


def load_grid(path: Path, dtype: Any) -> tuple[dict[str, float], np.ndarray]:
    metadata, _ = Mrunlib.read_header(str(path))
    return metadata, np.loadtxt(str(path), dtype=dtype, skiprows=6)


def grid_value_at_xy(grid: np.ndarray, metadata: dict[str, float], x: float, y: float) -> float:
    cellsize = float(metadata["cellsize"])
    col = int((x - float(metadata["xllcorner"])) // cellsize)
    row_from_bottom = int((y - float(metadata["yllcorner"])) // cellsize)
    row = int(metadata["nrows"]) - 1 - row_from_bottom
    if row < 0 or col < 0 or row >= grid.shape[0] or col >= grid.shape[1]:
        raise ValueError(f"Coordinate ({x}, {y}) lies outside grid")
    value = float(grid[row, col])
    if np.isclose(value, float(metadata["nodata_value"])):
        raise ValueError(f"Coordinate ({x}, {y}) is nodata in grid")
    return value


def yearless_date(value: str, year: int = 2001) -> date:
    _, month, day = (int(part) for part in value.split("-"))
    return date(year, month, day)


def shifted_sowing_date(base: str, earliest: str, latest: str, shift_days: float) -> str:
    base_date = yearless_date(base)
    earliest_date = yearless_date(earliest)
    latest_date = yearless_date(latest)
    shifted = base_date + timedelta(days=int(round(shift_days)))

    if earliest_date <= latest_date:
        shifted = min(max(shifted, earliest_date), latest_date)
    else:
        # Circular windows are unusual for German sowing dates, but handle them safely.
        if latest_date < shifted < earliest_date:
            to_earliest = (earliest_date - shifted).days
            to_latest = (shifted - latest_date).days
            shifted = earliest_date if to_earliest <= to_latest else latest_date
    return f"0000-{shifted.month:02d}-{shifted.day:02d}"


def adjusted_latest_harvest(seed_data: dict[str, Any], sowing_date: str) -> str:
    sowing_doy = yearless_date(sowing_date).timetuple().tm_yday
    latest_harvest_doy = int(seed_data["latest-harvest-doy"])
    harvest_doy = min(latest_harvest_doy, sowing_doy - 1)
    harvest = date(2001, 1, 1) + timedelta(days=harvest_doy - 1)
    return f"0001-{harvest.month:02d}-{harvest.day:02d}"


def climate_subpath(setup: dict[str, Any], crow: int, ccol: int) -> str:
    parts = []
    for key in ("gcm", "rcm", "scenario", "ensmem", "version"):
        value = str(setup.get(key, "")).strip("/\\ ")
        if value:
            parts.append(value)
    parts.extend([str(crow), f"daily_mean_RES1_C{ccol}R{crow}.csv.gz"])
    return "/".join(parts)


DATE_PATTERN = re.compile(r"(?<!\d)((?:19|20)\d{2})[-/]?(\d{2})[-/]?(\d{2})(?!\d)")


def climate_date_extent(path: Path) -> tuple[date, date]:
    opener = gzip.open if path.suffix == ".gz" else open
    first = None
    last = None
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = DATE_PATTERN.search(line)
            if not match:
                continue
            parsed = date(*(int(part) for part in match.groups()))
            first = parsed if first is None else first
            last = parsed
    if first is None or last is None:
        raise ValueError("no YYYY-MM-DD or YYYYMMDD dates found")
    return first, last


def validate_climate_file(path: Path, start: date, end: date) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)
    first, last = climate_date_extent(path)
    if first > start or last < end:
        raise ValueError(f"coverage is {first} to {last}, expected at least {start} to {end}")


def load_base_environment(setup: dict[str, Any], config: dict[str, Any]) -> tuple[dict, dict]:
    with open(setup["sim.json"], encoding="utf-8") as handle:
        sim_json = json.load(handle)
    with open(setup["site.json"], encoding="utf-8") as handle:
        site_json = json.load(handle)
    with open(setup["crop.json"], encoding="utf-8") as handle:
        crop_json = json.load(handle)

    sim_json["climate.csv-options"]["start-date"] = config["start-date"]
    sim_json["climate.csv-options"]["end-date"] = config["end-date"]
    sim_json["output"]["events"] = sim_json["output"]["daily-yields"]
    sim_json["output"]["obj-outputs?"] = not setup["nc_mode"] and not setup["bgr"]

    crop_json["CropParameters"]["__enable_vernalisation_factor_fix__"] = setup.get(
        "use_vernalisation_fix", False
    )
    crop_json["cropRotation"][2] = setup["crop-id"]

    environment = monica_io3.create_env_json_from_json_config(
        {"crop": crop_json, "site": site_json, "sim": sim_json, "climate": ""}
    )
    environment["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
    if config["shared_id"]:
        environment["sharedId"] = config["shared_id"]

    environment["params"]["userCropParameters"][
        "__enable_T_response_leaf_expansion__"
    ] = setup["LeafExtensionModifier"]
    simulation = environment["params"]["simulationParameters"]
    simulation["UseNMinMineralFertilisingMethod"] = setup["fertilization"]
    simulation["NitrogenResponseOn"] = setup["NitrogenResponseOn"]
    simulation["WaterDeficitResponseOn"] = setup["WaterDeficitResponseOn"]
    simulation["EmergenceMoistureControlOn"] = setup["EmergenceMoistureControlOn"]
    simulation["EmergenceFloodingControlOn"] = setup["EmergenceFloodingControlOn"]
    return environment, sim_json


def normalise_interpolated_id(value: Any) -> int:
    array = np.asarray(value)
    if array.size != 1:
        raise ValueError(f"Expected one interpolated id, got {value}")
    return int(array.item())


def prepare_point_contexts(
    config: dict[str, Any],
    paths: dict[str, str],
    setup: dict[str, Any],
    points: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    data_dir = Path(paths["path-to-data-dir"])
    soil_meta, soil_grid = load_grid(data_dir / DATA_GRID_SOIL, np.int32)
    dem_meta, dem_grid = load_grid(data_dir / DATA_GRID_HEIGHT, float)
    slope_meta, slope_grid = load_grid(data_dir / DATA_GRID_SLOPE, float)

    wgs84 = CRS.from_epsg(4326)
    soil_crs = CRS.from_epsg(25832)
    climate_map = (
        Path(paths["path-to-climate-dir"])
        / str(setup["climate_path_to_latlon_file"])
        / "latlon_to_rowcol.json"
    )
    climate_coords: dict[tuple[int, int], tuple[float, float]] = {}
    climate_interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(
        str(climate_map), wgs84, soil_crs, climate_coords
    )

    crop_id = str(setup["crop-id"]).split("_")[0]
    ilr_data: defaultdict = defaultdict(
        lambda: {"interpolate": None, "data": defaultdict(dict), "is-winter-crop": None}
    )
    ilr_path = data_dir / TEMPLATE_PATH_HARVEST.format(crop_id=crop_id)
    Mrunlib.create_seed_harvest_geoGrid_interpolator_and_read_data(
        str(ilr_path), wgs84, soil_crs, ilr_data
    )

    start = date.fromisoformat(str(config["start-date"]))
    end = date.fromisoformat(str(config["end-date"]))
    check_climate = as_bool(config["check-climate-files"])
    checked_climate_files: set[Path] = set()
    soil_cache: dict[int, list[dict]] = {}
    contexts = []
    errors = []

    with sqlite3.connect(str(data_dir / DATA_SOIL_DB)) as soil_connection:
        for point in points:
            try:
                srow = int(point["srow"])
                scol = int(point["scol"])
                soil_id = int(point["soil_id"])
                grid_soil_id = int(soil_grid[srow, scol])
                if grid_soil_id != soil_id:
                    raise ValueError(f"soil CSV={soil_id}, grid={grid_soil_id}")

                if soil_id not in soil_cache:
                    soil_cache[soil_id] = soil_io3.soil_parameters(soil_connection, soil_id)
                if not soil_cache[soil_id]:
                    raise ValueError(f"soil id {soil_id} has no profile in buek200.sqlite")

                x = float(point["x"])
                y = float(point["y"])
                crow, ccol = climate_interpolator(x, y)
                crow = normalise_interpolated_id(crow)
                ccol = normalise_interpolated_id(ccol)
                climate_lat, climate_lon = climate_coords[(crow, ccol)]

                ilr_station = normalise_interpolated_id(ilr_data[crop_id]["interpolate"](x, y))
                seed_data = ilr_data[crop_id]["data"][ilr_station]
                if not seed_data:
                    raise ValueError(f"ILR station {ilr_station} has no {crop_id} dates")

                subpath = climate_subpath(setup, crow, ccol)
                flat_climate_dir = str(config["flat-climate-dir"]).strip()
                if flat_climate_dir:
                    host_climate_path = Path(flat_climate_dir) / Path(subpath).name
                else:
                    host_climate_path = (
                        Path(paths["path-to-climate-dir"])
                        / str(setup["climate_path_to_csvs"])
                        / Path(subpath)
                    )
                if check_climate and host_climate_path not in checked_climate_files:
                    validate_climate_file(host_climate_path, start, end)
                    checked_climate_files.add(host_climate_path)

                if flat_climate_dir:
                    monica_climate_path = str(host_climate_path.resolve())
                else:
                    monica_climate_path = posixpath.join(
                        paths["monica-path-to-climate-dir"].rstrip("/"),
                        str(setup["climate_path_to_csvs"]).strip("/"),
                        subpath,
                    )
                contexts.append(
                    {
                        **point,
                        "soil_profile": copy.deepcopy(soil_cache[soil_id]),
                        "height_nn": grid_value_at_xy(dem_grid, dem_meta, x, y),
                        "slope": grid_value_at_xy(slope_grid, slope_meta, x, y),
                        "crow": crow,
                        "ccol": ccol,
                        "climate_lat": float(climate_lat),
                        "climate_lon": float(climate_lon),
                        "host_climate_path": str(host_climate_path),
                        "monica_climate_path": monica_climate_path,
                        "base_sowing_date": seed_data["sowing-date"],
                        "earliest_sowing_date": seed_data["earliest-sowing-date"],
                        "latest_sowing_date": seed_data["latest-sowing-date"],
                        "seed_data": dict(seed_data),
                    }
                )
            except Exception as exc:
                errors.append(f"BKR {point['bkr_id']} row/col {point['srow']}/{point['scol']}: {exc}")

    if errors:
        joined = "\n  - ".join(errors)
        raise RuntimeError(f"SCR preflight failed:\n  - {joined}")
    return contexts


def set_quantity(container: dict, key: str, value: float) -> None:
    current = container.get(key)
    if isinstance(current, list) and len(current) == 2 and isinstance(current[1], str):
        current[0] = float(value)
    else:
        container[key] = float(value)


def apply_site_parameters(environment: dict, point: dict[str, Any], setup: dict[str, Any]) -> None:
    site = environment["params"]["siteParameters"]
    profile = copy.deepcopy(point["soil_profile"])
    for layer in profile:
        layer["SoilNitrate"] = 0
        layer["SoilAmmonium"] = 0
    site["SoilProfileParameters"] = profile

    if setup["elevation"]:
        set_quantity(site, "HeightNN", point["height_nn"])
    if setup["slope"]:
        site["Slope"] = float(point["slope"]) / 100.0
    if setup["latitude"]:
        site["Latitude"] = float(point["climate_lat"])

    environment_parameters = environment["params"]["userEnvironmentParameters"]
    if setup["groundwater-level"]:
        groundwater_depth = 20.0
        layer_depth = 0.0
        for layer in profile:
            if layer.get("is_in_groundwater", False):
                groundwater_depth = layer_depth
                break
            layer_depth += float(Mrunlib.get_value(layer["Thickness"]))
        environment_parameters["MinGroundwaterDepthMonth"] = 3
        environment_parameters["MinGroundwaterDepth"] = [max(0.0, groundwater_depth - 0.2), "m"]
        environment_parameters["MaxGroundwaterDepth"] = [groundwater_depth + 0.2, "m"]

    if setup["impenetrable-layer"]:
        impenetrable_depth = float(Mrunlib.get_value(environment_parameters["LeachingDepth"]))
        layer_depth = 0.0
        for layer in profile:
            if layer.get("is_impenetrable", False):
                impenetrable_depth = layer_depth
                break
            layer_depth += float(Mrunlib.get_value(layer["Thickness"]))
        environment_parameters["LeachingDepth"] = [impenetrable_depth, "m"]
        site["ImpenetrableLayerDepth"] = [impenetrable_depth, "m"]


def apply_crop_parameters(crop_parameters: dict, parameters: dict[str, Any]) -> None:
    cultivar = crop_parameters["cultivar"]
    species = crop_parameters["species"]

    stage_sums = cultivar["StageTemperatureSum"][0]
    for stage in range(5):
        stage_sums[stage] = float(parameters[f"stage_temperature_sum_{stage}"])
    cultivar["VernalisationRequirement"][1] = float(parameters["vernalisation_requirement"])
    for stage in (1, 2, 3):
        cultivar["DaylengthRequirement"][0][stage] = float(parameters["daylength_requirement"])
    cultivar["MaxAssimilationRate"] = float(parameters["max_assimilation_rate"])
    cultivar["SpecificLeafArea"][0] = [
        float(value) * float(parameters["specific_leaf_area_scale"])
        for value in cultivar["SpecificLeafArea"][0]
    ]
    species["RootPenetrationRate"] = float(parameters["root_penetration_rate"])
    cultivar["CropSpecificMaxRootingDepth"] = float(parameters["max_rooting_depth"])
    shift = float(parameters["drought_threshold_shift"])
    cultivar["DroughtStressThreshold"] = [
        min(1.0, max(0.0, float(value) + shift))
        for value in cultivar["DroughtStressThreshold"]
    ]


def apply_management(environment: dict, point: dict[str, Any], parameters: dict[str, Any]) -> tuple[str, str]:
    worksteps = environment["cropRotation"][0]["worksteps"]
    sowing = next(step for step in worksteps if step["type"].endswith("Sowing"))
    harvest = next(step for step in worksteps if step["type"].endswith("Harvest"))
    sowing_date = shifted_sowing_date(
        point["base_sowing_date"],
        point["earliest_sowing_date"],
        point["latest_sowing_date"],
        parameters["sowing_shift_days"],
    )
    latest_harvest = adjusted_latest_harvest(point["seed_data"], sowing_date)
    sowing["date"] = sowing_date
    harvest["latest-date"] = latest_harvest

    crop_parameters = sowing["crop"]["cropParams"]
    apply_crop_parameters(crop_parameters, parameters)

    nitrogen_scale = float(parameters["n_fertilizer_scale"])
    for step in worksteps:
        if step["type"].endswith("MineralFertilization"):
            step["amount"][0] = float(step["amount"][0]) * nitrogen_scale

    simulation = environment["params"]["simulationParameters"]
    nmin = simulation["NMinUserParams"]
    nmin["min"] = float(nmin["min"]) * nitrogen_scale
    nmin["max"] = float(nmin["max"]) * nitrogen_scale
    simulation["UseAutomaticIrrigation"] = bool(parameters["irrigation_enabled"])
    irrigation = simulation["AutoIrrigationParams"]
    irrigation["amount"][0] = float(parameters["irrigation_amount_mm"])
    irrigation["threshold"] = float(parameters["irrigation_threshold"])
    return sowing_date, latest_harvest


def build_environment(
    base_environment: dict,
    point: dict[str, Any],
    parameters: dict[str, Any],
    setup: dict[str, Any],
    config: dict[str, Any],
) -> dict:
    environment = copy.deepcopy(base_environment)
    apply_site_parameters(environment, point, setup)
    sowing_date, latest_harvest = apply_management(environment, point, parameters)
    environment["pathToClimateCSV"] = [point["monica_climate_path"]]

    start_year = date.fromisoformat(str(config["start-date"])).year
    expected_seasons = int(config["expected-seasons"])
    sid = scenario_id(point, parameters)
    environment["customId"] = {
        "pipeline": PIPELINE_ID,
        "scenario_id": sid,
        "setup_id": int(config["setup-id"]),
        "bkr_id": int(point["bkr_id"]),
        "bkr_name": point["name"],
        "param_id": parameters["param_id"],
        "srow": int(point["srow"]),
        "scol": int(point["scol"]),
        "crow": int(point["crow"]),
        "ccol": int(point["ccol"]),
        "soil_id": int(point["soil_id"]),
        "expected_seasons": expected_seasons,
        "first_growing_year": start_year + 1,
        "last_growing_year": start_year + expected_seasons,
        "sowing_date": sowing_date,
        "latest_harvest_date": latest_harvest,
        "sowing_shift_days": int(round(float(parameters["sowing_shift_days"]))),
        "n_fertilizer_scale": float(parameters["n_fertilizer_scale"]),
        "irrigation_enabled": bool(parameters["irrigation_enabled"]),
        "irrigation_amount_mm": float(parameters["irrigation_amount_mm"]),
        "irrigation_threshold": float(parameters["irrigation_threshold"]),
    }
    return environment


def atomic_write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def create_run_signature(config: dict[str, Any], setup: dict[str, Any]) -> dict[str, Any]:
    source_paths = {
        "points_sha256": config["points-file"],
        "parameters_sha256": config["params-file"],
        "setups_sha256": config["setups-file"],
        "sim_sha256": setup["sim.json"],
        "crop_sha256": setup["crop.json"],
        "site_sha256": setup["site.json"],
    }
    return {
        "pipeline": PIPELINE_ID,
        "setup_id": int(config["setup-id"]),
        "start_date": config["start-date"],
        "end_date": config["end-date"],
        "expected_seasons": int(config["expected-seasons"]),
        "flat_climate_dir": str(config["flat-climate-dir"]),
        **{key: file_sha256(path) for key, path in source_paths.items()},
    }


def ensure_compatible_output_directory(output_dir: Path, signature: dict[str, Any]) -> None:
    manifest_path = output_dir / "run_manifest.json"
    completed_outputs = list(output_dir.glob("BKR_*/BKR*_P*_daily.csv")) if output_dir.exists() else []
    if not manifest_path.is_file():
        if completed_outputs:
            raise RuntimeError(
                f"{output_dir} contains scenario CSV files but no run_manifest.json; "
                "use a new output directory"
            )
        return
    with manifest_path.open("r", encoding="utf-8") as handle:
        existing = json.load(handle).get("run_signature", {})
    differences = [key for key, value in signature.items() if existing.get(key) != value]
    if differences and completed_outputs:
        raise RuntimeError(
            f"Output directory {output_dir} belongs to different inputs ({', '.join(differences)}); "
            "use a new out= directory"
        )


def main() -> None:
    config = parse_key_value_args(DEFAULTS, sys.argv[1:])
    if config["mode"] not in PATHS:
        raise ValueError(f"Unknown mode: {config['mode']}")
    paths = PATHS[config["mode"]]
    points = limited(read_points(config["points-file"]), config["max-points"])
    parameter_groups = limited(read_parameter_groups(config["params-file"]), config["max-params"])

    setups = Mrunlib.read_sim_setups(config["setups-file"])
    setup_id = int(config["setup-id"])
    if setup_id not in setups:
        raise ValueError(f"Setup {setup_id} not found in {config['setups-file']}")
    setup = setups[setup_id]
    output_dir = Path(config["out"])
    run_signature = create_run_signature(config, setup)
    ensure_compatible_output_directory(output_dir, run_signature)

    print(f"Preparing {len(points)} SCR points and {len(parameter_groups)} parameter groups")
    point_contexts = prepare_point_contexts(config, paths, setup, points)
    base_environment, _ = load_base_environment(setup, config)

    resume = as_bool(config["resume"])
    pending = []
    skipped = 0
    for point in point_contexts:
        for parameters in parameter_groups:
            output_path = scenario_output_path(config["out"], point, parameters)
            if resume and is_complete_output(output_path):
                skipped += 1
            else:
                pending.append((point, parameters))

    manifest = {
        "pipeline": PIPELINE_ID,
        "run_signature": run_signature,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "setup_id": setup_id,
        "start_date": config["start-date"],
        "end_date": config["end-date"],
        "expected_seasons_per_scenario": int(config["expected-seasons"]),
        "points": len(point_contexts),
        "parameter_groups": len(parameter_groups),
        "total_scenarios": len(point_contexts) * len(parameter_groups),
        "skipped_complete": skipped,
        "pending_scenarios": len(pending),
        "check_climate_files": as_bool(config["check-climate-files"]),
        "scenarios": [scenario_id(point, parameters) for point, parameters in pending],
    }
    atomic_write_json(output_dir / "run_manifest.json", manifest)

    print(f"Preflight passed for {len(point_contexts)} SCR points")
    print(f"Scenarios total={manifest['total_scenarios']} skipped={skipped} pending={len(pending)}")
    if as_bool(config["preflight-only"]):
        print("Preflight-only mode: no MONICA environments were sent")
        return
    if not pending:
        print("All scenario CSV files already exist; nothing to send")
        return

    context = zmq.Context()
    socket = context.socket(zmq.PUSH)
    socket.setsockopt(zmq.LINGER, 10000)
    socket.connect(f"tcp://{config['server']}:{config['server-port']}")
    progress_every = max(1, int(config["progress-every"]))
    started = time.perf_counter()
    try:
        for index, (point, parameters) in enumerate(pending, start=1):
            environment = build_environment(base_environment, point, parameters, setup, config)
            socket.send_json(environment)
            if index % progress_every == 0 or index == len(pending):
                elapsed = time.perf_counter() - started
                print(f"Sent {index}/{len(pending)} environments in {elapsed:.1f} s")
    finally:
        socket.close()
        context.term()


if __name__ == "__main__":
    main()
