#!/usr/bin/env python3
"""Create reproducible SCR representative points and LHS parameter groups."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import numpy as np
import shapefile
from matplotlib.path import Path as PolygonPath
from pyproj import CRS, Transformer
from scipy.stats import qmc

import monica_run_lib as Mrunlib
from scr_simulation_lib import as_bool, parse_key_value_args


DEFAULTS = {
    "shape-file": "data/germany/BKR_fixed_utm32.shp",
    "soil-grid": "data/germany/buek200_1000_25832_etrs89-utm32n.asc",
    "crop-grid": "data/germany/crops-all2017-2019_1000_25832_etrs89-utm32n.asc",
    "points-file": "scr_points.csv",
    "params-file": "scr_parameter_groups.csv",
    "parameter-groups": "20",
    "lhs-seed": "42",
    "include-baseline": "false",
}


PARAMETER_SPECS = [
    ("stage_temperature_sum_0", 116.8, 175.2),
    ("stage_temperature_sum_1", 220.8, 331.2),
    ("stage_temperature_sum_2", 268.8, 403.2),
    ("stage_temperature_sum_3", 89.6, 134.4),
    ("stage_temperature_sum_4", 468.8, 703.2),
    ("vernalisation_requirement", 40.0, 60.0),
    ("daylength_requirement", 17.0, 23.0),
    ("max_assimilation_rate", 24.0, 36.0),
    ("specific_leaf_area_scale", 0.85, 1.15),
    ("root_penetration_rate", 0.0007, 0.0013),
    ("max_rooting_depth", 1.0, 1.6),
    ("drought_threshold_shift", -0.15, 0.10),
    ("sowing_shift_days", -14.0, 14.0),
    ("n_fertilizer_scale", 0.75, 1.25),
    ("irrigation_amount_mm", 0.0, 8.0),
    ("irrigation_threshold", 0.075, 0.14),
]


BASELINE = {
    "stage_temperature_sum_0": 146.0,
    "stage_temperature_sum_1": 276.0,
    "stage_temperature_sum_2": 336.0,
    "stage_temperature_sum_3": 112.0,
    "stage_temperature_sum_4": 586.0,
    "vernalisation_requirement": 50.0,
    "daylength_requirement": 20.0,
    "max_assimilation_rate": 30.0,
    "specific_leaf_area_scale": 1.0,
    "root_penetration_rate": 0.0011,
    "max_rooting_depth": 1.3,
    "drought_threshold_shift": 0.0,
    "sowing_shift_days": 0.0,
    "n_fertilizer_scale": 1.0,
    "irrigation_amount_mm": 0.0,
    "irrigation_threshold": 0.35,
}


def geometry_contains(shape: shapefile.Shape, points: np.ndarray) -> np.ndarray:
    geometry = shape.__geo_interface__
    polygons = [geometry["coordinates"]] if geometry["type"] == "Polygon" else geometry["coordinates"]
    inside = np.zeros(len(points), dtype=bool)
    for rings in polygons:
        in_polygon = PolygonPath(np.asarray(rings[0])).contains_points(points, radius=1e-9)
        for hole in rings[1:]:
            in_polygon &= ~PolygonPath(np.asarray(hole)).contains_points(points, radius=1e-9)
        inside |= in_polygon
    return inside


def atomic_write_rows(path: str | Path, fieldnames: list[str], rows: list[dict]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, output)


def create_points(config: dict[str, str]) -> list[dict]:
    soil_metadata, _ = Mrunlib.read_header(config["soil-grid"])
    crop_metadata, _ = Mrunlib.read_header(config["crop-grid"])
    for key in ("ncols", "nrows", "xllcorner", "yllcorner", "cellsize"):
        if not np.isclose(float(soil_metadata[key]), float(crop_metadata[key])):
            raise ValueError(f"Soil and crop grids differ at {key}")

    soil = np.loadtxt(config["soil-grid"], dtype=np.int32, skiprows=6)
    crop = np.loadtxt(config["crop-grid"], dtype=np.int16, skiprows=6)
    nodata = int(soil_metadata["nodata_value"])
    valid_rows, valid_cols = np.nonzero((soil != nodata) & (crop == 1))

    nrows = int(soil_metadata["nrows"])
    cellsize = float(soil_metadata["cellsize"])
    xll = float(soil_metadata["xllcorner"])
    yll = float(soil_metadata["yllcorner"])
    xs = xll + (valid_cols + 0.5) * cellsize
    ys = yll + (nrows - valid_rows - 0.5) * cellsize
    candidate_points = np.column_stack([xs, ys])

    shape_path = Path(config["shape-file"])
    prj_text = shape_path.with_suffix(".prj").read_text(encoding="utf-8", errors="ignore")
    if "UTM_Zone_32N" not in prj_text and "25832" not in prj_text:
        raise ValueError("SCR shapefile is not recognised as ETRS89 / UTM zone 32N")

    to_wgs84 = Transformer.from_crs(CRS.from_epsg(25832), CRS.from_epsg(4326), always_xy=True)
    reader = shapefile.Reader(str(shape_path), encoding="utf-8")
    output_rows = []

    for shape_record in reader.iterShapeRecords():
        record = shape_record.record.as_dict()
        bkr_id = int(float(record["BKR10_ID"]))
        min_x, min_y, max_x, max_y = shape_record.shape.bbox
        bbox_mask = (xs >= min_x) & (xs <= max_x) & (ys >= min_y) & (ys <= max_y)
        candidate_indices = np.flatnonzero(bbox_mask)
        local_inside = geometry_contains(shape_record.shape, candidate_points[candidate_indices])
        region_indices = candidate_indices[local_inside]
        if len(region_indices) == 0:
            raise ValueError(f"BKR {bkr_id} has no valid soil/crop candidate cells")

        region_soils = soil[valid_rows[region_indices], valid_cols[region_indices]]
        soil_ids, counts = np.unique(region_soils, return_counts=True)
        modal_soil_id = int(soil_ids[np.argmax(counts)])
        modal_indices = region_indices[region_soils == modal_soil_id]

        centroid_x = float(xs[region_indices].mean())
        centroid_y = float(ys[region_indices].mean())
        distance2 = (xs[modal_indices] - centroid_x) ** 2 + (ys[modal_indices] - centroid_y) ** 2
        selected = int(modal_indices[np.argmin(distance2)])
        lon, lat = to_wgs84.transform(float(xs[selected]), float(ys[selected]))

        output_rows.append(
            {
                "bkr_id": bkr_id,
                "name": record["NAME"],
                "srow": int(valid_rows[selected]),
                "scol": int(valid_cols[selected]),
                "soil_id": modal_soil_id,
                "x": f"{xs[selected]:.3f}",
                "y": f"{ys[selected]:.3f}",
                "lat": f"{lat:.8f}",
                "lon": f"{lon:.8f}",
                "candidate_cells": int(len(region_indices)),
                "modal_soil_share": f"{counts.max() / len(region_indices):.6f}",
            }
        )

    output_rows.sort(key=lambda row: row["bkr_id"])
    if len(output_rows) != 50:
        raise ValueError(f"Expected 50 SCR regions, found {len(output_rows)}")
    atomic_write_rows(config["points-file"], list(output_rows[0]), output_rows)
    return output_rows


def create_parameter_groups(config: dict[str, str]) -> list[dict]:
    total_groups = int(config["parameter-groups"])
    include_baseline = as_bool(config["include-baseline"])
    sampled_groups = total_groups - int(include_baseline)
    if sampled_groups < 1:
        raise ValueError("parameter-groups must leave room for at least one LHS group")

    sampler = qmc.LatinHypercube(d=len(PARAMETER_SPECS), seed=int(config["lhs-seed"]), optimization="random-cd")
    sample = sampler.random(sampled_groups)
    scaled = qmc.scale(
        sample,
        [spec[1] for spec in PARAMETER_SPECS],
        [spec[2] for spec in PARAMETER_SPECS],
    )

    rows = []
    if include_baseline:
        rows.append({"param_id": "P001", "irrigation_enabled": False, **BASELINE})

    for values in scaled:
        group_number = len(rows) + 1
        row = {"param_id": f"P{group_number:03d}", "irrigation_enabled": False}
        for (name, _, _), value in zip(PARAMETER_SPECS, values):
            row[name] = int(round(value)) if name == "sowing_shift_days" else float(value)
        row["irrigation_enabled"] = row["irrigation_amount_mm"] >= 4.0
        if not row["irrigation_enabled"]:
            row["irrigation_amount_mm"] = 0.0
        rows.append(row)

    fieldnames = ["param_id", "irrigation_enabled"] + [spec[0] for spec in PARAMETER_SPECS]
    atomic_write_rows(config["params-file"], fieldnames, rows)
    return rows


def main() -> None:
    config = parse_key_value_args(DEFAULTS, sys.argv[1:])
    points = create_points(config)
    groups = create_parameter_groups(config)
    print(f"Wrote {len(points)} SCR representative points to {config['points-file']}")
    print(f"Wrote {len(groups)} parameter groups to {config['params-file']}")
    print(f"MONICA environments: {len(points) * len(groups)}")
    print(f"Expected crop-season samples: {len(points) * len(groups) * 40}")
    irrigated = sum(bool(group["irrigation_enabled"]) for group in groups)
    print(f"Management groups: rainfed={len(groups) - irrigated}, irrigated={irrigated}")


if __name__ == "__main__":
    main()
