#!/usr/bin/env python3
"""Receive SCR ensemble results and atomically write one CSV per scenario."""

from __future__ import annotations

import csv
import json
import os
import sys
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import zmq

from scr_simulation_lib import (
    PIPELINE_ID,
    file_sha256,
    is_complete_output,
    limited,
    parse_key_value_args,
    read_parameter_groups,
    read_points,
    scenario_id,
    scenario_output_path,
)


DEFAULTS = {
    "server": "login01.cluster.zalf.de",
    "port": "7777",
    "shared_id": "scr-lhs",
    "points-file": "scr_points.csv",
    "params-file": "scr_parameter_groups.csv",
    "out": "scr-csv-out",
    "expected-seasons": "40",
    "max-points": "-1",
    "max-params": "-1",
    "timeout": "3600000",
    "progress-every": "25",
}


OUTPUT_FIELDS = [
    "ScenarioId",
    "BKR_ID",
    "ParameterId",
    "GrowingYear",
    "SeasonIndex",
    "ValidLength",
    "Date",
    "DOY",
    "DaysAfterSowing",
    "Srow",
    "Scol",
    "Crow",
    "Ccol",
    "SoilId",
    "CM-count",
    "Crop",
    "Tmin",
    "Tavg",
    "Tmax",
    "Precip",
    "Globrad",
    "Relhumid",
    "Wind",
    "Sunhours",
    "Irrig",
    "ET0",
    "Yield",
    "AbBiom",
    "LAI",
    "SowingDate",
    "SowingShiftDays",
    "NFertilizerScale",
    "IrrigationEnabled",
    "IrrigationAmountMm",
    "IrrigationThreshold",
]

STATUS_FIELDS = ["TimestampUTC", "ScenarioId", "Status", "Rows", "Seasons", "Message"]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def parse_result_date(value: Any) -> date:
    text = str(value).strip()
    if len(text) < 10:
        raise ValueError(f"Invalid MONICA date: {value}")
    return date.fromisoformat(text[:10])


def active_crop(value: Any) -> bool:
    if value is None:
        return False
    text = str(value).strip().lower()
    return text not in {"", "0", "none", "null", "nan"}


def extract_daily_rows(message: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for section in message.get("data", []):
        if str(section.get("origSpec", "")).strip('"') != "daily":
            continue
        for result in section.get("results", []):
            if "Date" in result:
                rows.append(dict(result))
    if not rows:
        raise ValueError("MONICA response contains no daily rows")
    rows.sort(key=lambda row: parse_result_date(row["Date"]))
    dates = [parse_result_date(row["Date"]) for row in rows]
    if len(set(dates)) != len(dates):
        raise ValueError("MONICA response contains duplicate daily dates")
    return rows


def split_crop_cycles(rows: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    cycles = []
    current = []
    previous_date = None
    for row in rows:
        current_date = parse_result_date(row["Date"])
        if active_crop(row.get("Crop")):
            if current and current_date != previous_date + timedelta(days=1):
                cycles.append(current)
                current = []
            current.append(row)
            previous_date = current_date
        elif current:
            cycles.append(current)
            current = []
            previous_date = None
    if current:
        cycles.append(current)
    return cycles


def select_complete_cycles(
    rows: list[dict[str, Any]], custom_id: dict[str, Any], expected_seasons: int
) -> list[list[dict[str, Any]]]:
    cycles = split_crop_cycles(rows)
    if len(cycles) < expected_seasons:
        raise ValueError(f"found only {len(cycles)} crop cycles, expected {expected_seasons}")

    selected = cycles[:expected_seasons]
    expected_years = list(
        range(int(custom_id["first_growing_year"]), int(custom_id["last_growing_year"]) + 1)
    )
    actual_years = [parse_result_date(cycle[-1]["Date"]).year for cycle in selected]
    if actual_years != expected_years:
        raise ValueError(f"growing years are {actual_years}, expected {expected_years}")

    if len(cycles) == expected_seasons:
        final_output_date = parse_result_date(rows[-1]["Date"])
        final_crop_date = parse_result_date(selected[-1][-1]["Date"])
        if final_crop_date >= final_output_date:
            raise ValueError("last expected crop is still active on the final simulation day")

    for index, cycle in enumerate(selected, start=1):
        dates = [parse_result_date(row["Date"]) for row in cycle]
        for previous, current in zip(dates, dates[1:]):
            if current != previous + timedelta(days=1):
                raise ValueError(f"season {index} contains a date gap after {previous}")
    return selected


def output_rows(
    cycles: list[list[dict[str, Any]]], custom_id: dict[str, Any]
) -> list[dict[str, Any]]:
    rows = []
    for season_index, cycle in enumerate(cycles, start=1):
        growing_year = parse_result_date(cycle[-1]["Date"]).year
        sowing_template = custom_id.get("sowing_date")
        if sowing_template:
            _, sowing_month, sowing_day = (int(part) for part in sowing_template.split("-"))
            sowing = date(growing_year - 1, sowing_month, sowing_day)
        else:
            sowing = parse_result_date(cycle[0]["Date"])
        valid_length = len(cycle)
        for result in cycle:
            current_date = parse_result_date(result["Date"])
            rows.append(
                {
                    "ScenarioId": custom_id["scenario_id"],
                    "BKR_ID": custom_id["bkr_id"],
                    "ParameterId": custom_id["param_id"],
                    "GrowingYear": growing_year,
                    "SeasonIndex": season_index,
                    "ValidLength": valid_length,
                    "Date": current_date.isoformat(),
                    "DOY": current_date.timetuple().tm_yday,
                    "DaysAfterSowing": (current_date - sowing).days,
                    "Srow": custom_id["srow"],
                    "Scol": custom_id["scol"],
                    "Crow": custom_id["crow"],
                    "Ccol": custom_id["ccol"],
                    "SoilId": custom_id["soil_id"],
                    "CM-count": result.get("CM-count", ""),
                    "Crop": result.get("Crop", ""),
                    "Tmin": result.get("Tmin", ""),
                    "Tavg": result.get("Tavg", ""),
                    "Tmax": result.get("Tmax", ""),
                    "Precip": result.get("Precip", ""),
                    "Globrad": result.get("Globrad", ""),
                    "Relhumid": result.get("Relhumid", ""),
                    "Wind": result.get("Wind", ""),
                    "Sunhours": result.get("Sunhours", ""),
                    "Irrig": result.get("Irrig", ""),
                    "ET0": result.get("ET0", ""),
                    "Yield": result.get("Yield", ""),
                    "AbBiom": result.get("AbBiom", ""),
                    "LAI": result.get("LAI", ""),
                    "SowingDate": sowing.isoformat(),
                    "SowingShiftDays": custom_id["sowing_shift_days"],
                    "NFertilizerScale": custom_id["n_fertilizer_scale"],
                    "IrrigationEnabled": custom_id["irrigation_enabled"],
                    "IrrigationAmountMm": custom_id["irrigation_amount_mm"],
                    "IrrigationThreshold": custom_id["irrigation_threshold"],
                }
            )
    return rows


def scenario_path_from_custom_id(output_dir: str | Path, custom_id: dict[str, Any]) -> Path:
    return (
        Path(output_dir)
        / f"BKR_{int(custom_id['bkr_id']):03d}"
        / f"{custom_id['scenario_id']}_daily.csv"
    )


def atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            temporary.unlink()
        except OSError:
            pass
        raise


def append_status(path: Path, row: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    needs_header = not path.is_file() or path.stat().st_size == 0
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATUS_FIELDS)
        if needs_header:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def write_scenario_result(
    message: dict[str, Any], output_dir: str | Path, expected_seasons: int
) -> tuple[str, int, int]:
    custom_id = message.get("customId", {})
    sid = str(custom_id.get("scenario_id", ""))
    if custom_id.get("pipeline") != PIPELINE_ID or not sid:
        raise ValueError("response does not belong to the SCR LHS pipeline")
    errors = message.get("errors", [])
    if errors:
        raise RuntimeError("; ".join(str(error) for error in errors))
    if int(custom_id.get("expected_seasons", -1)) != expected_seasons:
        raise ValueError("producer and consumer expected-seasons differ")

    daily_rows = extract_daily_rows(message)
    cycles = select_complete_cycles(daily_rows, custom_id, expected_seasons)
    flattened = output_rows(cycles, custom_id)
    output_path = scenario_path_from_custom_id(output_dir, custom_id)
    atomic_write_csv(output_path, OUTPUT_FIELDS, flattened)
    return sid, len(flattened), len(cycles)


def expected_scenarios(config: dict[str, Any]) -> tuple[set[str], dict[str, Path]]:
    points = limited(read_points(config["points-file"]), config["max-points"])
    parameters = limited(read_parameter_groups(config["params-file"]), config["max-params"])
    ids = set()
    paths = {}
    for point in points:
        for group in parameters:
            sid = scenario_id(point, group)
            if sid in ids:
                raise ValueError(f"Duplicate scenario id: {sid}")
            ids.add(sid)
            paths[sid] = scenario_output_path(config["out"], point, group)
    return ids, paths


def write_missing(path: Path, scenario_ids: set[str], reason: str) -> None:
    rows = [{"ScenarioId": sid, "Reason": reason} for sid in sorted(scenario_ids)]
    atomic_write_csv(path, ["ScenarioId", "Reason"], rows)


def validate_run_manifest(config: dict[str, Any], output_paths: dict[str, Path]) -> None:
    output_dir = Path(config["out"])
    manifest_path = output_dir / "run_manifest.json"
    completed = any(is_complete_output(path) for path in output_paths.values())
    if not manifest_path.is_file():
        if completed:
            raise RuntimeError(
                f"{output_dir} contains scenario CSV files but no run_manifest.json; "
                "use the matching producer or a new output directory"
            )
        return
    with manifest_path.open("r", encoding="utf-8") as handle:
        signature = json.load(handle).get("run_signature", {})
    expected = {
        "pipeline": PIPELINE_ID,
        "expected_seasons": int(config["expected-seasons"]),
        "points_sha256": file_sha256(config["points-file"]),
        "parameters_sha256": file_sha256(config["params-file"]),
    }
    differences = [key for key, value in expected.items() if signature.get(key) != value]
    if differences:
        raise RuntimeError(
            f"run_manifest.json does not match the consumer inputs ({', '.join(differences)})"
        )


def main() -> int:
    config = parse_key_value_args(DEFAULTS, sys.argv[1:])
    expected_ids, output_paths = expected_scenarios(config)
    validate_run_manifest(config, output_paths)
    already_complete = {sid for sid, path in output_paths.items() if is_complete_output(path)}
    awaiting = expected_ids - already_complete
    expected_seasons = int(config["expected-seasons"])
    status_path = Path(config["out"]) / "run_status.csv"
    missing_path = Path(config["out"]) / "missing_scenarios.csv"

    print(
        f"Consumer expects {len(expected_ids)} scenarios: "
        f"complete={len(already_complete)} awaiting={len(awaiting)}"
    )
    if not awaiting:
        write_missing(missing_path, set(), "")
        print("All scenario CSV files are already complete")
        return 0

    context = zmq.Context()
    if config["shared_id"]:
        socket = context.socket(zmq.DEALER)
        socket.setsockopt(zmq.IDENTITY, str(config["shared_id"]).encode("utf-8"))
    else:
        socket = context.socket(zmq.PULL)
    socket.setsockopt(zmq.RCVTIMEO, int(config["timeout"]))
    socket.setsockopt(zmq.LINGER, 0)
    socket.connect(f"tcp://{config['server']}:{config['port']}")

    failures = set()
    received = 0
    progress_every = max(1, int(config["progress-every"]))
    try:
        while awaiting:
            try:
                message = socket.recv_json()
            except zmq.error.Again:
                write_missing(missing_path, awaiting, "timeout waiting for MONICA response")
                print(f"Timed out with {len(awaiting)} scenarios still missing")
                return 2

            custom_id = message.get("customId", {})
            sid = str(custom_id.get("scenario_id", ""))
            if custom_id.get("pipeline") != PIPELINE_ID or sid not in expected_ids:
                print(f"Ignoring unrelated response customId={custom_id}")
                continue
            if sid not in awaiting:
                print(f"Ignoring duplicate response for {sid}")
                continue

            try:
                _, row_count, season_count = write_scenario_result(
                    message, config["out"], expected_seasons
                )
                status = "complete"
                detail = ""
            except Exception as exc:
                row_count = 0
                season_count = 0
                status = "error"
                detail = str(exc)
                failures.add(sid)
                print(f"Failed {sid}: {detail}")

            append_status(
                status_path,
                {
                    "TimestampUTC": utc_now(),
                    "ScenarioId": sid,
                    "Status": status,
                    "Rows": row_count,
                    "Seasons": season_count,
                    "Message": detail,
                },
            )
            awaiting.remove(sid)
            received += 1
            if received % progress_every == 0 or not awaiting:
                print(
                    f"Processed {received} responses; "
                    f"remaining={len(awaiting)} errors={len(failures)}"
                )
    except KeyboardInterrupt:
        write_missing(missing_path, awaiting, "consumer interrupted")
        print(f"Interrupted with {len(awaiting)} scenarios still missing")
        return 130
    finally:
        socket.close()
        context.term()

    write_missing(missing_path, failures, "MONICA or output validation error")
    if failures:
        print(f"Finished with {len(failures)} failed scenarios; rerun producer and consumer to retry")
        return 1
    print(f"Finished successfully: {len(expected_ids)} scenario CSV files are complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
