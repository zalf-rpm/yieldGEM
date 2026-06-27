#!/usr/bin/env python3
"""Shared helpers for the SCR parameter-ensemble MONICA pipeline."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from typing import Any, Iterable


PIPELINE_ID = "scr-lhs-v1"


def parse_key_value_args(defaults: dict[str, Any], args: Iterable[str]) -> dict[str, Any]:
    config = dict(defaults)
    for arg in args:
        if "=" not in arg:
            raise ValueError(f"Expected key=value argument, got: {arg}")
        key, value = arg.split("=", maxsplit=1)
        if key in config:
            config[key] = value
    return config


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", ""}:
        return False
    raise ValueError(f"Cannot interpret as boolean: {value}")


def limited(records: list[dict[str, Any]], maximum: Any) -> list[dict[str, Any]]:
    limit = int(maximum)
    return records if limit < 0 else records[:limit]


def read_csv_records(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_points(path: str | Path) -> list[dict[str, Any]]:
    points = []
    for row in read_csv_records(path):
        points.append(
            {
                **row,
                "bkr_id": int(float(row["bkr_id"])),
                "srow": int(row["srow"]),
                "scol": int(row["scol"]),
                "soil_id": int(row["soil_id"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "lat": float(row["lat"]),
                "lon": float(row["lon"]),
            }
        )
    if not points:
        raise ValueError(f"No SCR points found in {path}")
    return points


def read_parameter_groups(path: str | Path) -> list[dict[str, Any]]:
    groups = []
    for row in read_csv_records(path):
        parsed: dict[str, Any] = dict(row)
        for key, value in row.items():
            if key == "param_id":
                continue
            if key == "irrigation_enabled":
                parsed[key] = as_bool(value)
            else:
                parsed[key] = float(value)
        groups.append(parsed)
    if not groups:
        raise ValueError(f"No parameter groups found in {path}")
    return groups


def scenario_id(point: dict[str, Any], parameters: dict[str, Any]) -> str:
    return f"BKR{int(point['bkr_id']):03d}_{parameters['param_id']}"


def scenario_output_path(
    output_dir: str | Path,
    point: dict[str, Any],
    parameters: dict[str, Any],
) -> Path:
    sid = scenario_id(point, parameters)
    return Path(output_dir) / f"BKR_{int(point['bkr_id']):03d}" / f"{sid}_daily.csv"


def is_complete_output(path: str | Path) -> bool:
    output = Path(path)
    return output.is_file() and output.stat().st_size > 0


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
