#!/usr/bin/env python3
"""List the DWD climate files required by the SCR representative points."""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

from pyproj import CRS

import monica_run_lib as Mrunlib
from scr_simulation_lib import parse_key_value_args, read_points


DEFAULTS = {
    "points-file": "scr_points.csv",
    "mapping-file": "data/dwd/csvs/latlon_to_rowcol.json",
    "climate-root": "data/dwd/csvs/germany_ubn_1951-01-01_to_2024-08-30",
    "flat-climate-root": "",
    "output": "scr_required_climate_files.csv",
}


def scalar_int(value) -> int:
    try:
        return int(value.item())
    except AttributeError:
        return int(value)


def main() -> None:
    config = parse_key_value_args(DEFAULTS, sys.argv[1:])
    mapping_path = Path(config["mapping-file"])
    if not mapping_path.is_file():
        raise FileNotFoundError(
            f"Missing {mapping_path}. Copy latlon_to_rowcol.json from the HPC climate archive first."
        )

    climate_coordinates = {}
    interpolator = Mrunlib.create_climate_geoGrid_interpolator_from_json_file(
        str(mapping_path),
        CRS.from_epsg(4326),
        CRS.from_epsg(25832),
        climate_coordinates,
    )

    cell_to_bkr = defaultdict(list)
    for point in read_points(config["points-file"]):
        crow, ccol = interpolator(float(point["x"]), float(point["y"]))
        cell_to_bkr[(scalar_int(crow), scalar_int(ccol))].append(int(point["bkr_id"]))

    climate_root = Path(config["climate-root"])
    flat_climate_root = str(config["flat-climate-root"]).strip()
    rows = []
    for (crow, ccol), bkr_ids in sorted(cell_to_bkr.items()):
        relative_path = Path(str(crow)) / f"daily_mean_RES1_C{ccol}R{crow}.csv.gz"
        full_path = (
            Path(flat_climate_root) / relative_path.name
            if flat_climate_root
            else climate_root / relative_path
        )
        rows.append(
            {
                "crow": crow,
                "ccol": ccol,
                "relative_path": relative_path.as_posix(),
                "full_path": str(full_path),
                "bkr_ids": ",".join(str(value) for value in sorted(bkr_ids)),
                "bkr_count": len(bkr_ids),
                "file_exists": full_path.is_file(),
            }
        )

    output_path = Path(config["output"])
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    missing = sum(not row["file_exists"] for row in rows)
    print(f"SCR points: {sum(len(ids) for ids in cell_to_bkr.values())}")
    print(f"Unique climate files: {len(rows)}")
    print(f"Missing climate files: {missing}")
    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
