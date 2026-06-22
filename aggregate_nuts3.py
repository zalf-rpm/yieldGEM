#!/usr/bin/env python

import argparse
import csv
import re
from pathlib import Path

import numpy as np
import pandas as pd
import shapefile
from matplotlib.path import Path as MplPath
from pyproj import CRS, Transformer


COL_FILE_RE = re.compile(r"^col-(\d+)\.csv$", re.IGNORECASE)
DEFAULT_EXCLUDED_COLUMNS = {"Date", "Crop", "CM-count"}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate MONICA row/col daily CSV files to NUTS3 daily means."
    )
    parser.add_argument("--input-dir", required=True, help="Setup output directory containing row folders.")
    parser.add_argument("--nuts-shp", required=True, help="Path to the NUTS shapefile.")
    parser.add_argument("--grid-asc", required=True, help="Original ESRI ASCII grid used by MONICA.")
    parser.add_argument("--output", required=True, help="Output NUTS3 daily mean CSV.")
    parser.add_argument("--id-field", default="NUTS_ID", help="NUTS identifier field.")
    parser.add_argument("--level-field", default="LEVL_CODE", help="Optional NUTS level field.")
    parser.add_argument("--grid-crs", default="EPSG:25832", help="CRS of the MONICA grid.")
    parser.add_argument(
        "--nuts-crs",
        help="NUTS CRS if the .prj file is missing, for example EPSG:3035.",
    )
    parser.add_argument(
        "--variables",
        nargs="+",
        help="Variables to average. By default all numeric columns except CM-count are used.",
    )
    parser.add_argument(
        "--mapping-output",
        help="Optional CSV path for the Srow/Scol to NUTS3 mapping.",
    )
    return parser.parse_args()


def read_ascii_header(path):
    metadata = {}
    with open(path, encoding="utf-8") as handle:
        for _ in range(6):
            parts = handle.readline().split()
            if len(parts) >= 2:
                metadata[parts[0].lower()] = float(parts[1])

    required = {"nrows", "ncols", "xllcorner", "yllcorner", "cellsize"}
    missing = required.difference(metadata)
    if missing:
        raise ValueError(f"Missing ASC header fields: {sorted(missing)}")
    return metadata


def discover_cell_files(input_dir):
    cells = []
    for path in Path(input_dir).rglob("col-*.csv"):
        match = COL_FILE_RE.match(path.name)
        if not match:
            continue
        try:
            srow = int(path.parent.name)
        except ValueError:
            continue
        cells.append({"Srow": srow, "Scol": int(match.group(1)), "Path": str(path)})

    if not cells:
        raise FileNotFoundError(f"No row/col CSV files found below {input_dir}")
    return pd.DataFrame(cells).sort_values(["Srow", "Scol"]).reset_index(drop=True)


def add_cell_centers(cells, grid_meta):
    cellsize = grid_meta["cellsize"]
    nrows = int(grid_meta["nrows"])
    cells = cells.copy()
    cells["X"] = grid_meta["xllcorner"] + (cells["Scol"] + 0.5) * cellsize
    cells["Y"] = grid_meta["yllcorner"] + (nrows - cells["Srow"] - 0.5) * cellsize
    return cells


def read_nuts_crs(shp_path, explicit_crs):
    if explicit_crs:
        return CRS.from_user_input(explicit_crs)
    prj_path = Path(shp_path).with_suffix(".prj")
    if not prj_path.exists():
        raise FileNotFoundError("NUTS .prj file not found; provide --nuts-crs.")
    return CRS.from_wkt(prj_path.read_text(encoding="utf-8"))


def polygon_contains_points(polygon_coords, points):
    if not polygon_coords:
        return np.zeros(len(points), dtype=bool)
    inside = MplPath(np.asarray(polygon_coords[0])).contains_points(points, radius=1e-10)
    for hole in polygon_coords[1:]:
        inside &= ~MplPath(np.asarray(hole)).contains_points(points, radius=1e-10)
    return inside


def geometry_contains_points(geometry, points):
    if geometry["type"] == "Polygon":
        polygons = [geometry["coordinates"]]
    elif geometry["type"] == "MultiPolygon":
        polygons = geometry["coordinates"]
    else:
        return np.zeros(len(points), dtype=bool)

    inside = np.zeros(len(points), dtype=bool)
    for polygon in polygons:
        inside |= polygon_contains_points(polygon, points)
    return inside


def assign_nuts3(cells, shp_path, grid_crs, nuts_crs, id_field, level_field):
    transformer = Transformer.from_crs(grid_crs, nuts_crs, always_xy=True)
    x, y = transformer.transform(cells["X"].to_numpy(), cells["Y"].to_numpy())
    points = np.column_stack([x, y])
    assignment = np.full(len(cells), None, dtype=object)

    reader = shapefile.Reader(shp_path)
    field_names = [field[0] for field in reader.fields[1:]]
    if id_field not in field_names:
        raise KeyError(f"Field {id_field!r} not found. Available fields: {field_names}")
    id_index = field_names.index(id_field)
    level_index = field_names.index(level_field) if level_field in field_names else None

    for shape_record in reader.iterShapeRecords():
        record = shape_record.record
        if level_index is not None and int(record[level_index]) != 3:
            continue

        xmin, ymin, xmax, ymax = shape_record.shape.bbox
        candidates = np.flatnonzero(
            (assignment == None)  # noqa: E711
            & (points[:, 0] >= xmin)
            & (points[:, 0] <= xmax)
            & (points[:, 1] >= ymin)
            & (points[:, 1] <= ymax)
        )
        if len(candidates) == 0:
            continue

        geometry = shape_record.shape.__geo_interface__
        inside = geometry_contains_points(geometry, points[candidates])
        assignment[candidates[inside]] = str(record[id_index])

    mapped = cells.copy()
    mapped["NUTS_ID"] = assignment
    return mapped


def read_daily_section(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        for row in reader:
            if row and row[0].strip().strip('"') == "daily":
                header = next(reader, None)
                next(reader, None)  # units row
                if not header:
                    return pd.DataFrame()
                rows = []
                for values in reader:
                    if not values or not any(value.strip() for value in values):
                        break
                    rows.append(values)
                return pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header)
    return pd.DataFrame()


def discover_variables(sample_file, requested):
    frame = read_daily_section(sample_file)
    if frame.empty:
        raise ValueError(f"No daily data found in {sample_file}")
    if requested:
        missing = [name for name in requested if name not in frame.columns]
        if missing:
            raise KeyError(f"Requested variables missing from CSV: {missing}")
        return requested

    variables = []
    for column in frame.columns:
        if column in DEFAULT_EXCLUDED_COLUMNS:
            continue
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if numeric.notna().any():
            variables.append(column)
    if not variables:
        raise ValueError("No numeric daily variables found.")
    return variables


def aggregate_region(files, variables):
    total_sum = None
    total_count = None

    for path in files:
        frame = read_daily_section(path)
        if frame.empty or "Date" not in frame.columns:
            continue

        available = [name for name in variables if name in frame.columns]
        numeric = frame[available].apply(pd.to_numeric, errors="coerce").replace(-9999, np.nan)
        numeric = numeric.reindex(columns=variables)
        dates = frame["Date"].astype(str)
        cell_sum = numeric.fillna(0).groupby(dates).sum()
        cell_count = numeric.notna().astype(np.int64).groupby(dates).sum()

        if total_sum is None:
            total_sum = cell_sum
            total_count = cell_count
        else:
            total_sum = total_sum.add(cell_sum, fill_value=0)
            total_count = total_count.add(cell_count, fill_value=0)

    if total_sum is None:
        return pd.DataFrame(columns=["Date"] + variables)

    means = total_sum.divide(total_count.replace(0, np.nan)).sort_index()
    means.index.name = "Date"
    return means.reset_index()


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cells = discover_cell_files(input_dir)
    cells = add_cell_centers(cells, read_ascii_header(args.grid_asc))
    nuts_crs = read_nuts_crs(args.nuts_shp, args.nuts_crs)
    cells = assign_nuts3(
        cells,
        args.nuts_shp,
        CRS.from_user_input(args.grid_crs),
        nuts_crs,
        args.id_field,
        args.level_field,
    )

    if args.mapping_output:
        mapping_path = Path(args.mapping_output)
        mapping_path.parent.mkdir(parents=True, exist_ok=True)
        cells[["Srow", "Scol", "X", "Y", "NUTS_ID", "Path"]].to_csv(mapping_path, index=False)

    mapped = cells.dropna(subset=["NUTS_ID"])
    if mapped.empty:
        raise RuntimeError("No MONICA grid-cell centers fall inside the supplied NUTS3 polygons.")

    variables = discover_variables(mapped.iloc[0]["Path"], args.variables)
    print(f"Mapped {len(mapped)} of {len(cells)} cell files to NUTS3 regions.")
    print("Averaging variables:", ", ".join(variables))

    first_write = True
    for nuts_id, region_cells in mapped.groupby("NUTS_ID", sort=True):
        result = aggregate_region(region_cells["Path"], variables)
        if result.empty:
            continue
        result.insert(0, "NUTS_ID", nuts_id)
        result.to_csv(output_path, mode="w" if first_write else "a", header=first_write, index=False)
        first_write = False
        print(f"Wrote {nuts_id}: {len(region_cells)} cells, {len(result)} dates")

    if first_write:
        raise RuntimeError("No daily records were aggregated.")
    print(f"Finished: {output_path}")


if __name__ == "__main__":
    main()
