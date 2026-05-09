#!/usr/bin/python
# -*- coding: UTF-8

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/. */

# Authors:
# Michael Berg-Mohnicke <michael.berg@zalf.de>
#
# Maintainers:
# Currently maintained by the authors.
#
# This file has been created at the Institute of
# Landscape Systems Analysis at the ZALF.
# Copyright (C: Leibniz Centre for Agricultural Landscape Research (ZALF)

from collections import defaultdict, OrderedDict
import csv
import numpy as np
import os
from pyproj import CRS, Transformer
import sys
import zmq

import monica_io3
import monica_run_lib as Mrunlib

PATHS = {
    "re-local-remote": {
        "path-to-data-dir": "data/",
        "path-to-output-dir": "D:/projects/KlimErtrag/out_remote_local/",
        "path-to-csv-output-dir": "D:/projects/KlimErtrag/out_remote_local/"
    },
    "mbm-local-remote": {
        "path-to-data-dir": "data/",
        "path-to-output-dir": "out/",
        "path-to-csv-output-dir": "csv-out/"
    },
    "remoteConsumer-remoteMonica": {
        "path-to-data-dir": "./data/",
        "path-to-output-dir": "/out/out/",
        "path-to-csv-output-dir": "/out/csv-out/"
    }
}

TEMPLATE_SOIL_PATH = "{local_path_to_data_dir}germany/buek200_1000_25832_etrs89-utm32n.asc"
# TEMPLATE_LANDUSE_PATH = "{local_path_to_data_dir}germany/landuse_1000_31469_gk5.asc"
DATA_SOIL_DB = "germany/buek200.sqlite"
USE_LANDUSE = False


def create_output(msg):
    cm_count_to_vals = defaultdict(dict)
    for data in msg.get("data", []):
        results = data.get("results", [])

        is_daily_section = data.get("origSpec", "") == '"daily"'

        for vals in results:
            if "CM-count" in vals:
                cm_count_to_vals[vals["CM-count"]].update(vals)
            elif is_daily_section and "Date" in vals:
                cm_count_to_vals[vals["Date"]].update(vals)

    if not cm_count_to_vals:
        return cm_count_to_vals

    cmcs = sorted(cm_count_to_vals.keys())
    last_cmc = cmcs[-1]
    if "Year" not in cm_count_to_vals[last_cmc]:
        cm_count_to_vals.pop(last_cmc, None)

    return cm_count_to_vals


def write_row_to_grids(row_col_data, row, ncols, header, path_to_output_dir, path_to_csv_output_dir, setup_id):
    """write grids row by row"""

    if not hasattr(write_row_to_grids, "list_of_output_files"):
        write_row_to_grids.list_of_output_files = defaultdict(list)

    if not hasattr(write_row_to_grids, "cmc_to_crop"):
        write_row_to_grids.cmc_to_crop = defaultdict(dict)

    if not hasattr(write_row_to_grids, "file_rows_written"):
        write_row_to_grids.file_rows_written = defaultdict(int)

    cmc_to_crop = write_row_to_grids.cmc_to_crop[setup_id]

    make_dict_nparr = lambda: defaultdict(lambda: np.full((ncols,), -9999, dtype=float))

    output_grids = {
        "Yield": {"data": make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "AbBiom": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "LAI": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "EffRootDep": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # # "NFert": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "N": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "N_30": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "N_90": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "Nstress": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "Nmin": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "iniNmin": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "finalNmin": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},         
        # "SumNUp": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "AbBiomNc_last": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "AbBiomNc_max": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "NLeach": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
            
        # "PASW": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "PASW_30": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "PASW_90": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Mois": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Mois_30": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Mois_90": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "WaterContent": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # "Act_ET": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Evaporation": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Tra": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "TraDef": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # "Precip": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "Drain": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},            

        # # "Nstress_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # # "SumNUp_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # "PASW_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "PASW_30_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "PASW_90_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "WaterContent_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # "Act_ET_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Evaporation_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "Tra_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},
        # "TraDef_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1},

        # "Precip_spring": {"data" : make_dict_nparr(), "cast-to": "float", "digits": 1}
    }
    output_keys = list(output_grids.keys())

    def ensure_file_ready(path_to_file: str, current_row: int):
        if not os.path.isfile(path_to_file):
            with open(path_to_file, "w") as f:
                f.write(header)
            write_row_to_grids.list_of_output_files[setup_id].append(path_to_file)
            write_row_to_grids.file_rows_written[path_to_file] = 0

        already = write_row_to_grids.file_rows_written[path_to_file]
        missing = current_row - already
        if missing > 0:
            nodata_line = " ".join(["-9999"] * ncols) + "\n"
            with open(path_to_file, "a") as f:
                for _ in range(missing):
                    f.write(nodata_line)
            write_row_to_grids.file_rows_written[path_to_file] = current_row

    # skip this part if we write just a nodata line
    if row in row_col_data:
        for col in range(0, ncols):
            if col not in row_col_data[row]:
                continue

            rcd_val = row_col_data[row][col]

            if rcd_val == -9999:
                continue

            cmc_and_year_to_vals = defaultdict(lambda: defaultdict(list))

            for cell_data in rcd_val:
                # if we got multiple datasets per cell, iterate over them and aggregate them in the following step
                for cm_count, data in cell_data.items():
                    if "Crop" in data:
                        c = str(data["Crop"]).strip()
                        if c:
                            cmc_to_crop[cm_count] = c

                    year = data.get("Year", None)
                    if year is None:
                        continue

                    for key in output_keys:
                        # only further process/store data we actually received
                        if key in data:
                            v = data[key]
                            if isinstance(v, list):
                                for i, v_ in enumerate(v):
                                    cmc_and_year_to_vals[(cm_count, year)][f"{key}_{i + 1}"].append(v_)
                            else:
                                cmc_and_year_to_vals[(cm_count, year)][key].append(v)
                        # if a key is missing, because that monica event was never raised/reached, create the empty list
                        # so a no-data value is being produced
                        else:
                            cmc_and_year_to_vals[(cm_count, year)][key]

            # potentially aggregate multiple data per cell and finally store them for this row
            for (cm_count, year), key_to_vals in cmc_and_year_to_vals.items():
                for key, vals in key_to_vals.items():
                    if key not in output_grids:
                        continue
                    out = output_grids[key]["data"]
                    out[(cm_count, year)][col] = (sum(vals) / len(vals)) if vals else -9999

    # iterate over all prepared data for a single row and write row
    for key, y2d_ in output_grids.items():
        y2d = y2d_["data"]
        digits = y2d_.get("digits", 0)

        mold = (lambda x: str(round(float(x), digits)))

        for (cm_count, year), row_arr in y2d.items():
            crop = str(cmc_to_crop.get(cm_count, "none")).strip() or "none"
            crop = crop.replace("/", "").replace(" ", "")
            key2 = key.replace("/", "_")
            path_to_file = f"{path_to_output_dir}{crop}_{key2}_{year}_{cm_count}.asc"

            ensure_file_ready(path_to_file, row)

            rowstr = " ".join(["-9999" if int(x) == -9999 else mold(x) for x in row_arr])
            with open(path_to_file, "a") as f:
                f.write(rowstr + "\n")

            write_row_to_grids.file_rows_written[path_to_file] += 1

    if row in row_col_data:
        del row_col_data[row]


def finalize_outputs(setup_id: int, total_rows: int, ncols: int):
    if not hasattr(write_row_to_grids, "list_of_output_files"):
        return
    if not hasattr(write_row_to_grids, "file_rows_written"):
        return

    nodata_line = " ".join(["-9999"] * ncols) + "\n"

    for path_to_file in write_row_to_grids.list_of_output_files.get(setup_id, []):
        already = write_row_to_grids.file_rows_written.get(path_to_file, 0)
        missing = total_rows - already
        if missing > 0:
            with open(path_to_file, "a") as f:
                for _ in range(missing):
                    f.write(nodata_line)
            write_row_to_grids.file_rows_written[path_to_file] = total_rows


def write_daily_csv(daily_data_dict, path_to_csv_output_dir):
    """write daily yields data to CSV file"""

    print(f"[DEBUG] write_daily_csv called with data entries: {len(daily_data_dict)}")

    if not daily_data_dict:
        print("[DEBUG] No daily data to write")
        return

    # Create output directory if needed
    if not os.path.exists(path_to_csv_output_dir):
        try:
            os.makedirs(path_to_csv_output_dir)
            print(f"[DEBUG] Created directory: {path_to_csv_output_dir}")
        except OSError:
            print("c: Couldn't create dir:", path_to_csv_output_dir, "! Exiting.")
            return

    # daily_data_dict structure: {(crop, cm_count): {date: data_dict, ...}}
    for (crop, cm_count), date_to_data in daily_data_dict.items():
        print(f"[DEBUG] Processing crop={crop}, cm_count={cm_count}, dates={len(date_to_data)}")
        file_path = f"{path_to_csv_output_dir}{crop}_daily_yields_{cm_count}.csv"

        # Sort by date and get all unique field names
        sorted_dates = sorted(date_to_data.keys())
        if not sorted_dates:
            continue

        # Collect all field names
        all_fields = set()
        for data in date_to_data.values():
            all_fields.update(data.keys())

        fieldnames = ["Date"] + sorted([f for f in all_fields if f != "Date" and f != "CM-count"])

        try:
            with open(file_path, "w", newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                for date_str in sorted_dates:
                    row = {"Date": date_str}
                    row.update({k: v for k, v in date_to_data[date_str].items() if k in fieldnames})
                    writer.writerow(row)

            print(f"Wrote daily CSV: {file_path}")
        except Exception as e:
            print(f"Error writing daily CSV {file_path}: {e}")


def run_consumer(leave_after_finished_run=True, server={"server": None, "port": None}, shared_id=None):
    """collect data from workers"""

    config = {
        "mode": "re-local-remote",  # "mbm-local-remote",
        "port": server["port"] if server["port"] else "7778",
        "server": server["server"] if server["server"] else "login01.cluster.zalf.de",
        "start-row": "0",
        "end-row": "-1",
        "shared_id": shared_id,
        "timeout": 600000  # 10 minutes
    }

    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=")
            if k in config:
                config[k] = v

    paths = PATHS[config["mode"]]

    if not "out" in config:
        config["out"] = paths["path-to-output-dir"]
    if not "csv-out" in config:
        config["csv-out"] = paths["path-to-csv-output-dir"]

    print("consumer config:", config)

    context = zmq.Context()
    if config["shared_id"]:
        socket = context.socket(zmq.DEALER)
        socket.setsockopt(zmq.IDENTITY, config["shared_id"])
    else:
        socket = context.socket(zmq.PULL)

    socket.connect("tcp://" + config["server"] + ":" + config["port"])
    socket.RCVTIMEO = config["timeout"]
    leave = False
    write_normal_output_files = False

    path_to_soil_grid = TEMPLATE_SOIL_PATH.format(local_path_to_data_dir=paths["path-to-data-dir"])
    soil_epsg_code = int(path_to_soil_grid.split("/")[-1].split("_")[2])
    soil_crs = CRS.from_epsg(soil_epsg_code)
    soil_metadata, header = Mrunlib.read_header(path_to_soil_grid)
    soil_grid_template = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)

    scols = int(soil_metadata["ncols"])
    srows = int(soil_metadata["nrows"])
    scellsize = int(soil_metadata["cellsize"])
    xllcorner = int(soil_metadata["xllcorner"])
    yllcorner = int(soil_metadata["yllcorner"])
    nodata_value = int(soil_metadata["nodata_value"])

    if USE_LANDUSE:
        path_to_landuse_grid = TEMPLATE_LANDUSE_PATH.format(local_path_to_data_dir=paths["path-to-data-dir"])
        landuse_epsg_code = int(path_to_landuse_grid.split("/")[-1].split("_")[2])
        landuse_crs = CRS.from_epsg(landuse_epsg_code)
        landuse_transformer = Transformer.from_crs(soil_crs, landuse_crs)
        landuse_meta, _ = Mrunlib.read_header(path_to_landuse_grid)
        landuse_grid = np.loadtxt(path_to_landuse_grid, dtype=int, skiprows=6)
        landuse_interpolate = Mrunlib.create_ascii_grid_interpolator(landuse_grid, landuse_meta)

        for srow in range(0, srows):
            # print(srow)
            for scol in range(0, scols):
                soil_id = soil_grid_template[srow, scol]
                if soil_id == -9999:
                    continue

                # get coordinate of clostest climate element of real soil-cell
                sh = yllcorner + (scellsize / 2) + (srows - srow - 1) * scellsize
                sr = xllcorner + (scellsize / 2) + scol * scellsize

                # check if current grid cell is used for agriculture                
                lur, luh = landuse_transformer(sh, sr)
                landuse_id = landuse_interpolate(lur, luh)
                if landuse_id not in [2, 3, 4]:
                    soil_grid_template[srow, scol] = -9999

        print("filtered through CORINE")

    # set all data values to one, to count them later
    soil_grid_template[soil_grid_template != nodata_value] = 1
    # set all no-data values to 0, to ignore them while counting
    soil_grid_template[soil_grid_template == nodata_value] = 0

    # count cols in rows
    datacells_per_row = np.sum(soil_grid_template, axis=1)

    start_row = int(config["start-row"])
    end_row = int(config["end-row"])
    ncols = int(soil_metadata["ncols"])
    setup_id_to_data = defaultdict(lambda: {
        "start_row": start_row,
        "end_row": end_row,
        "nrows": end_row - start_row + 1 if start_row > 0 and end_row >= start_row else int(soil_metadata["nrows"]),
        "ncols": ncols,
        "header": header,
        "out_dir_exists": False,
        "row-col-data": defaultdict(lambda: defaultdict(list)),
        "daily-data": defaultdict(lambda: defaultdict(dict)),
        "datacell-count": datacells_per_row.copy(),
        "next-row": start_row
    })

    def process_message(msg):
        if len(msg["errors"]) > 0:
            cid = msg.get("customId", {})
            print("There were errors for", cid, "errors:", msg["errors"])

            setup_id = cid.get("setup_id")
            row = cid.get("srow")
            col = cid.get("scol")
            if setup_id is not None and row is not None and col is not None:
                data = setup_id_to_data[setup_id]
                data["row-col-data"][row][col] = -9999
                data["datacell-count"][row] -= 1
            return False

        if not hasattr(process_message, "wnof_count"):
            process_message.wnof_count = 0
            process_message.setup_count = 0

        leave = False

        if not write_normal_output_files:
            custom_id = msg["customId"]
            setup_id = custom_id["setup_id"]
            is_nodata = custom_id["nodata"]

            data = setup_id_to_data[setup_id]

            row = custom_id["srow"]
            col = custom_id["scol"]
            # crow = custom_id.get("crow", -1)
            # ccol = custom_id.get("ccol", -1)
            # soil_id = custom_id.get("soil_id", -1)

            debug_msg = "received work result " + str(process_message.received_env_count) + " customId: " + str(
                msg.get("customId", "")) \
                        + " next row: " + str(data["next-row"]) \
                        + " cols@row to go: " + str(data["datacell-count"][row]) + "@" + str(
                row) + " cells_per_row: " + str(datacells_per_row[row])  # \
            # + " rows unwritten: " + str(data["row-col-data"].keys())
            # print(debug_msg)
            # debug_file.write(debug_msg + "\n")
            if is_nodata:
                data["row-col-data"][row][col] = -9999
            else:
                output = create_output(msg)
                # Check if this is daily data by looking at message origSpec
                is_daily = any(d.get("origSpec", "") == '"daily"' for d in msg.get("data", []))

                origSpecs = [d.get("origSpec", "") for d in msg.get("data", [])]
                print(f"[DEBUG] origSpecs: {origSpecs}, is_daily: {is_daily}")

                if is_daily:
                    # Store daily data separately (key: date, value: output dict)
                    print(f"[DEBUG] Processing daily data, output keys: {list(output.keys())}")
                    for date_key, date_data in output.items():
                        print(f"[DEBUG] date_data type: {type(date_data)}, content: {date_data}")
                        crop = str(date_data.get("Crop", "unknown")).strip()
                        crop = crop.replace("/", "").replace(" ", "") or "unknown"
                        cm_count = date_data.get("CM-count")
                        print(f"[DEBUG] Date: {date_key}, Crop: {crop}, CM-count: {cm_count}, all_keys: {list(date_data.keys())}")
                        if cm_count:
                            key = (crop, cm_count)
                            data["daily-data"][key][date_key] = date_data
                            print(f"[DEBUG] Stored daily data for key: {key}")
                else:
                    # Store regular (event-based) data
                    data["row-col-data"][row][col].append(output)
            data["datacell-count"][row] -= 1

            process_message.received_env_count = process_message.received_env_count + 1

            while (data["next-row"] in data["row-col-data"] and data["datacell-count"][data["next-row"]] == 0) \
                    or (
                    len(data["datacell-count"]) > data["next-row"] and data["datacell-count"][data["next-row"]] == 0):

                path_to_out_dir = config["out"] + str(setup_id) + "/"
                path_to_csv_out_dir = config["csv-out"] + str(setup_id) + "/"
                print(path_to_out_dir)
                if not data["out_dir_exists"]:
                    if os.path.isdir(path_to_out_dir) and os.path.exists(path_to_out_dir):
                        data["out_dir_exists"] = True
                    else:
                        try:
                            os.makedirs(path_to_out_dir)
                            data["out_dir_exists"] = True
                        except OSError:
                            print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                            exit(1)
                    if os.path.isdir(path_to_csv_out_dir) and os.path.exists(path_to_csv_out_dir):
                        data["out_dir_exists"] = True
                    else:
                        try:
                            os.makedirs(path_to_csv_out_dir)
                            data["out_dir_exists"] = True
                        except OSError:
                            print("c: Couldn't create dir:", path_to_csv_out_dir, "! Exiting.")
                            exit(1)

                # write_row_to_grids(data["row-col-data"], data["next-row"], data["ncols"], data["header"],
                #                    path_to_out_dir, path_to_csv_out_dir, setup_id)

                # debug_msg = "wrote row: " + str(data["next-row"]) + " next-row: " + str(
                #     data["next-row"] + 1) + " rows unwritten: " + str(list(data["row-col-data"].keys()))
                # print(debug_msg)
                # debug_file.write(debug_msg + "\n")

                data["next-row"] += 1  # move to next row (to be written)

                if leave_after_finished_run \
                        and ((data["end_row"] < 0 and data["next-row"] > data["nrows"] - 1)
                             or (0 <= data["end_row"] < data["next-row"])):
                    process_message.setup_count += 1

        elif write_normal_output_files:
            if msg.get("type", "") in ["jobs-per-cell", "no-data", "setup_data"]:
                # print "ignoring", result.get("type", "")
                return

            print("received work result ", process_message.received_env_count, " customId: ",
                  str(msg.get("customId", "")))

            custom_id = msg["customId"]
            is_nodata = custom_id["nodata"]
            if is_nodata:
                return leave
            setup_id = custom_id["setup_id"]
            row = custom_id["srow"]
            col = custom_id["scol"]
            # crow = custom_id.get("crow", -1)
            # ccol = custom_id.get("ccol", -1)
            # soil_id = custom_id.get("soil_id", -1)

            process_message.wnof_count += 1

            path_to_out_dir = config["out"] + str(setup_id) + "/" + str(row) + "/"
            print(path_to_out_dir)
            if not os.path.exists(path_to_out_dir):
                try:
                    os.makedirs(path_to_out_dir)
                except OSError:
                    print("c: Couldn't create dir:", path_to_out_dir, "! Exiting.")
                    exit(1)

            # with open("out/out-" + str(i) + ".csv", 'wb') as _:
            with open(path_to_out_dir + "col-" + str(col) + ".csv", "w", newline='') as _:

                writer = csv.writer(_, delimiter=",")
                for data_ in msg.get("data", []):
                    results = data_.get("results", [])
                    orig_spec = data_.get("origSpec", "")
                    output_ids = data_.get("outputIds", [])

                    if len(results) > 0:
                        writer.writerow([orig_spec.replace("\"", "")])
                        for row in monica_io3.write_output_header_rows(output_ids,
                                                                       include_header_row=True,
                                                                       include_units_row=True,
                                                                       include_time_agg=False):
                            writer.writerow(row)

                        # for row in monica_io3.write_output(output_ids, results):
                        #     writer.writerow(row)

                        for result in results:
                            row = []
                            for output_id in output_ids:
                                field_name = output_id["name"]
                                row.append(result.get(field_name, ""))
                            writer.writerow(row)

                writer.writerow([])

            process_message.received_env_count = process_message.received_env_count + 1

        return leave

    process_message.received_env_count = 1

    while not leave:
        try:
            # start_time_recv = timeit.default_timer()
            msg = socket.recv_json()  # encoding="latin-1"
            # elapsed = timeit.default_timer() - start_time_recv
            # print("time to receive message" + str(elapsed))
            # start_time_proc = timeit.default_timer()
            print(f"[DEBUG] Received message #{process_message.received_env_count}")
            leave = process_message(msg)
            # elapsed = timeit.default_timer() - start_time_proc
            # print("time to process message" + str(elapsed))
        except zmq.error.Again as _e:
            print('no response from the server (with "timeout"=%d ms) ' % socket.RCVTIMEO)
            for setup_id, data in setup_id_to_data.items():
                print(f"[DEBUG] Setup {setup_id} - daily-data keys: {list(data.get('daily-data', {}).keys())}")
                print(f"[DEBUG] Setup {setup_id} - daily-data total entries: {sum(len(v) for v in data.get('daily-data', {}).values())}")
                # finalize_outputs(setup_id, data["nrows"], data["ncols"])
                # Write daily CSV data if available
                if data.get("daily-data"):
                    path_to_csv_out_dir = config["csv-out"] + str(setup_id) + "/"
                    print(f"[DEBUG] Writing daily CSV to: {path_to_csv_out_dir}")
                    write_daily_csv(dict(data["daily-data"]), path_to_csv_out_dir)
                else:
                    print(f"[DEBUG] No daily-data found for setup {setup_id}")
            return
        except Exception as e:
            print("Exception:", e)
            import traceback
            traceback.print_exc()
            # continue

    # Write daily CSV data when exiting normally
    for setup_id, data in setup_id_to_data.items():
        print(f"[DEBUG] Normal exit - Setup {setup_id} - daily-data keys: {list(data.get('daily-data', {}).keys())}")
        if data.get("daily-data"):
            path_to_csv_out_dir = config["csv-out"] + str(setup_id) + "/"
            print(f"[DEBUG] Writing daily CSV to: {path_to_csv_out_dir}")
            write_daily_csv(dict(data["daily-data"]), path_to_csv_out_dir)

    print("exiting run_consumer()")
    # debug_file.close()


if __name__ == "__main__":
    run_consumer()
