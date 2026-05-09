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

from datetime import datetime, timedelta
import json
import monica_io3
import numpy as np
import soil_io3
from pyproj import CRS, Transformer
from scipy.interpolate import NearestNDInterpolator
import sqlite3
import sys
import zmq

def create_ascii_grid_interpolator(grid, meta_data, ignore_nodata=True):
    "read an ascii grid into a map, without the no-data values"
    "grid - 2D array of values"

    rows, cols = grid.shape

    cellsize = int(meta_data["cellsize"])
    xll = int(meta_data["xllcorner"])
    yll = int(meta_data["yllcorner"])
    nodata_value = meta_data["nodata_value"]

    xll_center = xll + cellsize // 2
    yll_center = yll + cellsize // 2
    yul_center = yll_center + (rows - 1)*cellsize

    points = []
    values = []

    for row in range(rows):
        for col in range(cols):
            value = grid[row, col]
            if ignore_nodata and value == nodata_value:
                continue
            r = xll_center + col * cellsize
            h = yul_center - row * cellsize
            points.append([r, h])
            values.append(value)

    return NearestNDInterpolator(np.array(points), np.array(values))


def read_header(path_to_ascii_grid_file):
    "read metadata from esri ascii grid file"
    metadata = {}
    header_str = ""
    with open(path_to_ascii_grid_file) as _:
        for i in range(0, 6):
            line = _.readline()
            header_str += line
            sline = [x for x in line.split() if len(x) > 0]
            if len(sline) > 1:
                metadata[sline[0].strip().lower()] = float(sline[1].strip())
    return metadata, header_str


def run_producer(server = {"server": None, "port": None}):

    context = zmq.Context()
    socket = context.socket(zmq.PUSH) # pylint: disable=no-member

    config = {
        "port": server["port"] if server["port"] else "6666",
        "server": server["server"] if server["server"] else "localhost",
        "out": "/out",
        "path_to_data_dir": "D:/Local/pythonProject/YieldGEM_all/data/",
        "monica_path_to_data_dir": "D:/Local/pythonProject/YieldGEM_all/data/",
        # "path_to_data_dir": "data/",
        # "monica_path_to_data_dir": "/home/berg/Desktop/ehsan/monica_course/data/",
        "writenv": False,
    }
    # read commandline args only if script is invoked directly from commandline
    if len(sys.argv) > 1 and __name__ == "__main__":
        for arg in sys.argv[1:]:
            k, v = arg.split("=", maxsplit=1)
            if k in config:
                config[k] = v.lower() == "true" if v.lower() in ["true", "false"] else v 
    print("config:", config)
    
    socket.connect("tcp://" + config["server"] + ":" + config["port"])

    wgs84_crs = CRS.from_epsg(4326)
    utm32_crs = CRS.from_epsg(25832)
    latlon_to_utm32n_transformer = Transformer.from_crs(wgs84_crs, utm32_crs, always_xy=True)

    soil_db_con = sqlite3.connect(config["path_to_data_dir"] + "soil/buek200.sqlite")
    path_to_soil_grid = config["path_to_data_dir"] + "soil/buek200_1000_25832_etrs89-utm32n.asc"
    soil_metadata, _ = read_header(path_to_soil_grid)
    soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
    soil_interpolate = create_ascii_grid_interpolator(soil_grid, soil_metadata)
    print("read: ", path_to_soil_grid)

    # height data for germany
    path_to_dem_grid = config["path_to_data_dir"] + "dem_1000_25832_etrs89-utm32n.asc" 
    dem_metadata, _ = read_header(path_to_dem_grid)
    dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=6)
    dem_interpolate = create_ascii_grid_interpolator(dem_grid, dem_metadata)
    print("read: ", path_to_dem_grid)

    # slope data
    path_to_slope_grid = config["path_to_data_dir"] + "slope_1000_25832_etrs89-utm32n.asc"
    slope_metadata, _ = read_header(path_to_slope_grid)
    slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
    slope_interpolate = create_ascii_grid_interpolator(slope_grid, slope_metadata)
    print("read: ", path_to_slope_grid)




    with open("site.json") as _:
        site_json = json.load(_)

    # with open("crop.json") as _:
    #     crop_json = json.load(_)

    soil_data = []
    with open("Soil_coordinates_1km.csv") as _:
        for line in _.readlines()[1:]:
            es = line.split(",")
            soil_data.append({
                "id": int(es[0]),
                "filename": es[1],
                "lat": float(es[3]), 
                "lon": float(es[2]),
                "field_id": float(es[4]),
                "sowing_day": float(es[5]),
                "nfertilizer": (es[6]),
                "irrigation": float(es[7]),
                "sow_Y": int(es[8]),
                "sow_M": int(es[9]),
                "sow_D": int(es[10]),
                "crop_id": es[12].strip()
            })

    # sowing_doy = [244]
    nfertilizer_amounts = [100]
    irrigation_amounts = [100]

    # env = monica_io3.create_env_json_from_json_config({
    #     "crop": crop_json,
    #     "site": site_json,
    #     "sim": sim_json,
    #     "climate": "" #climate_csv
    # })

    env_count = 0
    for data in soil_data:
        crop_path = data["crop_id"]
        # crop_path = "data/crops/crop.json"
        with open(crop_path) as _:
            crop_json = json.load(_)

        sim_path = crop_path.replace("crops", "sims")
        # sim_path = "sim.json"
        with open(sim_path) as _:
            sim_json = json.load(_)

        env = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""  # climate_csv
        })
        env["UseAutomaticIrrigation"] = False
        env["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
        env["pathToClimateCSV"] = config["monica_path_to_data_dir"] + "climate1/" + data["filename"].replace("txt", "csv")

        sr, sh = latlon_to_utm32n_transformer.transform(data["lon"], data["lat"])

        soil_id = int(soil_interpolate(sr, sh))

        soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
        if len(soil_profile) == 0:
            print("no soil profile found for soil id: ", soil_id)

            continue
        env["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

        env["params"]["siteParameters"]["Latitude"] = data["lat"]

        height_nn = float(dem_interpolate(sr, sh))
        env["params"]["siteParameters"]["heightNN"] = height_nn

        slope = float(slope_interpolate(sr, sh))
        env["params"]["siteParameters"]["slope"] = slope / 100.0

        sow_Y = data["sow_Y"]
        sow_M = data["sow_M"]
        sow_D = data["sow_D"]
  
        sd = datetime(sow_Y, sow_M, sow_D)
        lhd = datetime(sow_Y+1, sow_M, sow_D)
        
        env["cropRotation"][0]["worksteps"][0]["date"] = "{:02d}-{:02d}-{:02d}".format(sd.year,sd.month, sd.day)
        env["cropRotation"][0]["worksteps"][-1]["latest-date"] = "{:02d}-{:02d}-{:02d}".format(lhd.year,lhd.month, lhd.day)

        


        env["customId"] = {
            "env-id": env_count,
            "id": data["id"],
            # "sowing_doy": sdoy,
        }

        socket.send_json(env)

        env_count += 1
        print("sent env:", env_count)

    print("done")

if __name__ == "__main__":
    run_producer()