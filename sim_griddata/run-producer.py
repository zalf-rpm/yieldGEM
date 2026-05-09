#!/usr/bin/python
# -*- coding: UTF-8


from datetime import datetime, timedelta
from pathlib import Path
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

def set_initial_n(soil_profile, nh4, no3):
    """设置土壤层的初始NH4和NO3浓度"""
    for layer in soil_profile:
        if "SoilAmmonium" not in layer:
            layer["SoilAmmonium"] = [float(nh4), "kg NH4-N m-3"]
        if "SoilNitrate" not in layer:
            layer["SoilNitrate"] = [float(no3), "kg NO3-N m-3"]
    return soil_profile

def set_initial_water(soil_profile, fc_percent):
    """
    设置土壤初始含水量（相对于田间持水量的百分比）

    Args:
        soil_profile: 土壤层列表
        fc_percent: 0-100 的百分比值
                   - 30-50: 干旱条件
                   - 60-80: 正常田间含水量（推荐）
                   - 90-100: 饱和或近饱和（积水风险）

    Monica 计算公式: vs_SoilMoisture_m3 = FieldCapacity * fc_percent / 100.0
    """
    for layer in soil_profile:
        if "SoilMoisturePercentFC" not in layer:
            layer["SoilMoisturePercentFC"] = [float(fc_percent), "%"]
    return soil_profile

def run_producer(server = {"server": None, "port": None}):

    context = zmq.Context()
    socket = context.socket(zmq.PUSH) # pylint: disable=no-member

    script_dir = Path(__file__).parent
    data_dir = script_dir / "data"

    import os
    os.chdir(script_dir)

    config = {
        "port": server["port"] if server["port"] else "6666",
        "server": server["server"] if server["server"] else "localhost",
        "out": "/out",
        "path_to_data_dir": str(data_dir) + "/",
        "monica_path_to_data_dir": str(data_dir) + "/",

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



    soil_data = []
    setup_file = script_dir / "Set_up.csv"
    with open(setup_file) as _:
        for line in _.readlines()[1:]:
            es = line.split(",")
            soil_data.append({
                "id": int(es[0]),
                "filename": es[1],
                "lat": float(es[3]), 
                "lon": float(es[2]),
                "field_id": int(es[4]),
                "sowing_day": float(es[5]),
                "nfertilizer": (es[6]),
                "irrigation": float(es[7]),
                "sow_Y": int(es[8]),
                "sow_M": int(es[9]),
                "sow_D": int(es[10]),
                "crop_id": es[12].strip()
            })


    env_count = 0
    for data in soil_data:
        crop_path = Path(script_dir) / data["crop_id"]
        with open(crop_path) as _:
            crop_json = json.load(_)

        sim_path = Path(str(crop_path).replace("crops", "sims"))
        with open(sim_path) as _:
            sim_json = json.load(_)

        site_path = Path(str(crop_path).replace("crops", "sites"))
        with open(site_path) as _:
            site_json = json.load(_)

        sim_json["UseAutomaticIrrigation"] = False
        sim_json["NitrogenResponseOn"] = True
        sim_json["WaterDeficitResponseOn"] = True
        sim_json["EmergenceMoistureControlOn"] = False
        sim_json["EmergenceFloodingControlOn"] = False



        env = monica_io3.create_env_json_from_json_config({
            "crop": crop_json,
            "site": site_json,
            "sim": sim_json,
            "climate": ""  # climate_csv
        })
        
        env["csvViaHeaderOptions"] = sim_json["climate.csv-options"]
        env["pathToClimateCSV"] = config["monica_path_to_data_dir"] + "climate/" + data["filename"].replace("txt", "csv")
    
        sr, sh = latlon_to_utm32n_transformer.transform(data["lon"], data["lat"])


        # 设置初始 NH4 和 NO3（单位：kg NH4-N/NO3-N m-3）
        # 典型范围：NH4=1-3，NO3=5-15（根据土壤肥力调整）
        set_initial_n(env["params"]["siteParameters"]["SoilProfileParameters"], nh4=0.001, no3=0.001)

        # 设置初始土壤含水量（相对于田间持水量的百分比）
        # 60-80%: 适合萌发的最佳条件
        # set_initial_water(env["params"]["siteParameters"]["SoilProfileParameters"], fc_percent=80)


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