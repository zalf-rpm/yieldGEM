import json
import sys
import monica_io3
import zmq
import csv
import os
from datetime import date,datetime
import collections
import threading
from threading import Thread
from collections import defaultdict
import monica_io3
import numpy as np
import soil_io3
from pyproj import CRS, Transformer
from scipy.interpolate import NearestNDInterpolator
import sqlite3

class monica_adapter(object):



    def __init__(self, exp_maps, obslist):

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


        #for multi-experiment: create a M-2 relationship between exp_IDs and param files
        self.IDs_paramspaths = {}
        for exp_map in exp_maps:
            self.IDs_paramspaths[exp_map["exp_ID"]] = {}
            self.IDs_paramspaths[exp_map["exp_ID"]]["species"] = exp_map["species_file"]
            self.IDs_paramspaths[exp_map["exp_ID"]]["cultivar"] = exp_map["cultivar_file"]

        #observations data structures
        self.observations = [] #for spotpy
        for record in obslist:
            self.observations.append(record["value"])

        self.species_params={} #map to store different species params sets avoiding repetition
        self.cultivar_params={} #map to store different cultivar params sets avoiding repetition

        #create envs
        self.envs = []
        for exp_map in exp_maps:

            config = {
                "path_to_data_dir": "D:/Local/pythonProject/spotpy_phenology_all/data/",
                "monica_path_to_data_dir": "D:/Local/pythonProject/spotpy_phenology_all/data/",
            } 


            wgs84_crs = CRS.from_epsg(4326)
            utm32_crs = CRS.from_epsg(25832)
            latlon_to_utm32n_transformer = Transformer.from_crs(wgs84_crs, utm32_crs, always_xy=True)

            soil_db_con = sqlite3.connect(config["path_to_data_dir"] + "soil/buek200.sqlite")
            path_to_soil_grid = config["path_to_data_dir"] + "soil/buek200_1000_25832_etrs89-utm32n.asc"
            soil_metadata, _ = read_header(path_to_soil_grid)
            soil_grid = np.loadtxt(path_to_soil_grid, dtype=int, skiprows=6)
            soil_interpolate = create_ascii_grid_interpolator(soil_grid, soil_metadata)

            # height data for germany
            path_to_dem_grid = config["path_to_data_dir"] + "dem_1000_25832_etrs89-utm32n.asc" 
            dem_metadata, _ = read_header(path_to_dem_grid)
            dem_grid = np.loadtxt(path_to_dem_grid, dtype=float, skiprows=6)
            dem_interpolate = create_ascii_grid_interpolator(dem_grid, dem_metadata)


            # slope data
            path_to_slope_grid = config["path_to_data_dir"] + "slope_1000_25832_etrs89-utm32n.asc"
            slope_metadata, _ = read_header(path_to_slope_grid)
            slope_grid = np.loadtxt(path_to_slope_grid, dtype=float, skiprows=6)
            slope_interpolate = create_ascii_grid_interpolator(slope_grid, slope_metadata)
     
                
            with open(exp_map["sim_file"]) as simfile:
                sim = json.load(simfile)
                sim["crop.json"] = exp_map["crop_file"]
                sim["site.json"] = exp_map["site_file"]
                #sim["climate.csv"] = exp_map["climate_file"]

            with open(exp_map["site_file"]) as sitefile:
                site = json.load(sitefile)

            with open(exp_map["crop_file"]) as cropfile:
                crop = json.load(cropfile)


            env = monica_io3.create_env_json_from_json_config({
                "crop": crop,
                "site": site,
                "sim": sim,
                "climate": ""
            })
            env["UseAutomaticIrrigation"] = False
            #climate is read by the server
            env["csvViaHeaderOptions"] = sim["climate.csv-options"]
            env["csvViaHeaderOptions"]["start-date"] = sim["climate.csv-options"]["start-date"]
            env["csvViaHeaderOptions"]["end-date"] = sim["climate.csv-options"]["end-date"]
            env["pathToClimateCSV"] = config["monica_path_to_data_dir"] + "climate" + exp_map["climate_file"].replace("txt", "csv")

            print("climate", env["pathToClimateCSV"])
            sr, sh = latlon_to_utm32n_transformer.transform(exp_map["lon"], exp_map["lat"])
            soil_id = int(soil_interpolate(sr, sh))
            soil_profile = soil_io3.soil_parameters(soil_db_con, soil_id)
            if len(soil_profile) == 0:
                print("no soil profile found for soil id: ", soil_id)
                continue
            env["params"]["siteParameters"]["SoilProfileParameters"] = soil_profile

            env["params"]["siteParameters"]["Latitude"] = exp_map["lat"]

            height_nn = float(dem_interpolate(sr, sh))
            env["params"]["siteParameters"]["heightNN"] = height_nn

            slope = float(slope_interpolate(sr, sh))
            env["params"]["siteParameters"]["slope"] = slope / 100.0

            sow_Y = exp_map["sow_Y"]
            sow_M = exp_map["sow_M"]
            sow_D = exp_map["sow_D"]

            sd = datetime(sow_Y, sow_M, sow_D)
            lhd = datetime(sow_Y+1, sow_M, sow_D)
            env["cropRotation"][0]["worksteps"][0]["date"] = "{:02d}-{:02d}-{:02d}".format(sd.year,sd.month, sd.day)
            env["cropRotation"][0]["worksteps"][-1]["latest-date"] = "{:02d}-{:02d}-{:02d}".format(lhd.year,lhd.month, lhd.day)

            # irrigation_amount = exp_map["irrigation"]
            # nfertilizer_amount_list = exp_map["nfertilizer"].strip().split(" ")
            # for i in range(0, len(nfertilizer_amount_list)):
            #      env["cropRotation"][0]["worksteps"][i+1]["amount"][0] = float(nfertilizer_amount_list[i])

            # env["cropRotation"][0]["worksteps"][len(nfertilizer_amount_list)+1]["amount"][0] = irrigation_amount

            position = int(exp_map["where_in_rotation"][0])

            for position in exp_map["where_in_rotation"]:
                for ws in env["cropRotation"][position]["worksteps"]:
                    if ws["type"] == "Seed" or ws["type"] == "Sowing":
                        self.species_params[exp_map["species_file"]] = ws["crop"]["cropParams"]["species"]
                        self.cultivar_params[exp_map["cultivar_file"]] = ws["crop"]["cropParams"]["cultivar"]
                        break

            #monica_io3.add_climate_data_to_env(env, sim) this does not work anymore properly
            
            env["customId"] = exp_map["exp_ID"]
            env["where_in_rotation"] = exp_map["where_in_rotation"]
            self.envs.append(env)

        self.context = zmq.Context()
        self.socket_producer = self.context.socket(zmq.PUSH)
        #self.socket_producer.connect("tcp://cluster2:6666")
        self.socket_producer.connect("tcp://localhost:6666")

    def run(self,args):
        return self._run(*args)

    def _run(self,vector, user_params):

        evallist = []
        self.out = {}

        def seek_set_param(par, p_value, model_params):
           
            p_name = par["name"]
            array = par["array"]
            add_index = False
            if isinstance(model_params[p_name], int) or isinstance(model_params[p_name], float):
                add_index = False
            elif len(model_params[p_name]) > 1 and isinstance(model_params[p_name][1], str):
                add_index = True #the param contains text (e.g., units)
            if array.upper() == "FALSE":
                if add_index:
                    model_params[p_name][0] = p_value
                else:
                    model_params[p_name] = p_value
            else: #param is in an array (possibly nested)
                array = array.split("_") #nested array
                if add_index:
                    array = [0] + array
                if len(array) == 1:
                    model_params[p_name][int(array[0])] = p_value
                elif len(array) == 2:
                    model_params[p_name][int(array[0])][int(array[1])] = p_value
                elif len(array) == 3:
                    model_params[p_name][int(array[0])][int(array[1])][int(array[2])] = p_value
                else:
                    print("param array too nested, contact developers")
            
        #set params according to spotpy sampling. Update all the species/cultivar available
        for i in range(len(user_params)):                        #loop on the user params
            for s in self.species_params:               #loop on the species
                if user_params[i]["name"] in self.species_params[s]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.species_params[s]) if "derive_function" in user_params[i] else vector[i],
                    self.species_params[s])
                else:
                    break                                   #break loop on species if the param is not there
            for cv in self.cultivar_params:                 #loop on the cultivars
                if user_params[i]["name"] in self.cultivar_params[cv]:
                    seek_set_param(user_params[i],
                    user_params[i]["derive_function"](vector, self.cultivar_params[cv]) if "derive_function" in user_params[i] else vector[i],
                    self.cultivar_params[cv])
                else:
                    break
        
        # customize events to get the desired output (DOY BBCH30 and BBCH50, AgMIP-calibration speciic)
        # cv_key = self.cultivar_params.keys()[0]
        # cv_key = list(self.cultivar_params.keys())[0]
        # TSUMS_cv = self.cultivar_params[cv_key]["StageTemperatureSum"][0]

        # Tm1 = TSUMS_cv[1]
        # Tm2 = TSUMS_cv[2]
        # Tm3 = TSUMS_cv[3] + TSUMS_cv[4]

        # dr2stemel = Tm2 * 0.25

        # em2stemel = Tm1 + dr2stemel
 

        # for env in self.envs:

        #     env["events"][0]["while"][2] = Tm1
        #     env["events"][2]["while"][2] = em2stemel
        #     env["events"][4]["while"][2] = Tm1 + Tm2 
        #     env["events"][6]["while"][2] = Tm1 + Tm2 + Tm3 

        #launch parallel thread for the collector
        collector = Thread(target=self.collect_results)
        collector.daemon = True
        collector.start()

        #send jobs to the MONICA server
        for env in self.envs:
            species = self.species_params[self.IDs_paramspaths[env["customId"]]["species"]]
            cultivar = self.cultivar_params[self.IDs_paramspaths[env["customId"]]["cultivar"]]
            for crop_to_cal in env["where_in_rotation"]:
            #if the crop appears more than once in the rotation, the same params will be set
                for ws in env["cropRotation"][int(crop_to_cal)]["worksteps"]:
                    if ws["type"] == "Seed" or ws["type"] == "Sowing":
                        ws["crop"]["cropParams"]["species"] = species
                        ws["crop"]["cropParams"]["cultivar"] = cultivar
                        break

            self.socket_producer.send_json(env)

        #wait until the collector finishes
        collector.join()
        
        #build the evaluation list for spotpy        
        ordered_out = collections.OrderedDict(sorted(self.out.items()))
        for k, v in ordered_out.items():
            for value in v:
                evallist.append(float(value))
           
        return evallist

        
    def collect_results(self):
        socket_collector = self.context.socket(zmq.PULL)
        #socket_collector.connect("tcp://cluster2:7777")
        socket_collector.connect("tcp://localhost:7777")
        received_results = 0
        leave = False
        while not leave:
            try:
                #Start consumer here and save to json output
                rec_msg = socket_collector.recv_json()
            except:
                continue            

            results_rec = []
            for res in rec_msg["data"]:
                
                try:
                    
                    results_rec.append(res["results"][0][0])
                    
                except:
                    print("no results in custom id " + rec_msg["customId"])
            self.out[int(rec_msg["customId"])] = results_rec
            # print (rec_msg["customId"], results_rec)
            received_results += 1
            #print("total received: " + str(received_results))

            if received_results == len(self.envs):
                leave = True
