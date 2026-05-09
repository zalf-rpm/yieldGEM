import os
import spotpy
import spotpy_setup_MONICA
import csv
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import colorscali
from collections import defaultdict
#change cur_dir
#start in powershell: 
#> monica-zmq-server -bi -i tcp://*:6666 -bo -o tcp://*:7777
#cur_dir="D:\Projekte_ongoing\Monica\AgMIP_Calibration-master\phase2\step1\MultiExpCalibrator_parallel"
#os.chdir(cur_dir)

font = {'family' : 'calibri',
    'weight' : 'normal',
    'size'   : 18}

def make_lambda(excel):
    return lambda v, p: eval(excel)


cultivar = "bermude"#"apache" #"bermude"

crop_sim_site_MAP = "crop_sim_site_MAP_bermude _original_group1.csv"
observations = "obs_sim_BBCHdoys_bermude_group1.csv"

#cur_dir = os.getcwd() 
os.chdir(os.getcwd())
#read general settings
exp_maps = []
basepath = os.path.dirname(os.path.abspath(__file__))
#cur_dir = os.getcwd() 
with open(os.path.join(basepath, crop_sim_site_MAP)) as exp_mapfile:
    dialect = csv.Sniffer().sniff(exp_mapfile.read(), delimiters=';,\t')
    exp_mapfile.seek(0)
    reader = csv.reader(exp_mapfile, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        exp_map = {}
        exp_map["exp_ID"] = row[0]
        exp_map["sim_file"] = os.path.join(basepath,  row[1])
        exp_map["crop_file"] = os.path.join(basepath,  row[2])
        exp_map["site_file"] = os.path.join(basepath,  row[3])
        exp_map["climate_file"] = os.path.join("/",row[4])
        exp_map["species_file"] = os.path.join(basepath,  row[5])
        exp_map["cultivar_file"] = os.path.join(basepath,  row[6])
        exp_map["where_in_rotation"] = [int(x) for x in row[7].split("-")]
        exp_map["crop_ID"] = row[8]
        exp_map["lat"] = float(row[10])
        exp_map["lon"] =float(row[9])
        exp_map["nfertilizer"] = row[13]
        exp_map["irrigation"] = float(row[14])
        exp_map["sow_Y"] = int(row[15])
        exp_map["sow_M"] = int(row[16])
        exp_map["sow_D"] = int(row[17])
        exp_maps.append(exp_map)

#read observations for which the likelihood of parameter space is calculated
obslist = [] #for envs (outputs)
with open(os.path.join(basepath, observations)) as obsfile:
    dialect = csv.Sniffer().sniff(obsfile.read(), delimiters=';,\t')
    obsfile.seek(0)
    reader = csv.reader(obsfile, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        if row[3].upper() == "Y":
            record = {}
            record["exp_ID"] = row[0]
            record["value"] = float(row[1])
            obslist.append(record) #TODO:Add weight here?

#order obslist by exp_id to avoid mismatch between observation/evaluation lists
def getKey(record):
    return int(record["exp_ID"])
obslist = sorted(obslist, key=getKey)

#read parameters which are to be calibrated
params = []
with open(os.path.join(basepath, "calibratethese_original.csv")) as paramscsv:
    dialect = csv.Sniffer().sniff(paramscsv.read(), delimiters=';,\t')
    paramscsv.seek(0)
    reader = csv.reader(paramscsv, dialect)
    next(reader, None)  # skip the header
    for row in reader:
        p={}
        p["name"] = row[0]
        p["array"] = row[1]
        p["low"] = row[2]
        p["high"] = row[3]
        p["stepsize"] = row[4]
        p["optguess"] = row[5]
        p["minbound"] = row[6]
        p["maxbound"] = row[7]
        if len(row) == 9 and row[8] != "":
            p["derive_function"] = make_lambda(row[8])
        params.append(p)


#Here, MONICA is initialized and a producer is started:
#Arguments are: Parameters, Sites, Observations
#Returns a ready made setup
spot_setup = spotpy_setup_MONICA.spot_setup(params, exp_maps, obslist)
#the same as for example: spot_setup = spot_setup(spotpy.objectivefunctions.rmse)
#Select maximum number of repititions


rep = 10000 #initial number was 10
results = []
#Set up the sampler with the model above
sampler = spotpy.algorithms.sceua(spot_setup, dbname='SCEUA_monica_results', dbformat='csv')
# sampler = sampler = spotpy.algorithms.mcmc(spot_setup, dbname='SCEUA_monica_results', dbformat='csv')

#Run the sampler to produce the paranmeter distribution 
#and identify optimal parameters based on objective function
#ngs = number of complexes
#kstop = max number of evolution loops before convergence
#peps = convergence criterion
#pcento = percent change allowed in kstop loops before convergence 
sampler.sample(rep, ngs=len(params)*2, peps=0.001, pcento=0.001)
# sampler.sample(rep)

#Extract the parameter samples from distribution
results = spotpy.analyser.load_csv_results("SCEUA_monica_results")
print(results)

# Plot how the objective function was minimized during sampling
# fig = plt.figure(1, figsize=(9, 6))
# #plt.plot(results["like1"],  marker='o')
# plt.plot(results["like1"],  'r+')
# plt.show()
# plt.ylabel("RMSE")
# plt.xlabel("Iteration")
# fig.savefig("SCEUA_objectivefunctiontrace_MONICA.png", dpi=150)


# print("sampler_MONICA.py finished")


