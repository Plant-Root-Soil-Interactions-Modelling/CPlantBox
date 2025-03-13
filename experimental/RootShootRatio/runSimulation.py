import sys;

CPBdir = "../.."
sys.path.append(CPBdir + "/src");
sys.path.append(CPBdir);
sys.path.append("../../..");
sys.path.append("..");
sys.path.append(CPBdir + "/src/python_modules");
sys.path.append("../build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../build-cmake/cpp/python_binding/")
sys.path.append("../modules/")  # python wrappers

import importlib
import pandas as pd

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis
import os
import numpy as np

def write_file_array(name, data, directoryN):
    name2 =  './results/'+directoryN+ name+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')
        
def write_file_float(name, data, directoryN):
    name2 =  './results/'+directoryN+  name+ '.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')
        
directoryN ="empMaiz/"

main_dir=os.environ['PWD']#dir of the file
results_dir = main_dir +"/results/"+directoryN
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
    import shutil
    shutil.rmtree(results_dir)
    os.makedirs(results_dir)
simT = 0          
simMax = 200  # simStartSim+ spellDuration
depth = 60
dt = 1# / 24  # 10min
verbose = True

N = int(simMax/dt)

# plant system
pl = pb.MappedPlant(seednum=2)  # pb.MappedRootSystem() #pb.MappedPlant()
# pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = CPBdir + "/modelparameter/structural/plant/"
name ="P3" # "Triticum_aestivum_test_2021_1cm"# "Triticum_aestivum_test_2021_1cm"

pl.readParameters( name + ".xml")#path +

# raise Exception
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth)

pl.setGeometry(sdf)  # creates soil space to stop roots from growing out of the soil

pl.initialize(verbose=True)  # , stochastic = False)

    



for i in range(N):
    print('loop',i+1,'/',N)
    pl.simulate(dt, False)
    simT += dt
    

    orgs_all = pl.getOrgans(-1, True)
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all])
    ot_orgs = np.array([org.organType() for org in orgs_all])
    st_orgs = np.array([org.getParameter("subType") for org in orgs_all])
    write_file_array("volOrg",  volOrg, directoryN) #with epsilon
    write_file_array("st_orgs", st_orgs, directoryN)
    write_file_array("ot_orgs", ot_orgs, directoryN)
    write_file_float("time", simT, directoryN)


print("fin")