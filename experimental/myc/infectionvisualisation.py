
import types
import importlib
import os
import sys
SRC_PATH = "../../src/"
sys.path.append("../.."); sys.path.append(SRC_PATH)

# Create a fake plantbox namespace
plantbox = types.SimpleNamespace()

# Automatically import all folders inside src and attach to plantbox
for name in os.listdir(SRC_PATH):
    folder_path = os.path.join(SRC_PATH, name)
    if os.path.isdir(folder_path) and not name.startswith('__'):
        try:
            module = importlib.import_module(name)
            setattr(plantbox, name, module)
            sys.modules[f'plantbox.{name}'] = module
        except ModuleNotFoundError:
            # skip folders that are not importable as modules
            pass
            
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import plantbox.visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()
path = "tomatoparameters/"
name = "TomatoJohanna_WildTypeTwoHyphaeTypes"

# name = "Heliantus_Pagès_2013"

animation = False
local = False
infbox = pb.SDF_PlantBox(4, 4, 4)
# infbox = pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, -10))
for i in range(1,5):
    infbox = pb.SDF_Union(infbox, pb.SDF_RotateTranslate(infbox, 0, 0, pb.Vector3d(0, 0, - i * 10)))

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = mycp.getOrganRandomParameter(pb.hyphae)

# for the first item of the list (0: undefined) 
# we will still have the default value, but the rest is updated
for hp in hyphae_parameter:    
    print(hp.hlt,'a',hp.a,
    'ln',hp.ln,
    'lb',hp.lb,
    'dx',hp.dx,
    'distTH',hp.distTH)

for hp in hyphae_parameter:    
    hp.a = 0.01
    hp.ln = 0.05
    hp.lb = 0.05
    hp.dx = 0.01
    hp.distTH = 0.05   # distance for anastomosis
    mycp.setOrganRandomParameter(hp)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1
    rp.highresolution = 0
    if local:
        rp.f_inf = pb.SoilLookUpSDF(infbox, 0.99, 0.0, 0.1)
    rp.dx = 0.2


mycp.initialize(True)
# print(mycp.toString())
# mycp.writeParameters(name + "_parameters.xml", 'plant', True)

simtime = 10
fps = 24
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "infection_" + str(simtime)
if animation:
    filename = "animation"
if not local:
    filename = filename + "_dispersed"
else:
    filename = filename + "_local"


if animation:
    for i in range(1, N):
        mycp.simulate(dt, True)
        ana = pb.SegmentAnalyser(mycp)
        ana.addData("infection", mycp.getNodeInfections(2))
        ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
        ana.addData("anastomosis", mycp.getAnastomosisPoints(5))
        ana.write("results/" +filename + "{:04d}".format(i) + ".vtp", ["radius", "subType", "creationTime", "organType", "infection", "infectionTime", "anastomosis"])
        print("Frame " + str(i) + " of " + str(N))

else:
    for i in range(1, N+1):
        print('step',i, '/',N)

        
        # print(hti)        
        mycp.simulate(dt, False)

    # vp.plot_plant(mycp, "organType")  
#     print('done')  
ana = pb.SegmentAnalyser(mycp)
ana.addData("infection", mycp.getNodeInfections(2))

vp.plot_plant(ana,"infection")
ana.write(filename + "_hyphalTrees" + ".vtp", ["radius", "subType", "creationTime","organType","hyphalTreeIndex"])
# ana.addData("infection", mycp.getNodeInfections(2))
#     ana.addData("infectionTime", mycp.getNodeInfectionTime(2))
# pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "infection", "infectionTime","organType"])
#     vp.plot_roots(ana, "infection")
#     # vp.plot_roots(ana, "infectionTime")
#     # 
# ana.write(filename + ".vtp", ["radius", "subType", "creationTime","organType"])# "infection", "infectionTime",