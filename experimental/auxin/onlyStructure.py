'''
todo: re-compute/check everzthing and write it all down
for structura clibration according to time for the shoot.

shoot created at day 5?
increase a bit the leaf delay?
until bigger leaf opens, the smaller one is hidden inside 
==> no photosynthesis
'''

###################################
# test structure of pea plant
###################################
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
            
import plantbox as pb  # |\label{l13:cplantbox}|
import plantbox.visualisation.vtk_plot as vp  # |\label{l13:vtk_plot}|
from plantbox.visualisation import figure_style
import matplotlib.pyplot as plt  # |\label{l2_2d:importStart}|
import numpy as np

plant = pb.MappedPlant()  # Create a new plant |\label{l13:plant}|

# Open plant and root parameter from a file
path = "./modelparameter/"
name = "UQ_1LeafRS"
plant.readParameters(path + name + ".xml")  # |\label{l13:readparameters}|

plant.initialize()  # Initialize |\label{l13:initialize}|

sim_time = 20  # days
numlats = []
lenmain1 = []
lenmain2 = []
lats = []
latsOT = []


min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])


# anim = vp.AnimateRoots(plant)
# anim.min = min_
# anim.max = max_
# anim.res = [1, 1, 1]
# anim.plant = True
# anim.start()

dt = 0.1
N = int(np.round(sim_time/dt))

for i in range(N):
    plant.simulate(dt)  # Simulate|\label{l13:simulate}|
    ana = pb.SegmentAnalyser(plant)
    vp.write_plant(f"results/UQ_1LeafRS{i:04d}", ana)
    # anim.update()

    numlats.append(plant.getOrgans(3)[0].getNumberOfChildren())
    lenmain1.append(plant.getOrgans(3)[0].getLength(True))
    lenmain2.append(plant.getOrgans(3)[0].getLength(False))
    for kid in range(plant.getOrgans(3)[0].getNumberOfChildren()):
        if kid < len(lats):
            lats[kid].append(plant.getOrgans(3)[0].getChild(kid).getLength(False))
        else:
            lats.append([plant.getOrgans(3)[0].getChild(kid).getLength(False)])
            latsOT.append(plant.getOrgans(3)[0].getChild(kid).organType())
    #print('lats',lats)
# fig, ax = figure_style.subplots12(1, 1)
# ax.plot([i for i in range(sim_time)], numlats)
# plt.show()
# fig, ax = figure_style.subplots12(1, 1)
# ax.plot([i for i in range(sim_time)], lenmain1, label='True')
# ax.plot([i for i in range(sim_time)], lenmain2, label='False')
# plt.show()
# fig, ax = figure_style.subplots12(1, 1)
# for kid in range(plant.getOrgans(3)[0].getNumberOfChildren()):
    # ax.plot([i for i in range(sim_time)], lats[kid], label=str(latsOT[kid]))
# plt.show()

# ana = pb.SegmentAnalyser(plant)
# ana.write("results/example_plant_segs.vtp")
# # Interactive plot, using vtk
# vp.plot_plant(plant, "organType")  # e.g. organType, subType, age |\label{l13:plot_plant}|
