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

plant.maxLBud = np.array([1.,1.,1.,1.])
plant.maxLBudDormant = np.array([0.15,0.15,0.15,0.15])
plant.budGR = 1.8

sim_time = 1400/20  # days
len_phytomeres = []
numP = []
len_leaves = []

min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])


anim = vp.AnimateRoots(plant)
anim.min = min_
anim.max = max_
anim.res = [1, 1, 1]
anim.plant = True
anim.start()

dt = 0.1
N = int(np.round(sim_time/dt))

for i in range(N):
    plant.simulate(dt)  # Simulate|\label{l13:simulate}|
    # ana = pb.SegmentAnalyser(plant)
    # vp.write_plant(f"results/UQ_1LeafRS{i:04d}", ana)
    anim.update()
    #print(dt * (i+1))
    stem = plant.getOrgans(3)[0]
    epsilonDxPerPhyto = np.array(stem.epsilonDxPerPhyto)
    lln = stem.getLlocalId_linking_nodes()
    print(len(epsilonDxPerPhyto), len(lln))
    lengthP = np.zeros(60)
    __ = np.diff([stem.getLength(lln_i) for lln_i in lln])
    if len(__) > 0:
        __ -= epsilonDxPerPhyto 
    lengthP[:len(__)] = __ /2.7
    lengthkids = np.zeros(60)
    __ = np.array([stem.getChild(np).getLength(False)/5 for np in range(stem.getNumberOfChildren()) if stem.getChild(np).organType() == pb.leaf])
    lengthkids[:len(__)] = __
    len_phytomeres.append(lengthP)
    len_leaves.append(lengthkids)
    print(dt * (i+1), 'sum(lengthP > 0.011)', sum(lengthP > 0.0))
    numP.append(sum(lengthP > 0.0))
    
    
        
len_phytomeres = np.array(len_phytomeres).T
len_leaves = np.array(len_leaves).T


# choose a colormap with enough distinct colors
cmap = plt.get_cmap('tab10')

fig, ax = figure_style.subplots12(1, 1)

# fig supplementary material of Barillotb
for idp, len_phytomere in enumerate(len_phytomeres):
    if idp < 12:
        color = cmap(idp % cmap.N)  # ensures cycling if more than 10
        x = [dt * (i+1) * 20 for i in range(N)]
        ax.plot(x, len_phytomere, color=color, label=f'phytomere {idp}')

for idp, len_leaf in enumerate(len_leaves):
    if idp < 9:
        color = cmap(idp % cmap.N)
        x = [dt * (i+1) * 20 for i in range(N)]
        ax.scatter(x, len_leaf, color=color, label=f'leaf {idp}')

ax.set_xlim([0, 310])
ax.legend()  # optional but useful now
plt.show()


# fig 7 barillot7
fig, ax = figure_style.subplots12(1, 1)
ax.plot([dt * (i+1) * 20 for i in range(N)], numP)
ax.set_xlim([0, 1410])
ax.set_ylim([0, 30])
plt.show()

# ana = pb.SegmentAnalyser(plant)
# ana.write("results/example_plant_segs.vtp")
# # Interactive plot, using vtk
# vp.plot_plant(plant, "organType")  # e.g. organType, subType, age |\label{l13:plot_plant}|
