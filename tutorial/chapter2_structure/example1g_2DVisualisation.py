import sys; sys.path.append("../.."); sys.path.append("../../src/") # |\label{l2_1g:importStart}|
import plantbox as pb
import matplotlib.pyplot as plt  
import numpy as np
from structural.MappedOrganism import MappedPlantPython # |\label{l2_1g:importEnd}|

simtime = 14  # [day]  # |\label{l2_1g:defineStart}|
pyPl = MappedPlantPython(pb.MappedPlant()) #for plant objects  |\label{l2_1g:MappedPlantPython}|
pl = pyPl.ms                                            # |\label{l2_1g:pyPlms}|
path = "../../modelparameter/structural/plant/" 
name = "fspm2023" 
pl.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500,500, 500, True) #to avoid root growing aboveground
pl.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil

verbose = False
pl.initialize(verbose)
pl.simulate(simtime,verbose) # |\label{l2_1g:defineEnd}|


fig, ax = plt.subplots() # |\label{l2_1g:plotStart}|
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for idx, ndType in enumerate([pb.root, pb.stem, pb.leaf]):
    organs = pl.getPolylines(ndType) # 3D vectors with coordinates of nodes regrouped per organs  |\label{l2_1g:getPolylines}|
    label_added = False
    for org in organs:
        org = pyPl.toNumpy(org) # 'plantbox.Vector3d' to 2D python array  |\label{l2_1g:toNumpy}|
        ax.plot(org[:,0], # x-axis                                         |\label{l2_1g:plot}|
                org[:,2], # z-axis
                c=color[idx],
                label=name[idx] if not label_added else None
                )
        label_added = True
            

ax.legend()
ax.grid(True)
plt.xlabel("x-axis (cm)")
plt.ylabel("Depth (cm)")
plt.title("2D representation")
plt.show() # |\label{l2_1g:plotEnd}|