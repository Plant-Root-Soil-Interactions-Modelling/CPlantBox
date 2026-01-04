import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import math
import visualisation.vtk_plot as vp
import numpy as np

mycp = pb.MycorrhizalPlant()
path = "tomatoparameters/"
name = "TomatoJohanna_WildType"

mycp.readParameters(path + name + ".xml", fromFile = True, verbose = True)

hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.dx = 0.05
hyphae_parameter.distTH = 0.1  # distance for anastomosis 
mycp.setOrganRandomParameter(hyphae_parameter)
# print(hyphae_parameter)

root = mycp.getOrganRandomParameter(pb.root)
for rp in root:
    rp.hyphalEmergenceDensity = 1
    rp.highresolution = 0
    rp.dx = 0.2


mycp.initialize(True)
# print(mycp.toString())
# mycp.writeParameters(name + "_parameters.xml", 'plant', True)

simtime = 50
fps = 2
anim_time = simtime
N = fps * anim_time
dt = simtime / N

filename = "anastomosis_" + str(simtime)


for i in range(0, N):
    print('step',i, '/',N)
    # print(hti)        
    mycp.simulate(dt, False)

ana = pb.SegmentAnalyser(mycp)
ana.addData("AnastomosisPoints", mycp.getAnastomosisPoints(5))
ana.write(filename + ".vtp", True)

# """2D representation of a plant using Matplotlib"""

# import matplotlib.pyplot as plt  # |\label{l2_1g:importStart}|

# import plantbox as pb
# from plantbox.structural.MappedOrganism import MappedPlantPython  # |\label{l2_1g:importEnd}|
# from plantbox.visualisation import figure_style 

# simtime = 14  # [day]  # |\label{l2_1g:defineStart}|
# plant = MappedPlantPython()  # |\label{l2_1g:MappedPlantPython}|
# path = "../../modelparameter/structural/plant/"
# name = "fspm2023"
# plant.readParameters(path + name + ".xml")

# soilSpace = pb.SDF_PlantContainer(500, 500, 500, True)  # to avoid root growing aboveground
# plant.setGeometry(soilSpace)  # creates soil space to stop roots from growing out of the soil

# verbose = False
# plant.initialize(verbose)
# plant.simulate(simtime, verbose)  # |\label{l2_1g:defineEnd}|

# fig, ax = figure_style .subplots11()  # |\label{l2_1g:plotStart}|
# name = ["root", "stem", "leaf", "root tips"]
# color = ["tab:red", "tab:orange", "tab:green", "tab:blue"]

# for idx, ot in enumerate([pb.root, pb.stem, pb.leaf]):
#     pl = plant.getPolylines(ot)  # 3D vectors with coordinates of nodes regrouped per organs  |\label{l2_1g:getPolylines}|
#     label_added = False
#     for node in pl:
#         node = plant.toNumpy(node)  # 'plantbox.Vector3d' to 2D python array  |\label{l2_1g:toNumpy}|
#         ax.plot(
#             node[:, 0],  # x-axis                                         |\label{l2_1g:plot}|
#             node[:, 2],  # z-axis
#             c=color[idx],
#             label=name[idx] if not label_added else None,
#         )
#         label_added = True

# root_tips = plant.get_root_tips()
# ax.scatter(root_tips[:, 0], root_tips[:, 2], c=color[3], label=name[3])

# ax.legend(bbox_to_anchor=(1, 0.5))
# ax.grid(True)
# plt.xlabel("X-axis (cm)")
# plt.ylabel("Depth (cm)")
# ax.relim()
# ax.set_aspect("equal", "box")
# plt.tight_layout()
# plt.savefig("results/example_2_4_2DVisualisation.png")
# plt.show()  # |\label{l2_1g:plotEnd}|