""" 2D representation of a plant using Matplotlib"""

import matplotlib.pyplot as plt  # |\label{l2_1g:importStart}|

import plantbox as pb
from plantbox.structural.MappedOrganism import MappedPlantPython  # |\label{l2_1g:importEnd}|

simtime = 14  # [day]  # |\label{l2_1g:defineStart}|
pl = MappedPlantPython()  # |\label{l2_1g:MappedPlantPython}|
path = "../../modelparameter/structural/plant/"
name = "fspm2023"
pl.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500, 500, 500, True)  # to avoid root growing aboveground
pl.setGeometry(soilSpace)  # creates soil space to stop roots from growing out of the soil

verbose = False
pl.initialize(verbose)
pl.simulate(simtime, verbose)  # |\label{l2_1g:defineEnd}|

fig, ax = plt.subplots()  # |\label{l2_1g:plotStart}|
name = ["root", "stem", "leaf", "root tips"]
color = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

for idx, ot in enumerate([pb.root, pb.stem, pb.leaf]):
    organs = pl.getPolylines(ot)  # 3D vectors with coordinates of nodes regrouped per organs  |\label{l2_1g:getPolylines}|
    label_added = False
    for org in organs:
        org = pl.toNumpy(org)  # 'plantbox.Vector3d' to 2D python array  |\label{l2_1g:toNumpy}|
        ax.plot(
            org[:, 0],  # x-axis                                         |\label{l2_1g:plot}|
            org[:, 2],  # z-axis
            c = color[idx],
            label = name[idx] if not label_added else None,
        )
        label_added = True

root_tips = pl.get_root_tips()
ax.scatter(root_tips[:, 0], root_tips[:, 2], c = color[3], label = name[3])

ax.legend(bbox_to_anchor = (1, 0.5))
ax.grid(True)
plt.xlabel("x-axis (cm)")
plt.ylabel("Depth (cm)")
plt.title("2D representation")
ax.relim()
ax.autoscale()
ax.set_aspect("equal", "box")
plt.show()  # |\label{l2_1g:plotEnd}|
