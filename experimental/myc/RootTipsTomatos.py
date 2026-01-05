
import matplotlib.pyplot as plt 
import sys; sys.path.append("../.."); sys.path.append("../../src/")
import plantbox as pb
from plantbox.structural.MappedOrganism import MappedPlantPython  
from plantbox.visualisation import figure_style 

simtime = 14  # [day]  
plant = MappedPlantPython()  
path = "tomatoparameters/"
name = "TomatoJohanna_WildType"
plant.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500, 500, 500, True) 
plant.setGeometry(soilSpace) 

verbose = False
plant.initialize(verbose)
plant.simulate(simtime, verbose) 

fig, ax = figure_style .subplots11()  
name = ["root", "stem", "leaf", "root tips"]
color = ["tab:red", "tab:orange", "tab:green", "tab:blue"]

for idx, ot in enumerate([pb.root, pb.stem, pb.leaf]):
    pl = plant.getPolylines(ot) 
    label_added = False
    for node in pl:
        node = plant.toNumpy(node)  
        ax.plot(
            node[:, 0],  # x-axis                                        
            node[:, 2],  # z-axis
            c=color[idx],
            label=name[idx] if not label_added else None,
        )
        label_added = True

root_tips = plant.get_root_tips()
ax.scatter(root_tips[:, 0], root_tips[:, 2], c=color[3], label=name[3])

ax.legend(bbox_to_anchor=(1, 0.5))
ax.grid(True)
plt.xlabel("X-axis (cm)")
plt.ylabel("Depth (cm)")
ax.relim()
ax.set_aspect("equal", "box")
plt.tight_layout()
plt.savefig("results/tomatoWildType_2DVisualisation.png")
plt.show()
