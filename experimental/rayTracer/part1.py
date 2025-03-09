""" something on PlantVisualiser, cpbvis """
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np

filename = "../../modelparameter/structural/plant/shootOnly_fspm2023.xml"
output = "./results/vis_plant"

time = 28
leaf_res = 30
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename)
vis = pb.PlantVisualiser(plant)


# Initialize
plant.initialize()
vis.SetGeometryResolution(8)
vis.SetLeafResolution(leaf_res)

# Simulate
plant.simulate(time, False)

vis.ResetGeometry()
vis.ComputeGeometryForOrganType(pb.stem, False)
vis.ComputeGeometryForOrganType(pb.leaf, False)

# Write the geometry to file#
data = cpbvis.PolydataFromPlantGeometry(vis)
cpbvis.WritePolydataToFile(data, output )

vp.plot_plant(plant, "subType")

#vp.write_plant("test", plant)  # will write test.vtp (with centerlines), and test_leafs.vtp (as polygones)
print("fin")
