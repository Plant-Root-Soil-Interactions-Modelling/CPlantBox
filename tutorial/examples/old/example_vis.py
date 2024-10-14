"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np

filename = "../../modelparameter/structural/plant/fspm2023.xml"
output = "./results/vis_plant"

time = 28
leaf_res = 30
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename)
# vis = pb.PlantVisualiser(plant)

# for p in plant.getOrganRandomParameter(pb.stem):
#   p.r = 0.758517633
#   p.lb = 4
#   # p.delayLat = 4
#   p.lmax = (time - 7) * p.r
#   p.dx = 0.1
# for p in plant.getOrganRandomParameter(pb.leaf):
#   p.la, p.lmax = 38.41053981, 38.41053981
#   # p.theta = 0.2 # 0.2
#   p.theta = 0.01
#   p.tropismT = 1  # 6: Anti-gravitropism to gravitropism
#   p.areaMax = 54.45388021  # cm2, area reached when length = lmax
#   phi = np.array([-90, -80, -45, 0., 45, 90]) / 180. * np.pi
#   l = np.array([38.41053981, 1 , 1, 0.3, 1, 38.41053981])  # distance from leaf center
#   p.leafGeometryPhi = phi
#   p.leafGeometryX = l
#   # p.tropismN = 5
#   p.tropismS = 0.08
#   p.tropismAge = 5  # < age at which tropism switch occures, only used if p.tropismT = 6
#   p.createLeafRadialGeometry(phi, l, leaf_res)

# Initialize
plant.initialize()
# vis.SetGeometryResolution(8)
# vis.SetLeafResolution(leaf_res)

# Simulate
plant.simulate(time, True)

# vis.ResetGeometry()
# vis.ComputeGeometryForOrganType(pb.stem, False)
# vis.ComputeGeometryForOrganType(pb.leaf, False)

# Write the geometry to file#
# data = cpbvis.PolydataFromPlantGeometry(vis)
# cpbvis.WritePolydataToFile(data, output + ".vtp")

# vp.plot_plant(plant, "subType")

vp.write_plant("test", plant)  # will write test.vtp (with centerlines), and test_leafs.vtp (as polygones)
print("fin")
