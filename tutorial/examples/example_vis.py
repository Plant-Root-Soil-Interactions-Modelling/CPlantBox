"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np
from tqdm import tqdm

filename = "../../experimental/photosynthesis/Maiz_PLevels/P3_plant.xml"

output = "./results/vis_grow"

time = 28
time_resolution = 24
leaf_res = 30
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename)
vis = pb.PlantVisualiser(plant)
vis.SetComputeMidlineInLeaf(True)

# Initialize
plant.initialize()

vis.SetGeometryResolution(8)
vis.SetLeafResolution(leaf_res)
for i in tqdm(range(time * time_resolution), desc="Sim+Vis", unit="day", unit_scale=True) :
  # Simulate
  plant.simulate(1. / time_resolution, False)
  vis.ResetGeometry()
  vis.ComputeGeometryForOrganType(pb.leaf, False)
  vis.ComputeGeometryForOrganType(pb.stem, False)
  # Write the geometry to file#
  data = cpbvis.PolydataFromPlantGeometry(vis)
  cpbvis.WritePolydataToFile(data, output + str(i) + ".vtp")

