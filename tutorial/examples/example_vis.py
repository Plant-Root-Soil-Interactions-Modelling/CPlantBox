"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np
from tqdm import tqdm
import vtk

filename = "../../experimental/photosynthesis/Maiz_PLevels/P3_plant.xml"

output = "./results/vis_grow"

time = 28
time_resolution = 1
leaf_res = 30
mode = 0 # 0 = time, 1 = end
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename, verbose = True)
vis = pb.PlantVisualiser(plant)
vis.SetComputeMidlineInLeaf(False)

for p in plant.getOrganRandomParameter(pb.seed) :
    p.delayDefinition = 0

# Initialize
plant.initialize()

vis.SetGeometryResolution(8)
vis.SetLeafResolution(leaf_res)

if mode == 0 :
  for i in tqdm(range(time * time_resolution), desc="Sim+Vis", unit="day", unit_scale=True) :
    # Simulate
    plant.simulate(1. / time_resolution, False)
    vis.ResetGeometry()
    vis.ComputeGeometryForOrganType(pb.leaf, False)
    vis.ComputeGeometryForOrganType(pb.stem, False)
    data = cpbvis.PolydataFromPlantGeometry(vis)
    cpbvis.WritePolydataToFile(data, output + str(i) + ".vtp")
  #endfor
elif mode == 1 :
  plant.simulate(time, False)
  vis.ResetGeometry()
  vis.ComputeGeometryForOrganType(pb.leaf, False)
  vis.ComputeGeometryForOrganType(pb.stem, False)
  #vis.ComputeGeometry()
  # Write the geometry to file#
  data = cpbvis.PolydataFromPlantGeometry(vis)
  cpbvis.WritePolydataToFile(data, output + ".vtp")

