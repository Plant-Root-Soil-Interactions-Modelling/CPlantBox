"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import os
import shutil
# save the current working directory
path = os.getcwd()
os.chdir("../../")
# make
os.system("make")
# change back to the working directory
os.chdir(path)

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np
from tqdm import tqdm
import vtk



filename = "../../experimental/photosynthesis/Maiz_PLevels/P3_plant.xml"
output_folder = "./results/"
output = output_folder+"vis_grow"

# remove all files in ./results/ first
import os
filelist = [ f for f in os.listdir(output_folder) if f.endswith(".vtp") and "vis_grow" in f ]
for f in filelist:
  os.remove(os.path.join(output_folder, f))

time = 80
time_resolution = 2
leaf_res = 41
mode = 1 # 0 = time, 1 = end
# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename, verbose = True)
vis = pb.PlantVisualiser(plant)
vis.SetComputeMidlineInLeaf(True)
vis.SetMidlineLeafDiameterFactor(0.5)
vis.SetAddVerticalLeafOffset(True)
vis.SetRestrictBladeWidthAtStem(True)
vis.SetRestrictBladeWidthAtTip(True)
vis.SetVerbose(True)
vis.SetTipWidth(10.0)

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
    vis.ComputeGeometry()
    #vis.ComputeGeometryForOrganType(pb.leaf, False)
    #vis.ComputeGeometryForOrganType(pb.stem, False)
    data = cpbvis.PolydataFromPlantGeometry(vis)
    cpbvis.WritePolydataToFile(data, output + str(i) + ".vtp")
  #endfor
elif mode == 1 :
  plant.simulate(time, False)
  vis.ResetGeometry()
  vis.ComputeGeometry()
  #vis.ComputeGeometryForOrganType(pb.leaf, False)
  #vis.ComputeGeometryForOrganType(pb.stem, False)
  #vis.ComputeGeometry()
  # Write the geometry to file#
  data = cpbvis.PolydataFromPlantGeometry(vis)
  cpbvis.WritePolydataToFile(data, output + ".vtp")

