"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/"); sys.path.append("./"); sys.path.append("./src/")

import os
import sys
# call cmake --build on ../../
os.system("cmake --build ../../ --config Release -j 4")

import plantbox as pb
import visualisation.vtk_plot as vp
import visualisation.vis_tools as cpbvis

import numpy as np

filename = "../../modelparameter/structural/plant/fspm2023.xml"
output = "./results/vis_plant"

# create a plant
plant = pb.MappedPlant()
plant.readParameters(filename)
vis = pb.PlantVisualiser(plant)


# Initialize
plant.initialize()
vis.SetGeometryResolution(8)
vis.SetLeafResolution(30)
vis.SetComputeMidlineInLeaf(False)

# Simulate
plant.simulate(28, False)

for o in plant.getOrgans() :
  if o.organType() == pb.leaf :
    print("Plant Organ: ", o.getId())

vis.ResetGeometry()
vis.ComputeGeometryForOrgan(2871)

# Write the geometry to file#
data = cpbvis.PolydataFromPlantGeometry(vis)

# view the vtk geometry
import vtk
renderer = vtk.vtkRenderer()
window = vtk.vtkRenderWindow()
window.AddRenderer(renderer)
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)

# Add the data to the renderer
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(data)
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(0.0, 1.0, 1.0)
renderer.AddActor(actor)

# Render
interactor.Initialize()
window.Render()
interactor.Start()


