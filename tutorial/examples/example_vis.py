"""root system length over time"""
#sys.path.append("../.."); sys.path.append("../../src/"); 
import sys; sys.path.append("./"); sys.path.append("./src/")

import os
import sys
# call cmake --build on ../../
os.system("cmake --build . --config Release -j 4")

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
plant.setSeed(2)
plant.initialize(False, True)
vis.SetGeometryResolution(8)
vis.SetLeafResolution(30)
vis.SetComputeMidlineInLeaf(False)
vis.SetVerbose(True)
vis.SetLeafMinimumWidth(0.001)
vis.SetRightPenalty(0.5)
vis.SetShapeFunction(lambda t : 2*((1 - t**2)**0.5))
vis.SetLeafWidthScaleFactor(1.0)

# Simulate
plant.simulate(50, False)

leaf_id = -1
leaf_lenght = -1
leaf_ptr = None

for o in plant.getOrgans() :
  if o.getLength(False) > leaf_lenght and o.organType() == pb.leaf :
    leaf_id = o.getId()
    leaf_lenght = o.getLength(True)
    leaf_ptr = o


vis.ResetGeometry()
vis.ComputeGeometryForOrgan(leaf_id)

uv = vis.GetGeometryTextureCoordinates()
points = vis.GetGeometry()
points = np.array(points)
points = points.reshape(-1, 3)
print("Shape of points:")
print(points.shape)

# Write the geometry to file#
data = cpbvis.PolydataFromPlantGeometry(vis)
cpbvis.WritePolydataToFile(data, "/mnt/f/Work/CPlantBox/dumux/CPlantBox/tutorial/examples/leaf.vtk")

# view the vtk geometry
import vtk
renderer = vtk.vtkRenderer()
window = vtk.vtkRenderWindow()
window.AddRenderer(renderer)
interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(window)
# make window large
window.SetSize(800, 800)

nodes = np.array([[v.x,v.y,v.z] for v in leaf_ptr.getNodes()])
spline = pb.CatmullRomSplineManager()
spline.setY(leaf_ptr.getNodes())

# function to add a sphere at a position
def add_sphere(renderer, position, color, radius=0.1):
  position = np.array([position.x, position.y, position.z]) if isinstance(position, pb.Vector3d) else position
  sphere = vtk.vtkSphereSource()
  sphere.SetCenter(position)
  sphere.SetRadius(radius)
  sphere.SetPhiResolution(100)
  sphere.SetThetaResolution(100)
  mapper = vtk.vtkPolyDataMapper()
  mapper.SetInputConnection(sphere.GetOutputPort())
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(color)
  renderer.AddActor(actor)

# make a sphere for each node
#for n in nodes:
#   add_sphere(renderer, n, [1.0, 0.0, 0.0], 0.2)

for t in np.linspace(0, 1, 100):
  n = spline(t)
  add_sphere(renderer, n, [0.0, 1.0, 0.0], 0.1)

add_sphere(renderer, spline.help_lower(), [0.5, 0.5, 0.0], 0.4)
add_sphere(renderer, spline.help_upper(), [0.0, 0.5, 0.5], 0.1)
num = spline.size()
add_sphere(renderer, spline.getControlPoint(num-1), [1.0, 0.0, 0.0], 0.2)
add_sphere(renderer, spline.getControlPoint(num-2), [1.0, 0.0, 0.0], 0.2)
add_sphere(renderer, spline.getControlPoint(num-3), [1.0, 0.0, 0.0], 0.2)

print("T Values of last spline: ", spline.getTValues(spline.splineSize() - 1))

Ts = np.array(spline.getT())
print("Ts: ", Ts)

# Add the data to the renderer
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(data)
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(0.0, 1.0, 1.0)
# make actor emissive
actor.GetProperty().SetAmbient(1)
renderer.AddActor(actor)
# show the wireframe
actor.GetProperty().SetRepresentationToWireframe()
actor.GetProperty().SetLineWidth(2)

# include uniform lighting
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(0, 0, 1)
renderer.AddLight(light)
# include uniform lighting
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(1, 0, 0)
renderer.AddLight(light)
# include uniform lighting
light = vtk.vtkLight()
light.SetFocalPoint(0, 0, 0)
light.SetPosition(0, 1, 0)
renderer.AddLight(light)


# Render
interactor.Initialize()
window.Render()
interactor.Start()


