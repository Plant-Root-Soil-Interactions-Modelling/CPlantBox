import plantbox as pb
from visualisation.vtk_tools import *

import vtk
import vtk.util.numpy_support as vtknp
import numpy as np

def PolydataFromPlantGeometry(vis : pb.PlantVisualiser ) :
  """Create a vtkPolyData object from a plantbox plant geometry"""
  if not vis.HasGeometry() :
    return None
  # end if
  pd = vtk.vtkPolyData()
  pd.Reset()
  points = vtk.vtkPoints()
  points.Reset()
  geom = np.array(vis.GetGeometry())
  geom = np.reshape(geom, (geom.shape[0]//3, 3))
  num_points = geom.shape[0]
  geom = vtknp.numpy_to_vtk(geom, deep=True)
  geom.SetName("points")
  points.SetData(geom)
  nodeids = np.array(vis.GetGeometryNodeIds())
  print(nodeids.shape)
  nodeids = vtknp.numpy_to_vtk(nodeids, deep=True)
  nodeids.SetName("nodeids")
  pd.GetPointData().AddArray(nodeids)
  texcoords = np.array(vis.GetGeometryTextureCoordinates())
  print(texcoords.shape)
  texcoords = np.reshape(texcoords, (texcoords.shape[0]//2, 2))
  texcoords = vtknp.numpy_to_vtk(texcoords, deep=True)
  texcoords.SetName("tcoords")
  pd.GetPointData().AddArray(texcoords)
  normals = np.array(vis.GetGeometryNormals())
  print(normals.shape)
  normals = np.reshape(normals, (normals.shape[0]//3, 3))
  normals = vtknp.numpy_to_vtk(normals, deep=True)
  normals.SetName("normals")
  pd.GetPointData().AddArray(normals)
  # end for
  cell_data = np.array(vis.GetGeometryIndices())
  cell_data = np.reshape(cell_data, (cell_data.shape[0]//3, 3))
  cells = vtk.vtkCellArray()
  for i in range(cell_data.shape[0]) :
    npcell = cell_data[i,:]
    if any (npcell < 0) or any (npcell >= num_points) :
      print("Error: cell data out of range")
    cells.InsertNextCell(3)
    cells.InsertCellPoint(cell_data[i, 0])
    cells.InsertCellPoint(cell_data[i, 1])
    cells.InsertCellPoint(cell_data[i, 2])
  # end for
  pd.GetPointData().SetTCoords(texcoords)
  pd.GetPointData().SetScalars(nodeids)
  pd.GetPointData().SetNormals(normals)
  pd.SetPoints(points)
  pd.SetPolys(cells)
  # Calculate surface tangents
  tangents = vtk.vtkPolyDataTangents()
  tangents.SetInputData(pd)
  tangents.Update()
  result = tangents.GetOutput()
  return result
# end def PolydataFromPlantGeometry

def WritePolydataToFile(pd : vtk.vtkPolyData, filename : str) :
  """Write a vtkPolyData object to file"""
  writer = vtk.vtkXMLPolyDataWriter()
  writer.SetFileName(filename)
  writer.SetInputData(pd)
  writer.Write()
# end def WritePolydataToFile

def WriteSimulationDataToFile(plant: pb.Plant, filename_prefix : str, times : list = range(30), vis: pb.PlantVisualiser = None, Organs : list = [1,2,3,4]) :
  """This function simulates the plant and writes each time step to a file"""
  # Create a plant visualiser if none is given
  if vis is None :
    vis = pb.PlantVisualiser(plant)
  # end if
  # assume that the plant is at least already initialized and has a parameter file attached
  time = 0
  for t in times :
    time_delta = t - time
    plant.simulate(time_delta)
    time = t
    vis.ResetGeometry()
    for o in Organs :
      vis.ComputeGeometryForOrganType(o)
    # Write the geometry to file
    data = PolydataFromPlantGeometry(vis)
    WritePolydataToFile(data, filename_prefix + "_geometry_" + str(t) + ".vtp")
  # end for
# end def WriteSimulationDataToFile

###########################################################################################
# Functions for visualisation of plant geometry in unreal engine
###########################################################################################

# note that this section requires the Synavis plugin for unreal engine and Synavis to be built
# and installed on the system

def CheckForSynavis(path : str = None) :
  """Check if the Synavis plugin is available"""
  # check if synavis is in the environment variables
  import os
  import sys
  if "SYNAVIS_PATH" in os.environ :
    sys.path.append(os.environ["SYNAVIS_PATH"])
  # end if
  if path is not None :
    sys.path.append(path)
  # end if
  try :
    import synavis
    return True
  except :
    return False
  # end try
# end def CheckForSynavis











