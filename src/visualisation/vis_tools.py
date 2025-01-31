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
  geom = np.reshape(geom, (-1, 3))
  num_points = geom.shape[0]
  geom = vtknp.numpy_to_vtk(geom, deep=True)
  geom.SetName("points")
  points.SetData(geom)
  nodeids = np.array(vis.GetGeometryNodeIds())
  nodeids = vtknp.numpy_to_vtk(nodeids, deep=True)
  nodeids.SetName("nodeids")
  pd.GetPointData().AddArray(nodeids)
  texcoords = np.array(vis.GetGeometryTextureCoordinates())
  texcoords = np.reshape(texcoords, (-1, 2))
  texcoords = vtknp.numpy_to_vtk(texcoords, deep=True)
  texcoords.SetName("tcoords")
  pd.GetPointData().AddArray(texcoords)
  normals = np.array(vis.GetGeometryNormals())
  normals = np.reshape(normals, (-1, 3))
  normals = vtknp.numpy_to_vtk(normals, deep=True)
  normals.SetName("normals")
  pd.GetPointData().AddArray(normals)
  # end for
  cell_data = np.array(vis.GetGeometryIndices())
  cell_data = np.reshape(cell_data, (-1, 3))
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

def WavefrontFromPlantGeometry(vis : pb.PlantVisualiser, plant : pb.Plant, filename: str, resolution = "organ", colour = None, fixedset = None, add_material_text = "") :
  """Create a wavefront object by generating a wavefront file from plant organs"""
  # Create a wavefront file
  wf = open(filename, "w")
  wf.write("# Wavefront object file generated from plant geometry\n")
  wf.write("o plant_geometry\n")
  if colour is not None :
    # cut .obj extension
    mtl = open(filename[:-4] + ".mtl", "w")
    wf.write("mtllib " + filename[:-4] + ".mtl\n")
    for c in colour :
      mtl.write("newmtl type_" + str(c) + add_material_text + "\n")
      mtl.write("Kd " + str(colour[c][0]) + " " + str(colour[c][1]) + " " + str(colour[c][2]) + "\n")
    # end for
    mtl.close()
  # end if
  #vis.SetVerbose(True)
  # Iterate through the organs
  geomset = []
  if resolution == "organ" :
    geomset = [o.getId() for o in plant.getOrgans() if o.organType() != 1][::-1]
  elif resolution == "type" :
    geomset = [4,3,2]
  # end if
  if fixedset is not None :
    geomset = fixedset
  # end if
  vertexcounter = 0 # for running index across sections
  for on in geomset :
    vis.ResetGeometry()
    try :
      if resolution == "organ" :
        vis.ComputeGeometryForOrgan(on)
      elif resolution == "type" :
        vis.ComputeGeometryForOrganType(on)
      #print(vis.SelfCheck())
    except:
      print("Error computing geometry for organ type ", on)
      vis.ResetGeometry()
      continue
    # flush the file to disk
    wf.flush()
    # write the vertices
    vertices = np.reshape(vis.GetGeometry(), (-1, 3))
    # make new group
    wf.write("\tg type_" + str(on) + add_material_text + "\n")
    # make a material entry such that all sections have their own material
    wf.write("\tusemtl type_" + str(on) + add_material_text + "\n")
    # comment metadata
    otype = ""
    if on == pb.stem :
      otype = "stem"
    elif on == pb.leaf :
      otype = "leaf"
    elif on == pb.seed :
      otype = "seed"
    elif on == pb.root :
      otype = "root"
    # end if
    wf.write("\t# Organ type: " + otype + "\n")
    for v in vertices :
      v = v #+ np.random.rand(3)
      wf.write("\t\tv " + str(v[0]) + " " + str(v[1]) + " " + str(v[2]))
      if colour is not None :
        pass
        #wf.write(" " + str(colour[on][0]) + " " + str(colour[on][1]) + " " + str(colour[on][2]))
      wf.write("\n")
    # end for
    #for v in vertices :
    #  wf.write("\t\tvc " + str(colour[on][0]) + " " + str(colour[on][1]) + " " + str(colour[on][2]) + "\n")
    # write the texture coordinates
    texcoords = np.reshape(vis.GetGeometryTextureCoordinates(), (-1, 2))
    for t in texcoords :
      wf.write("\t\tvt " + str(t[0]) + " " + str(t[1]) + "\n")
    # end for
    # write the normals
    normals = np.reshape(vis.GetGeometryNormals(), (-1, 3))
    for n in normals :
      wf.write("\t\tvn " + str(n[0]) + " " + str(n[1]) + " " + str(n[2]) + "\n")
    # end for
    # write the faces#
    indices = np.reshape(vis.GetGeometryIndices(), (-1, 3))
    for i in indices :
      i = i + vertexcounter
      # written as vertex/texcoord/normal
      wf.write("\t\tf " + str(i[0] + 1) + "/" + str(i[0] + 1) + "/" + str(i[0] + 1) + " " + str(i[1] + 1) + "/" + str(i[1] + 1) + "/" + str(i[1] + 1) + " " + str(i[2] + 1) + "/" + str(i[2] + 1) + "/" + str(i[2] + 1) + "\n")
    # end for
    vertexcounter += vertices.shape[0]
  # end for
  wf.close()
# end def WavefrontFromPlantGeometry

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











