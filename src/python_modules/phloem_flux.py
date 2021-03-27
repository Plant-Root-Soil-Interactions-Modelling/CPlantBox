import timeit
import math
import fipy as fp

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import plantbox as pb
from plantbox import XylemFlux
import rsml_reader as rsml  # todo
from Leuning import Leuning 

from fipy import Variable, FaceVariable, CellVariable, Grid1D, ExplicitDiffusionTerm, TransientTerm, DiffusionTerm, Viewer, VTKCellViewer
from fipy.tools import numerix as np
import matplotlib.pyplot as plt
import numpy as np
from fipy.meshes import mesh1D
from fipy.tools.numerix import MA
from Mesh1Dmod import Mesh1Dmod

class PhloemFluxPython(Leuning):
    """  solver for the plant phloem/sucrose flux
         with fipy
         wrapped around the Leuning module which solves the xylem flow and coupled assimilation/transpiration
    """
    def __init__(self, rs):
        """ @param rs is either a pb.MappedPlant or a string containing a rsml filename                     
        """
        super().__init__(rs)
        self.new2oldNodeID = {} #map to get MappedPlant node-ID from grid node-ID
        self.mesh = [] #to store grid
        self.phi = [] # to store solution

    def mp2mesh(self, segs):
        """ Converts a mappedPlant into a network of 1D grids            
            outputs: 
            creates mesh and stores it in self.mesh
            fills self.phi => allows for visualization of the mesh with vtk
        """
        vertices = np.array([]) #containes x,y,z coords of each vertex
        cells = np.array([]) #contains ID of the faces on each cell. here 1D => 1face == 1 vertex
        nodes = self.get_nodes() 
        newNodeID = 0
        for seg in segs:
            seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
            for j, n in enumerate(seg):
                oldNodeID = np.where(np.all(nodes == n, axis=1))[0][0] #find index of n in the node matrix == node ID in the MappedPlant object
                self.new2oldNodeID[newNodeID] =  oldNodeID
                if j > 0: #not first node of organ
                    if cells.size == 0 : #first organ of the plant
                        cells =  np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))) #nodes should be ordered. so seg n goes between nodes j-1 and j
                        
                    else:
                        cells = np.hstack((cells, np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))).reshape(3,1)))
                elif cells.size > 0:
                    if oldNodeID == 0: #if node 0 (root collar): only case where seg linked with to the x node of parent segment 
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[0]] == self.new2oldNodeID[newNodeID])[0][0]
                    else: #seg linked with to the y node of parent segment
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[1]] == self.new2oldNodeID[newNodeID])[0][0]
                    
                    cells[2,seg_parent] = newNodeID #add node to parent seg to link with child branch
                newNodeID += 1
                    
            nodes_x_coord = np.array([xi[0] for xi in seg], dtype=np.float64) 
            nodes_y_coord = np.array([xi[1] for xi in seg], dtype=np.float64) 
            nodes_z_coord = np.array([xi[2] for xi in seg], dtype=np.float64) 
            nodes_coord =  np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored 
            if vertices.size ==0:
                vertices = nodes_coord
            else:
                vertices = np.hstack((vertices, nodes_coord))
        
        faces = np.array((np.arange(0, vertices[0].size),), dtype=np.int64) 
        cells = MA.masked_values(cells, -1)
        self.mesh = Mesh1Dmod(vertexCoords=vertices, faceVertexIDs=faces, cellFaceIDs=cells)
        self.phi = CellVariable(name="solution variable", mesh=self.mesh, value=0.) #if we want to visualize the grid with vtk
        
        