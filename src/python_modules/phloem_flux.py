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
from scipy.linalg import norm

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
        self.old2newNodeID = {} #map to get grid node-ID from MappedPlant node-ID 
        self.mesh = [] #to store grid
        self.phi = [] # to store solution
        self.faces = []
        self.osmoCoeff = 1 #so for correspond to R*T * mean_phi
        self.D = 1
        self.exteriorCoeff =1.
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.k1 = 0.0001
        self.k2 = 0.01
        self.Gr = 0.05 #C-limited growth during time step
        self.GrEff = 0.5 #Growth efficiency ( = growth per unit of phi used)
        self.k_growth = 0.001
        self.phiThrGrowth = 0.001
        self.lengthThr = 20 #max length
        self.rhoC = 1 # C content per unit of volume
        self.volumes = 0. #volume of each cell/segment
        self.RmMax = 0.

    def mp2mesh(self, segs): #create an old2newNodeID map?
        """ Converts a mappedPlant into a network of 1D grids            
            outputs: 
            creates mesh and stores it in self.mesh
            fills self.phi => allows for visualization of the mesh with vtk
        """
        vertices = np.array([]) #containes x,y,z coords of each vertex
        cells = np.array([]) #contains ID of the faces on each cell. here 1D => 1face == 1 vertex
        nodes = self.get_nodes() 
        newNodeID = 0
        facesNorm = {} #map face to normal vector
        organTypes = self.get_organ_types()
        nOld = nodes[0]
        radii = np.array([])
        length = np.array([])
        for seg in segs:
            seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
            for j, n in enumerate(seg):
                oldNodeID = np.where(np.all(nodes == n, axis=1))[0][0] #find index of n in the node matrix == node ID in the MappedPlant object
                ot = 2
                if (oldNodeID > 0):
                    ot = int(organTypes[oldNodeID - 1]) 
                self.new2oldNodeID[newNodeID] =  oldNodeID
                temp = self.old2newNodeID.get(oldNodeID,np.array([], dtype=np.int64))
                self.old2newNodeID[oldNodeID] = np.append(temp, np.array([newNodeID], dtype=np.int64))
                normal = (nOld - n)/norm(nOld - n)
                #print('normal ', nOld, n, norm(nOld - n), normal)
                if j > 0: #not first node of organ
                    radii = np.concatenate((radii, [self.rs.radii[oldNodeID - 1]]))
                    length = np.concatenate( (length, [norm(nOld - n)]))
                    if cells.size == 0 : #first organ of the plant
                        cells =  np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))) #nodes should be ordered. so seg n goes between nodes j-1 and j
                        facesNorm[newNodeID] = normal
                        facesNorm[newNodeID - 1] = -normal #replace value for all node except for the one at the tip of organs
                        
                        #print(newNodeID - 1, -normal,newNodeID, normal)
                    elif(ot <3):
                        cells = np.hstack((cells, np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))).reshape(3,1)))
                        facesNorm[newNodeID - 1] = -normal #replace value for all node except for the one at the tip of organs
                        facesNorm[newNodeID] = normal
                        #print(newNodeID - 1, -normal,newNodeID, normal)
                    else: #stem of leaf
                        cells = np.hstack((cells, np.array((np.array([newNodeID]),np.array([newNodeID -1]),np.array([-1]))).reshape(3,1)))
                        facesNorm[newNodeID] = normal
                        #print(newNodeID , normal)
                elif cells.size > 0:
                    if oldNodeID == 0: #if node 0 (root collar): only case where seg linked with to the x node of parent segment 
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[0]] == self.new2oldNodeID[newNodeID])[0][0]
                        facesNorm[newNodeID] = facesNorm[0]
                        #print(newNodeID , facesNorm[0])
                    else: #seg linked with to the y node of parent segment
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[1]] == self.new2oldNodeID[newNodeID])[0][0]
                        if(ot >=3):
                            facesNorm[newNodeID] = facesNorm[cells[0,seg_parent]]
                            #print(newNodeID , facesNorm[newNodeID])
                    
                    cells[2,seg_parent] = newNodeID #add node to parent seg to link with child branch
                #print(self.facesNorm)
                newNodeID += 1
                nOld = n    
            nodes_x_coord = np.array([xi[0] for xi in seg], dtype=np.float64) 
            nodes_y_coord = np.array([xi[1] for xi in seg], dtype=np.float64) 
            nodes_z_coord = np.array([xi[2] for xi in seg], dtype=np.float64) 
            nodes_coord =  np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored 
            if vertices.size ==0:
                vertices = nodes_coord
            else:
                vertices = np.hstack((vertices, nodes_coord))
        self.faces = np.array((np.arange(0, vertices[0].size),), dtype=np.int64) 
        cells = MA.masked_values(cells, -1)
        self.mesh = Mesh1Dmod(facesNorm=facesNorm,radii = radii,length = length, vertexCoords=vertices, faceVertexIDs=self.faces, cellFaceIDs=cells)
        self.volumes = np.array([(np.pi * self.mesh.radii[xi]**2) *self.mesh.length[xi] for xi in range(self.mesh.numberOfCells) ])
        self.D = FaceVariable(name = 'diffusion coefficient', mesh = self.mesh, value = self.osmoCoeff)
        self.Source = CellVariable(name = 'An', mesh = self.mesh, value = 0.)
        self.Rm = CellVariable(name = 'Rm', mesh = self.mesh, value = 0.)
        self.GrSink = CellVariable(name = 'GrSink', mesh = self.mesh, value = 0.)
        self.GrC = CellVariable(name = 'GrC', mesh = self.mesh, value = 0.)
        self.GrEff= CellVariable(name = 'GrEff', mesh = self.mesh, value = 0.5)
        self.exteriorCoeff = FaceVariable(self.mesh,value=0.,rank=1)
        self.phi = CellVariable(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        self.setD() #no diffusion to/from tips 
        #self._setFreeOutflow() #=>  free flow BCs for root-tips
        self.setAnSource() #to reset each time we have new An value
        self.setRmSink()
        self.setGrSink()
        
    def _setFreeOutflow(self): #or set dirichlet of 0 on outer face? or fixed gradient?
        #self.exteriorCoeff.setValue(0.) #reset for security 
        tiproots, tipstem, tipleaf = self.get_organ_nodes_tips()
        tiproots_newID = np.array([self.old2newNodeID[xi] for xi in tiproots]).flatten()
        tiproots=  [np.any([tiproots_newID == yi]) for yi in self.mesh.faceVertexIDs[0]]
        print(len(tiproots), np.where(tiproots))
        self.phi.constrain(0., where = tiproots)

    def setAnSource(self): #sourceterm in leaf segment
        self.Source.setValue(0.) #= CellVariable(name = 'An', mesh = self.mesh, value = 0.) #reset for security
        # find seg index with leaf, should be ordered like the An vector
        nodeYLeaf = self.get_segments_index(4) + 1 #ynode_index of leaves
        nodeYLeaf_newID = np.array([self.old2newNodeID[xi] for xi in nodeYLeaf]).flatten() #new index
        for i, an in enumerate(self.An):
                self.Source.constrain(an, where= self.mesh.cellFaceIDs[0] == nodeYLeaf_newID[i]) #0 node = y node for stem and leaves

    def setD(self):
        self.D.setValue(self.osmoCoeff) 
        tiproots, tipstem, tipleaf = self.get_organ_nodes_tips()
        #tips = np.concatenate((tipstem, tipleaf)).flatten()
        tips = np.concatenate((self.get_organ_nodes_tips())).flatten()
        tips_newID = np.array([self.old2newNodeID[xi] for xi in tips]).flatten()
        tips=  np.array(([np.any([tips_newID == yi]) for yi in self.mesh.faceVertexIDs[0]])).flatten()
        self.D.constrain(0., where = tips) 
      
    def setRmSink(self): # Daudet 2002: Rm = (k1 + k2 * C)Sr. set k1 = 0
        self.RmMax = self.rhoC * self.volumes * (self.phi * self.k1 + self.k2 )#  * (self.phi - self.phiThrGrowth) error when use a k2
        self.CSat = self.RmMax < self.phi*0.9 #need to add 0.9 factor or phi goes < 0
        self.Rm = (self.RmMax * (self.RmMax < self.phi*0.9)  + 0.9*self.phi* (self.RmMax >= 0.9*self.phi)) 
        
    def setGrSink(self): #Rg + C allocation = Grsink = growth * (1/GrowEff)
        self.Gr = self.k_growth  * self.mesh.length * (-self.mesh.length + self.lengthThr) * (self.phi - self.phiThrGrowth) * (self.mesh.length < self.lengthThr)* (self.phi > self.phiThrGrowth) #adapted from Daudet volumic growth for next time step. To send to CPlantBox
        #extra length therm: then logarithmic instead of plateau
        self.GrSink =self.Gr * self.rhoC * (1/self.GrEff)#growth with C content
        
        
    def updateSourceSink(self): #usefull? will have to create a mesh at each timestep.
        self.setAnSource()
        self.setD()
        self.setRmSink()
        self.setGrSink()
        return 0.