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
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D

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
        self.new2organID ={} #for grwth rate evaluation
        self.organ2oldNodeID ={} #for growth rate evaluation
        self.oldNode2organID ={} #for growth rate evaluation
        self.mesh = [] #to store grid
        self.phi = [] # to store solution
        self.faces = []
        self.osmoCoeff = 10. #so for correspond to R*T * mean_phi
        self.intCoeff = 1
        self.outFlow = 1.
        self.radConductivity = 100
        self.satisfaction = 0
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.k1 = 10
        self.k2 = 10
        self.Gr = 0.05 #C-limited growth during time step
        self.Growth = 0.05 #C-limited growth during time step
        self.GrEff = 0.5 #Growth efficiency ( = growth per unit of phi used)
        self.k_growth = 1000
        self.phiThrMax = 1
        self.phiThrMin = 0.001
        self.rhoC = 1 # C content per unit of volume
        self.RmMax = 0.
        self.Px = 0. #xylem water potential
        self.cellsID = 0.
        self.rxCells= np.array([1])
        self.rxThrMax = -200 #cm
        self.rxThrMin = -800 #cm
        self.orgLength = np.array([])
        self.orgLengthThr = np.array([])
        self.orgGr = np.array([])
        self.phiFactor2= -1
        self.maxGrowth = np.array([])

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
        subTypes = self.get_subtypes()
        nOld = nodes[0]
        radiiCells = np.array([])
        radiiVertices = np.array([])
        length = np.array([])
        tempOrgGr = np.array([])
        tempOrgLength = np.array([])    
        
        
        
        for segnum, seg in enumerate(segs):
            radii_org = np.array([])
            length_org = np.array([])
            base = 0.
            seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
            if(segnum == 0):
                seg_b_root = seg
                
            for j, n in enumerate(seg):
                oldNodeID = np.where(np.all(nodes == n, axis=1))[0][0] #find index of n in the node matrix == node ID in the MappedPlant object
                if (oldNodeID > 0):
                    if(j == 0):
                        ot = int(organTypes[oldNodeID])
                        st = int(subTypes[oldNodeID])
                    else:
                        ot = int(organTypes[oldNodeID - 1]) 
                        st = int(subTypes[oldNodeID -1])
                if(oldNodeID == 0 and segnum >0):
                    ot = 3
                    st = 1
                if(oldNodeID == 0 and segnum ==0):
                    ot = 2
                    st = 1
                    
                    
                self.new2organID[newNodeID] = segnum    
                self.new2oldNodeID[newNodeID] =  oldNodeID
                temp = self.old2newNodeID.get(oldNodeID,np.array([], dtype=np.int64))
                self.old2newNodeID[oldNodeID] = np.append(temp, np.array([newNodeID], dtype=np.int64))
                normal = (nOld - n)/norm(nOld - n)
                
                if j > 0: #not first node of organ
                    radiiCells = np.concatenate((radiiCells, [self.rs.radii[oldNodeID - 1]]))
                    
                    self.orgLengthThr = np.concatenate((self.orgLengthThr, [self.rs.organParam[ot][st + 2*(ot==4)].getParameter('lmax')]))
                    self.maxGrowth = np.concatenate((self.maxGrowth, [self.rs.organParam[ot][st + 2*(ot==4)].getParameter('r')]))
                    
                    if(not isinstance(self.Px, float)):
                        self.rxCells = np.concatenate((self.rxCells[self.rxCells < 1], [np.minimum(self.rxThrMax, np.mean((self.Px[oldNodeID], self.Px[self.new2oldNodeID[newNodeID - 1]])))]))
                        
                    else:
                        self.rxCells = -800
                    if j == 1:
                        radiiVertices = np.concatenate((radiiVertices, [self.rs.radii[oldNodeID - 1],self.rs.radii[oldNodeID - 1]])) 
                    else:
                        radiiVertices = np.concatenate((radiiVertices, [self.rs.radii[oldNodeID - 1]]))
                    
                    length = np.concatenate( (length, [norm(nOld - n)]))           
                    length_org =  np.concatenate( (length_org, [norm(nOld - n)]))
                                        
                    if cells.size == 0 : #first organ of the plant
                        #cells =  np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1])))
                        cells =  np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))) #nodes should be ordered. so seg n goes between nodes j-1 and j
                        facesNorm[newNodeID] = normal
                        facesNorm[newNodeID - 1] = -normal #replace value for all node except for the one at the tip of organs
                        
                    elif(ot <3):
                        #cells = np.hstack((cells, np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))).reshape(3,1)))
                        cells = np.hstack((cells, np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))).reshape(3,1)))
                        facesNorm[newNodeID - 1] = -normal #replace value for all node except for the one at the tip of organs
                        facesNorm[newNodeID] = normal
                        
                    else: #stem of leaf
                        #cells = np.hstack((cells, np.array((np.array([newNodeID]),np.array([newNodeID -1]),np.array([-1]))).reshape(3,1)))
                        cells = np.hstack((cells, np.array((np.array([newNodeID - 1]),np.array([newNodeID]),np.array([-1]))).reshape(3,1)))
                        facesNorm[newNodeID] = normal
                elif segnum > 0:      
                    if oldNodeID == 0: #if node 0 (root collar): only case where seg linked with to the x node of parent segment 
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[0]] == self.new2oldNodeID[newNodeID])[0][0]
                        facesNorm[newNodeID] = facesNorm[0]
                        
                    else: #seg linked with to the y node of parent segment
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[1]] == self.new2oldNodeID[newNodeID])[0][0]
                        if(ot == 4):
                            facesNorm[newNodeID] = facesNorm[cells[0,seg_parent]]                       
                    
                    cells[2,seg_parent] = newNodeID #add node to parent seg to link with child branch
                    #cellsbis[2,seg_parent] = newNodeID #add node to parent seg to link with child branch
                
                newNodeID += 1
                nOld = n    
                
            tempOrgLength = np.append(tempOrgLength, sum(length_org))
            nodes_x_coord = np.array([xi[0] for xi in seg], dtype=np.float64) 
            nodes_y_coord = np.array([xi[1] for xi in seg], dtype=np.float64) 
            nodes_z_coord = np.array([xi[2] for xi in seg], dtype=np.float64) 
            nodes_coord =  np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored 
            if vertices.size ==0:
                vertices = nodes_coord
            else:
                vertices = np.hstack((vertices, nodes_coord))
                
        self.orgLength = np.array([tempOrgLength[self.new2organID[xi]] for xi in cells[1]])
        
        self.orgGr = np.full(segnum,-1)
        self.faces = np.array((np.arange(0, vertices[0].size),), dtype=np.int64) 
        #cells = MA.masked_values(cells, -1)
        cells = MA.masked_values(cells, -1)
        self.mesh = Mesh1Dmod(facesNorm=facesNorm, radiiVertices = radiiVertices, 
            radiiCells = radiiCells,length = length, vertexCoords=vertices, faceVertexIDs=self.faces, cellFaceIDs=cells)
                
        self.intCoeff = FaceVariable(name = 'diffusion coefficient', mesh = self.mesh, value = self.osmoCoeff)
        self.Source = CellVariable(name = 'An', mesh = self.mesh, value = 0.)
        self.phi = CellVariable(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        self.resetValues()
        
    def _setOutflow(self): #or set dirichlet of 0 on outer face? or fixed gradient?
        nodeYRoot = self.get_segments_index(2) + 1 #ynode_index of root
        nodeYRoot_newID = np.array([self.old2newNodeID[xi][0] for xi in nodeYRoot]).flatten() #new index
        rootCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in nodeYRoot_newID]) 
        outFlowMax = rootCells*self.phi*self.radConductivity * (self.mesh.radiiCells**2)*self.mesh.length* self.satisfaction #reset for security 
        self.outFlow = outFlowMax * (self.phi - outFlowMax > 0) * (self.phi > 0) + self.phi * (self.phi - outFlowMax <= 0)

    def setintCoeff(self): #use to get right value of k/mu* RT
        self.intCoeff.setValue(self.osmoCoeff) 

    def setAnSource(self): #sourceterm in leaf segment
        self.Source.setValue(0.)
        # find seg index with leaf, should be ordered like the An vector
        nodeYLeaf = self.get_segments_index(4) + 1 #ynode_index of leaves
        nodeYLeaf_newID = np.array([self.old2newNodeID[xi] for xi in nodeYLeaf]).flatten() #new index
        surfaceSide = 2*np.pi*self.mesh.radiiCells*self.mesh.length
        for i, an in enumerate(self.An):#An is per surface area and per second!
                sucrose = an * 1e9 *60*60 *24/ 12 #nano_mol sucrose m-2 d-1
                self.Source.constrain(sucrose* surfaceSide /self.mesh.cellVolumes, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #0 node = y node for stem and leaves
                #self.Source.constrain(an*self.mesh.length/self.mesh.cellVolumes, where= self.mesh.cellFaceIDs[0] == nodeYLeaf_newID[i])
      
    def setRmSink(self): # Daudet 2002: Rm = (k1 + k2 * C)Sr. 
        self.RmMax = self.rhoC * self.mesh.cellVolumes * (self.phi * self.k1 + self.k2 )
        self.Rm = (self.phi > 0)*( self.RmMax * (self.phi - self.RmMax >0 )+ (self.phi )* (self.phi - self.RmMax <= 0 ) )#*(self.phi >0)#*(self.phi> self.phiThrGrowth))
       
        self.CSat = (self.Rm == self.RmMax) 
       
    def setGrSink(self): #Rg + C allocation = Grsink = growth * (1/GrowEff)
    
        #length, phi and rx :  factor from 0 to 1
#set seglen/orglen * maxgrowth rate * rxf * lenfact * phifact        
        
        rxFactor = ((self.rxCells - self.rxThrMin)/(self.rxThrMax - self.rxThrMin)) *(self.rxCells > self.rxThrMin)
        lengthFactor = (self.orgLength*(-self.orgLength + self.orgLengthThr)/((self.orgLengthThr**2)/2)) * (self.orgLength < self.orgLengthThr)
        #self.phiFactor1 = (self.phi > self.phiThrMin)*((np.minimum(self.phi, self.phiThrMax) - self.phiThrMin)/(self.phiThrMax - self.phiThrMin))
    
        self.Gr =  ((self.maxGrowth*np.pi * self.mesh.radiiCells**2)*self.rhoC /self.GrEff)* (self.mesh.length/self.orgLength) *  rxFactor * lengthFactor * self.CSat #* abs(self.phiFactor1)
        
        self.phiFactor2 = ((self.phi - self.Gr) > self.phiThrMin)  
        self.GrSink =  (self.phi > self.phiThrMin)*((self.Gr  * self.phiFactor2 + (self.phi - self.phiThrMin) * (~self.phiFactor2)))
        
        self.Growth = self.GrSink / self.rhoC * self.GrEff
        self.satisfaction = self.CSat * (self.GrSink == self.Gr)
        
    def resetValues(self):
        self.setAnSource()
        self.setintCoeff()
        self.setRmSink()
        self.setGrSink()
        self._setOutflow()