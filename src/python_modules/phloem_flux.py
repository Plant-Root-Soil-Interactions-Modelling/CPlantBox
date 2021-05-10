import timeit
import math
import fipy as fp

from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
from fipy.tools import numerix 
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
from CellVariablemod import CellVariablemod
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
        self.newCell2organID ={} #for grwth rate evaluation
        self.organ2oldNodeID ={} #for growth rate evaluation
        self.organ2newCellID ={} #for growth rate evaluation
        self.oldNode2organID ={} #for growth rate evaluation
        self.mesh = [] #to store grid
        self.meshes = [] #to store grid
        self.phi = [] # to store solution
        self.faces = []
        self.outFlow = 1.
        self.radConductivity = 100
        self.satisfaction = 0
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.k1 = 0.001
        self.k2 = 0.004
        self.Gr = 0.05 #C-limited growth during time step
        self.Growth = 0.05 #C-limited growth during time step
        self.GrEff = 0.99 #Growth efficiency ( = growth per unit of phi used)
        #self.k_growth = 1000
        #self.phiThrMax = 1
        self.phiThrMin = 0.01
        self.rhoC = 1. # C content per unit of volume
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
        self.lengthFactor = -1
        self.maxGrowth = np.array([])
        self.phloemConductivity = 1/2.5 /(60*60)#*1e-6#cm3 s-1 Pa-1 cf daudet 2002
        self.phloemViscosity = 1e-3
        self.TairK = 293
        self.VolFractSucrose = np.array([])
        self.intCoeff = 100.
        self.intCoeff1 = 100.
        self.rootCells = 0.
        self.nodeYLeaf_newID = -1.

    def mp2mesh(self, segs, TairK = 293): #create an old2newNodeID map?
        """ Converts a mappedPlant into a network of 1D grids            
            outputs: 
            creates mesh and stores it in self.mesh
            fills self.phi => allows for visualization of the mesh with vtk"""
        self.TairK = TairK
        vertices = np.array([]) #containes x,y,z coords of each vertex
        cells = np.array([[],[],[] ,[]], dtype=np.int64) #contains ID of the faces on each cell. here 1D => 1face == 1 vertex
        #cells = [[bottom 1], [top 1], [bottom 2], [top 2]]
        #bottom 1 = x_node: link to segment bellow or parent seg bellow
        #top 1: y_node: link to segment above
        #bottom 2: x_node: link to child bellow or parent cell above
        #top 2: y_node: link to child above 
        #only bottom 1 and top 1 are used to define cellCenter
        nodes = self.get_nodes() 
        newNodeID = 0
        organTypes = self.get_organ_types()
        subTypes = self.get_subtypes()
        radiiCells = np.array([])
        radiiVertices = np.array([])
        length = np.array([])
        tempOrgLength = np.array([])    
        newNodeID_prev = -1
        nodes_x_coord = np.array([], dtype=np.float64)#np.array([xi[0] for xi in seg], dtype=np.float64) 
        nodes_y_coord = np.array([], dtype=np.float64)# np.array([xi[1] for xi in seg], dtype=np.float64) 
        nodes_z_coord = np.array([], dtype=np.float64)# np.array([xi[2] for xi in seg], dtype=np.float64) 
        self.tipmeshx = 0.
        self.tipmeshy = 0.
        self.tipplant = 0.
        
        for segnum, seg in enumerate(segs):
            length_org = np.array([])
            seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
                
            for j, n in enumerate(seg):
                oldNodeID = np.where(np.all(nodes == n, axis=1))[0][0] 
                                    
                if(oldNodeID ==0 and segnum != 0): #do not do it for first node of stem (newnodeId = )
                    newNodeID_prev = 0
                    nOld = n
                    continue
                self.new2oldNodeID[newNodeID] =  oldNodeID
                temp = self.old2newNodeID.get(oldNodeID,np.array([], dtype=np.int64))
                self.old2newNodeID[oldNodeID] = np.append(temp, np.array([newNodeID], dtype=np.int64))
                radiiVertices = np.concatenate((radiiVertices, [(0.2/np.pi)**0.5]))#[self.rs.radii[oldNodeID - 1]]))
                
                if j > 0: #not first node of organ
                    ot = int(organTypes[oldNodeID - 1]) 
                    st = int(subTypes[oldNodeID -1])
                    radiiCells = np.concatenate((radiiCells, [(0.2/np.pi)**0.5]))
                    
                    
                    self.newCell2organID[newNodeID] = segnum    
                    temp = self.organ2newCellID.get(segnum,np.array([], dtype=np.int64))
                    self.organ2newCellID[segnum ] = np.append(temp, np.array([newNodeID], dtype=np.int64))
                    
                    self.orgLengthThr = np.concatenate((self.orgLengthThr, [self.rs.organParam[ot][st + 2*(ot==4)].getParameter('lmax')]))
                    self.maxGrowth = np.concatenate((self.maxGrowth, [self.rs.organParam[ot][st + 2*(ot==4)].getParameter('r')]))
                    
                    if(not isinstance(self.Px, float)):
                        self.rxCells = np.concatenate((self.rxCells[self.rxCells < 1], [np.minimum(self.rxThrMax, np.mean((self.Px[oldNodeID], self.Px[self.new2oldNodeID[newNodeID_prev]])))]))
                        
                    
                    length = np.concatenate( (length, [norm(nOld - n)]))           
                    length_org =  np.concatenate( (length_org, [norm(nOld - n)]))
                    if (segnum == 0 or cells[1][len(cells[1]) - 1] != newNodeID): #not first cell of lateral stem, lat root or leaf
                        cells = np.hstack((cells, np.array((np.array([newNodeID_prev]),np.array([newNodeID]),np.array([-1]),np.array([-1]))).reshape(4,1)))  
                    
                elif segnum > 0:      
                    ot = int(organTypes[oldNodeID ]) 
                    seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[1]] == self.new2oldNodeID[newNodeID])[0][0]
                    if seg_parent >= 0: #not for base stem
                        cells[3,seg_parent] = newNodeID 
                        newNodeID_prev = newNodeID
                        newNodeID += 1
                        
                        nodes_x_coord = np.append( nodes_x_coord, n[0])#np.array([xi[0] for xi in seg], dtype=np.float64) 
                        nodes_y_coord = np.append( nodes_y_coord, n[1])# np.array([xi[1] for xi in seg], dtype=np.float64) 
                        nodes_z_coord = np.append( nodes_z_coord, n[2])# np.array([xi[2] for xi in seg], dtype=np.float64) 
                            
                        self.new2oldNodeID[newNodeID] =  oldNodeID
                        temp = self.old2newNodeID.get(oldNodeID,np.array([], dtype=np.int64))
                        self.old2newNodeID[oldNodeID] = np.append(temp, np.array([newNodeID], dtype=np.int64))
                        radiiVertices = np.concatenate((radiiVertices, [(0.2/np.pi)**0.5]))#[self.rs.radii[oldNodeID - 1]]))
                        
                        seg_parent = np.where([self.new2oldNodeID[xi] for xi in cells[0]] == self.new2oldNodeID[newNodeID])[0][0]
                        cells[2,seg_parent] = newNodeID #add node to parent seg to link with child branch: bot2
                        cells = np.hstack((cells, np.array((np.array([newNodeID_prev]),np.array([newNodeID + 1]),np.array([newNodeID]),np.array([-1]))).reshape(4,1))) 
                nodes_x_coord = np.append( nodes_x_coord, n[0])#np.array([xi[0] for xi in seg], dtype=np.float64) 
                nodes_y_coord = np.append( nodes_y_coord, n[1])# np.array([xi[1] for xi in seg], dtype=np.float64) 
                nodes_z_coord = np.append( nodes_z_coord, n[2])# np.array([xi[2] for xi in seg], dtype=np.float64) 
                newNodeID_prev = newNodeID
                newNodeID += 1
                
                nOld = n    
            
            tempOrgLength = np.append(tempOrgLength, sum(length_org))
                 
        vertices = np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored  
        self.orgGr = np.full(segnum +1, 0.)
        self.orgLength = np.array([tempOrgLength[self.newCell2organID[xi]] for xi in cells[1]])
        print("orrglength ", self.orgLength)
        self.faces = np.array((np.arange(0,  max(cells[1]) + 1),), dtype=np.int64) 
        cells = MA.masked_values(cells, -1)
        self.mesh = Mesh1Dmod( radiiVertices = radiiVertices, 
            radiiCells = radiiCells,length = length, vertexCoords=vertices, 
            faceVertexIDs=self.faces, cellFaceIDs=cells)
        
        self.phi = CellVariablemod(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        
    
    def _setOutflow(self): #only outflow when all segments have rm gr satisfaction?
        nodeYRoot = self.get_segments_index(2) + 1 #ynode_index of root
        nodeYRoot_newID = np.array([self.old2newNodeID[xi][0] for xi in nodeYRoot]).flatten() #new index
        self.rootCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in nodeYRoot_newID]) 
        outFlowMax = self.rootCells*self.phi*self.radConductivity * (self.mesh.radiiCells**2)*self.mesh.length* self.satisfaction #reset for security 
        self.outFlow = outFlowMax * (self.phi - outFlowMax > 0) * (self.phi > 0) + self.phi * (self.phi - outFlowMax <= 0)* (self.phi > 0)
        #print('ouflow ',outFlowMax, self.outFlow)

    def _setintCoeff(self): #use to get right value of k/mu* RT
        
        waterViscosity = 2.414e-5 * (10**(247.8/(self.TairK - 140))) #in [Pa.s = N-s/m^2] https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
        sucroseDensity = 1.59 #g/cm³, https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose#section=Density
        sucroseMolarMass = 342.3 #g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose
        sucroseMolarVolume = sucroseMolarMass/sucroseDensity #cm³/mol change to make variable volume?
        self.VolFractSucrose = self.phi.faceValue* sucroseMolarVolume#np.minimum(np.full(len(self.phi.faceValue),0.65),np.array(self.phi.faceValue* sucroseMolarVolume) ) #[mol/cm³] / [cm³/mol] = cm³/cm³ * sucroseMolarVolume
           
        #print(len(np.full(len(self.phi.faceValue),0.65)),len(self.phi.faceValue * sucroseMolarVolume),self.VolFractSucrose)
        self.phloemViscosity = waterViscosity * np.exp((4.68 * 0.956 * self.VolFractSucrose)/(1 - 0.956 * self.VolFractSucrose))
        R = self.Param['R'] * 1e6 # Pa * mL /(K mol)
        self.intCoeff =10#self.phloemConductivity/self.phloemViscosity * R * self.TairK #(self.osmoCoeff)
        #print('intcoeff definition',self.phi,self.phi.faceValue,self.VolFractSucrose, self.phloemViscosity,self.intCoeff)

    def _setAnSource(self): #sourceterm in leaf segment
        self.Source = CellVariable(name = 'An', mesh = self.mesh, value = 0.)
        # find seg index with leaf, should be ordered like the An vector
        nodeYLeaf = self.get_segments_index(4) + 1 #ynode_index of leaves
        nodeYLeaf_newID = np.array([self.old2newNodeID[xi] for xi in nodeYLeaf]).flatten() #new index
        surfaceSide = 2*np.pi*(self.mesh.radiiCells/100)*(self.mesh.length/100) #m2
        #print('co2 assimilation rate ', max(self.An)* 60*60 *24/ 12)
        for i, an in enumerate(self.An):#An is per surface area and per second!
                sucrose = an *1e6* 60*60 *24/ 12 #mol sucrose m-2 d-1
                self.Source.constrain(sucrose* surfaceSide /self.mesh.cellVolumes, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #0 node = y node for stem
        
        
        
        
    def _setRmSink(self): # Daudet 2002: Rm = (k1 + k2 * C)Sr. 
        self.RmMax = self.rhoC * self.mesh.cellVolumes * (self.phi * self.k1 + self.k2 )
        self.Rm = (self.phi > self.phiThrMin)*( self.RmMax * (self.phi - self.RmMax >self.phiThrMin )+ (self.phi )* (self.phi - self.RmMax <= self.phiThrMin ) )#*(self.phi >0)#*(self.phi> self.phiThrGrowth))
        self.CSat = (self.Rm == self.RmMax) 
        #print('rm sink ',self.RmMax, self.Rm)
       
    def _setGrSink(self): #Rg + C allocation = Grsink = growth * (1/GrowEff)
        ##only growth when all segments have rm satisfaction?
        #length and rx :  factor from 0 to 1              
        rxFactor = ((self.rxCells - self.rxThrMin)/(self.rxThrMax - self.rxThrMin)) *(self.rxCells > self.rxThrMin)
        self.lengthFactor = (self.orgLength*(-self.orgLength + self.orgLengthThr)/((self.orgLengthThr**2)/2)) * (self.orgLength < self.orgLengthThr)
        maxCconcentationNeeded = (((self.maxGrowth*np.pi * self.mesh.radiiCells**2)*self.rhoC /self.GrEff))\
        * (self.mesh.length/self.orgLength)/self.mesh.cellVolumes
        self.Gr =  maxCconcentationNeeded/12 *  rxFactor * self.lengthFactor * self.CSat #* abs(self.phiFactor1)
        
        self.phiFactor2 = ((self.phi - self.Gr) > self.phiThrMin)  
        self.GrSink =  (self.phi > self.phiThrMin)*((self.Gr  * self.phiFactor2 + (self.phi - self.phiThrMin) * (~self.phiFactor2)))
        
        self.Growth0 = self.GrSink / self.rhoC * self.GrEff * self.mesh.cellVolumes #get increase in cm³
        self.Growth = self.Growth0 /(np.pi * self.mesh.radiiCells**2)#get length increase (cm)
        self.satisfaction = self.CSat * (self.GrSink == self.Gr)
        #print('gr sink ',maxCconcentationNeeded, self.GrSink)

    @property        
    def organGr(self):#per cell not node
        for cellNum, cellID in enumerate(self.mesh.cellFaceIDs[1]):
            self.orgGr[self.newCell2organID[cellID]] += self.Growth[cellNum]
        return self.orgGr
    
    
        
    def resetValues(self):
        self._setAnSource()
        self._setintCoeff()
        self._setRmSink()
        self._setGrSink()
        self._setOutflow()
        

        
        
 