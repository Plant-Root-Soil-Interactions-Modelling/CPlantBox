import timeit
import math
import fipy as fp

from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
from fipy.tools import numerix 
import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA
from datetime import datetime, timedelta

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
        self.mesh = [] #to store grid
        self.phi = [] # to store solution
        self.outFlow = 1.
        self.radConductivity = 100
        self.satisfaction = 0
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.k1 = 1e-3 /(24*60*60)/12/12/10000 #mol sucroses cm2 s. [/12]: mol C to mol sucrose, [/12]: mol/g C, [10000] m2 to cm2
        self.k2 = 1e-3 #check from LEroux 2001 from range in table IV
        self.Gr = 0.05 #C-limited growth during time step
        self.Growth = 0.05 #C-limited growth during time step
        self.GrEff = 0.75 #Growth efficiency ( = growth per unit of phi used), value from LEroux 2001 from range in table IV
        #self.k_growth = 1000
        self.phiThrMax = 0.002
        self.phiThrMin = 1e-20
        self.rhoC = 1# C content per unit of volume
        self.RmMax = 0.
        self.Px = 0. #xylem water potential
        self.rxCells= np.array([1])
        self.rxThrMax = -200 #cm
        self.rxThrMin = -16300 #cm
        self.orgLength = np.array([])
        self.orgLengthThr = np.array([])
        self.orgGr = np.array([])
        self.phiFactor2= -1
        self.lengthFactor = -1
        self.maxGrowth = np.array([])
        #change again according to data of xiaoran
        self.phloemConductivity = 1.118e-12*1e4  #k in cm2 . see cpp PIA_munch file  1e-11#
        self.phloemViscosity = 1e-3
        self.TairK = 293
        self.VolFractSucrose = np.array([])
        self.intCoeff = 100.
        self.intCoeff1 = 100.
        self.rootCells = 0.
        self.nodeYLeaf_newID = -1.
        self.new2oldNodeID = {} #map to get MappedPlant node-ID from grid node-ID
        self.old2newNodeID = {} #map to get grid node-ID from MappedPlant node-ID 
        self.newCell2organID ={} #for grwth rate evaluation
        self.organ2oldNodeID ={} #for growth rate evaluation
        self.organ2newCellID ={} #for growth rate evaluation
        self.oldNode2organID ={} #for growth rate evaluation
        self.orgID=np.array([])
        self.VolumeSEG = 1
        self.StructC = 0

    def mp2mesh(self, segs, TairK = 293): #TODO: adapt to make more time efficient. do it in c++?
        """ Converts a mappedPlant into a network of 1D grids            
            outputs: 
            creates mesh and stores it in self.mesh
            fills self.phi => allows for visualization of the mesh with vtk"""
        self.TairK = TairK
        nodes = self.get_nodes() 
        
        #cells = [[bottom 1], [top 1], [bottom 2], [top 2]]
        #bottom 1 = x_node: link to segment bellow or parent seg bellow
        #top 1: y_node: link to segment above
        #bottom 2: x_node: link to child bellow or parent cell above
        #top 2: y_node: link to child above 
        #only bottom 1 and top 1 are used to define cellCenter
        
        for segnum, seg in enumerate(segs):#go through each organ
        
            seg = np.array(list(map(lambda x: np.array(x), seg))) # converts the list of Vector3d to a numpy array. shape: (num_nodes, 3)
            
            if segnum == 0: #basal root
                bot1 = np.array([s for s in range(len(seg)-1)], dtype=np.int64) #nodes put back in oder after growth?
                top1 = np.array([s for s in range(1,len(seg))]) #nodes put back in oder after growth?
                bot2 = np.full(len(seg)-1, -1)
                top2ID = np.array([])
                
            elif np.where(np.all(nodes == seg[0], axis=1))[0][0] == 0: #main stem
                start = max(top1) + 1
                #max val of nde instead of len seg
                bot1 = np.concatenate((bot1, np.array([0]), np.array([s + start for s in range(0,len(seg)-2)], dtype=np.int64))) #nodes put back in oder after growth?
                top1 =  np.concatenate((top1,np.array([s + start for s in range(0,len(seg)-1)]))) #nodes put back in oder after growth?
                bot2 = np.concatenate((bot2, np.full(len(seg)-1, -1)))
                
            else: #lateral 
                start = max(top1) +1
                bot1 =  np.concatenate((bot1, np.array([start]),np.array([s + start for s in range(2,len(seg))], dtype=np.int64))) #nodes put back in oder after growth?
                top1 = np.concatenate((top1,np.array([s + start for s in range(2,len(seg)+1)]))) #nodes put back in oder after growth?
                top2ID = np.concatenate((top2ID, [start]))
                bot2 = np.concatenate((bot2,[start + 1], np.full(len(seg)-2, -1)))
        
        oldbot1 = np.concatenate([np.array(list(map(lambda x: np.where(np.all(nodes == np.array(x), axis=1))[0][0], segs[i])))[:-1] for i in range(len(segs)) ])
        oldtop1 = np.concatenate([np.array(list(map(lambda x: np.where(np.all(nodes == np.array(x), axis=1))[0][0], segs[i])))[1:] for i in range(len(segs)) ])
        oldbot2 = oldbot1[np.where(bot2>=0)[0]] 
        newIDS = np.concatenate((bot1,top1, bot2[np.where(bot2>=0)[0]]))
        oldIDS  = np.concatenate((oldbot1,oldtop1, oldbot2))
        self.new2oldNodeID = dict(zip(newIDS,oldIDS))
        oldIDSu = np.array(list(set(oldIDS)))
        newIDS = [np.array(list(set(newIDS[np.where(oldIDS == ui)[0]]))) for ui in oldIDSu]
        self.old2newNodeID = dict(zip(oldIDSu,newIDS ))
        
        self.orgID = np.concatenate([np.full(len(seg)-1, segnum) for segnum, seg in enumerate(segs)])# organ id for each cell
        faces = np.array((np.arange(0,  max(top1) + 1),), dtype=np.int64) #index of verticies on each face
        
        nodes_x_coord = [nodes[self.new2oldNodeID[ni]][0] for ni in faces[0]]     
        nodes_y_coord = [nodes[self.new2oldNodeID[ni]][1] for ni in faces[0]]     
        nodes_z_coord = [nodes[self.new2oldNodeID[ni]][2] for ni in faces[0]]     
        vertices = np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored 
        
        top2 = np.full(len(bot1), -1)
        top2Idx = np.where([sum(oldbot1== t1) > 1 for t1 in oldtop1])[0]#np.where(oldtop1== b1)#np.concatenate((top2Idx,np.where()))
        top2[top2Idx] = top2ID
        bot2[top2Idx + 1] =  top2ID + 1
        cells =  np.vstack((bot1,top1, bot2, top2))
        
        if(not isinstance(self.Px, float)):#get mean rx of the segment (used by GrSink)
            self.rxCells = np.mean([[self.Px[self.new2oldNodeID[xi]] for xi in cells[0]], 
                                   [self.Px[self.new2oldNodeID[xi]] for xi in cells[1]]], axis=0)
            #rxCell >= rxThrMax -> rxFactor = 1
        else:
            self.rxCells = np.minimum(self.rxThrMax,self.Px)
            
        length = np.array([norm(nodes[self.new2oldNodeID[cells[0][i]]]-nodes[self.new2oldNodeID[cells[1][i]]]) for i in range(len(cells[0]))])
        
        
        cells = MA.masked_values(cells, -1)
        radiiVertices = np.array([self.rs.radii[max((self.new2oldNodeID[xi]-1,0))] for xi in faces[0]]) #radius for each node, seed node take radius of 1st segment  of root  
        radiiCells = np.array([self.rs.radii[self.new2oldNodeID[xi]-1] for xi in cells[1]])
        
        self.mesh = Mesh1Dmod( radiiVertices = radiiVertices*1e-2, 
            radiiCells = radiiCells*1e-2,length = length, vertexCoords=vertices, 
            faceVertexIDs=faces, cellFaceIDs=cells) #1e-2 to go from radius of segment to radius of phloem
        #take out radiiCells and length: compute durectly in Mesh1Dmod
        
                                                
        #store max growth rate and max organ length for computaiton of the CW-limited growth
        self.orgLength = np.array([sum(length[np.where(self.orgID == oid )[0]]) for oid in self.orgID])
        organTypes = self.get_organ_types()
        subTypes = self.get_subtypes()
        self.orgLengthThr = np.array( [self.rs.organParam[organTypes[self.new2oldNodeID[xi] - 1]][subTypes[self.new2oldNodeID[xi] -1] + 2*(organTypes[self.new2oldNodeID[xi] - 1]==4)].getParameter('lmax') for xi in cells[1]])
        self.maxGrowth = np.array( [self.rs.organParam[organTypes[self.new2oldNodeID[xi] - 1]][subTypes[self.new2oldNodeID[xi] -1]+ 2*(organTypes[self.new2oldNodeID[xi] - 1]==4)].getParameter('r')/(24*60*60) for xi in cells[1]]) #maxgrwth rate from cm/d to cm/s
                      
        
        self.phi = CellVariablemod(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        self.Rm = CellVariablemod(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        
    
    def _setOutflow(self): #only outflow when all segments have rm gr satisfaction?
        nodeYRoot = self.get_segments_index(2) + 1 #ynode_index of root
        nodeYRoot_newID = np.array([self.old2newNodeID[xi][0] for xi in nodeYRoot]).flatten() #new index
        self.rootCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in nodeYRoot_newID]) 
        #outFlowMax = self.rootCells*self.phi*self.radConductivity * ((self.mesh.radiiCells)**2)*self.mesh.length#* self.satisfaction #reset for security 
        #here rad of sefg or of phloem?
        self.outFlow =(self.phi > 1e-10) *(self.phi - 1e-10)*self.rootCells*self.radConductivity * ((self.mesh.radiiCells)**2)*self.mesh.length#*self.phi* ((np.minimum(self.phi, self.phiThrMax)- self.phiThrMin)/(self.phiThrMax- self.phiThrMin))#*((self.phi -  np.min(self.phi, self.phiThrMin))/(0.003- self.phiThrMin))#(outFlowMax * (self.phi - outFlowMax > self.phiThrMin) + self.phi * (self.phi - outFlowMax <= self.phiThrMin))* (self.phi > self.phiThrMin)
        #print('ouflow ',outFlowMax, self.outFlow)

    def _setintCoeff(self): #use to get right value of k/mu* RT
        # R=8.314⋅m3⋅Pa⋅K−1mol−1 = 8.314⋅e6 mL⋅Pa⋅K−1mol−1
        #for fipy, use resistivity (conductivity) and not resistance (conductance). See: https://www.ctcms.nist.gov/fipy/documentation/numerical/discret.html#diffusion-term-nabla-cdot-left-gamma-1-nabla-phi-right
        waterViscosity = 2.414e-5 * (10**(247.8/(self.TairK - 140))) #in [Pa.s = N-s/m^2] https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
        sucroseDensity = 1.59 #g/cm³, https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose#section=Density
        sucroseMolarMass = 342.3 #g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose
        sucroseMolarVolume = sucroseMolarMass/sucroseDensity #cm³/mol like the one used by lacointe 2019
        self.VolFractSucrose = self.phi.faceValue* sucroseMolarVolume #  #[mol/cm³] / [cm³/mol] = cm³/cm³ * sucroseMolarVolume
           
        #valid for VolFractSucrose <= 0.65 <-> [sucrose] <= 0.003019281 mol/cm3
        self.phloemViscosity = waterViscosity * np.exp((4.68 * 0.956 * self.VolFractSucrose)/(1 - 0.956 * self.VolFractSucrose)) #in [Pa.s = N-s/m^2]
        R = self.Param['R'] * 1e6 # Pa * mL /(K mol)    
        self.intCoeff =self.phloemConductivity/self.phloemViscosity * R * self.TairK 
        # flow = (k/mu) * RT * d(C)/dx = [cm2]/[Pa s] * [Pa * mL /(K mol)] * K * d[mol/mL] 
        # = [cm2 Pa-1 s-1] * d[Pa ] = cm2 s-1 (?)

    def _setAnSource(self): #sourceterm in leaf segment
        self.Source = CellVariablemod(name = 'An', mesh = self.mesh, value = 0.)
        # find seg index with leaf, should be ordered like the An vector
        nodeYLeaf = self.get_segments_index(4) + 1 #ynode_index of leaves
        nodeYLeaf_newID = np.array([self.old2newNodeID[xi] for xi in nodeYLeaf]).flatten() #new index
        surfaceSide = (2*np.pi*((self.mesh.radiiCells)/100)*(self.mesh.length/100)) #m2. 1e2 : to go from phloem radius to leaf radius
        for i, an in enumerate(self.An):
                sucrose = 1.8e-12#an / 12 *0.00001#mol sucrose m-2 s-1, redue by 0.001 just to get ok sucrose increase 1.8e-10 #mol/s-1#
                self.Source.constrain(sucrose, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #(sucrose* surfaceSide /self.mesh.cellVolumes, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #0 node = y node for stem
        #print('source ',self.An,surfaceSide , self.Source)
        
        
        
    def _setRmSink(self): # Daudet 2002: Rm = (k1 + k2 * C)Sr. 
        self.volumesSEG = np.array([(np.pi * (self.mesh.radiiCells[xi]*1e2)**2) *self.mesh.length[xi] for xi in range(self.mesh.numberOfCells) ])
        #1e2 to go from phloem radius to plant seg radius
        self.RmMax = (self.rhoC * self.volumesSEG * ( self.k2))#(self.phi * self.k1)# + self.k2 )
        self.StructC = CellVariablemod(value = self.rhoC * self.volumesSEG , mesh = self.mesh)
        self.Rm = CellVariablemod(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        #self.Rm = self.RmMax * (self.RmMax/self.phi <= 0.8) + 0.8 *  (self.RmMax/self.phi >= 0.8)
        #(self.phi > self.phiThrMin)* self.rhoC * volumesSEG * (self.phi * self.k1 + self.k2 )
        #*((self.phi - np.min(self.phi, self.phiThrMin))/(0.003- self.phiThrMin))
        #(self.phi > self.phiThrMin)*( self.RmMax * (self.phi - self.RmMax >self.phiThrMin )+ (self.phi - self.phiThrMin)* (self.phi - self.RmMax <= self.phiThrMin ) )#*(self.phi >0)#*(self.phi> self.phiThrGrowth))
        #self.CSat = (self.Rm == self.RmMax) 
        #print('rm sink ',self.RmMax, self.Rm)
       
    def _setGrSink(self): #Rg + C allocation = Grsink = growth * (1/GrowEff)
        ##only growth when all segments have rm satisfaction?
        #length and rx :  factor from 0 to 1              
        volumesSEG = np.array([(np.pi * (self.mesh.radiiCells[xi]*1e2)**2) *self.mesh.length[xi] for xi in range(self.mesh.numberOfCells) ])
        rxFactor = ((self.rxCells - self.rxThrMin)/(self.rxThrMax - self.rxThrMin)) *(self.rxCells > self.rxThrMin)
        self.lengthFactor = (self.orgLength*(-self.orgLength + self.orgLengthThr)/((self.orgLengthThr/2)**2)) * (self.orgLength < self.orgLengthThr)
        maxCconcentationNeeded = (((self.maxGrowth*np.pi * (self.mesh.radiiCells*1e2)**2)*self.rhoC /self.GrEff))\
        * (self.mesh.length/self.orgLength)/self.mesh.cellVolumes
        self.Gr =  maxCconcentationNeeded/12 *  rxFactor * self.lengthFactor* (1/self.GrEff)# * self.CSat 
        
        #self.phiFactor2 = ((self.phi - self.Gr) > self.phiThrMin)  
        #self.GrSink = self.GrSink * (self.GrSink/self.phi <= 0.8) + 0.8 *  (self.GrSink/self.phi >= 0.8)
        #self.Gr * ((np.minimum(self.phi, self.phiThrMax)- self.phiThrMin)/(self.phiThrMax- self.phiThrMin))#*self.phi #(self.phi > self.phiThrMin)*(self.Gr  * self.phiFactor2 + (self.phi - self.phiThrMin) * (~self.phiFactor2))
        
        #self.Growth0 = self.GrSink* self.mesh.cellVolumes *12/ self.rhoC * self.GrEff  #get increase in cm³
        #self.Growth = self.Growth0 /(np.pi * (self.mesh.radiiCells*1e2)**2)#get length increase (cm)
        #self.satisfaction = self.CSat * (self.GrSink == self.Gr)
        #print('gr sink ',maxCconcentationNeeded, self.lengthFactor ,rxFactor)

    @property        
    def organGr(self):#per cell not node
    #ATT: changed to dummy value for testing
        self.GrSink = self.phi * self.Rm #- self.RmMax
        #print('Gr ',self.totalSinkFactor,self.RmMax,  self.GrSink, self.orgID, np.array([max(0,xi) for xi in self.GrSink]))
        self.GrSink = np.array([max(0,xi) for xi in self.GrSink])
        self.Growth0 = self.GrSink* self.mesh.cellVolumes *12/ self.rhoC * self.GrEff  #get increase in cm³
        self.Growth = self.Growth0 /(np.pi * (self.mesh.radiiCells*1e2)**2)#get length increase (cm)
        self.orgGr = np.array([0.01 for xi in range(max(self.orgID )+1)]) #np.array([sum(self.Growth[np.where(self.orgID == xi)[0]])*1e25 for xi in range(max(self.orgID )+1)]) 
        
        return self.orgGr

        
    def resetValues(self):
        self._setAnSource()
        self._setintCoeff()
        self._setRmSink()
        self._setGrSink()
        self._setOutflow()#+ self.outFlow
        
        self.totalSinkFactor = (self.RmMax + self.Gr)/(self.phi+(self.RmMax + self.Gr)*1e3)#)/self.phi )* ((self.RmMax + self.Gr )/self.phi <= 0.8) + 0.8 * ((self.RmMax + self.Gr )/self.phi >= 0.8))*(self.phi>0)
        

def GetTime(seconds):
    sec = timedelta(seconds=int(seconds))
    d = datetime(1,1,1) + sec
    return(f'{d.day-1}d {d.hour}hr {d.minute}mn {d.second}s')        
        
 