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
        self.plant = rs
        self.mesh = [] #to store grid
        self.phi = [] # to store solution
        self.outFlow = 1.
        self.radConductivity = 1#00*10197.2*(60*60) #MPa hr/ ml to cm/(s ml)#0.2851281
        self.satisfaction = 0
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.k1 = 2e-4# 4e-10#0.04/(60*60)*1e-2 # from xiaoran Zhou h to s, 1e-3 /(24*60*60)/12/12/10000 #mol sucroses cm2 s. [/12]: mol C to mol sucrose, [/12]: mol/g C, [10000] m2 to cm2
        self.k2 = 4e-4#4e-10#3.3e-9#1e-6#XZ 1e-3 #check from LEroux 2001 from range in table IV
        self.Gr = 0.05 #C-limited growth during time step
        self.Growth = 0.05 #C-limited growth during time step
        self.GrEff = 0.75 #Growth efficiency ( = growth per unit of phi used), value from LEroux 2001 from range in table IV
        #self.k_growth = 1000
        self.phiThrMax = 0.002
        self.phiThrMin = 1e-20
        self.rhoC = 0.07*0.4*0.4/12 
        #straw bulk density: mean(24, 111) kg m3 => 0.07 g cm-3 #http://biomasslogistics.org/Publications/22lam.pdf
        # mol C/cm3 fresh plant: [0.07 g/cm3freshmatter] [0.4 gC/gdryplant (cf elemental analysis ackelysimeter)]*[0.4 gDryplant/gwetplant (cf erntetabelle ackelysimeter)]\
        #*[1/12 molC/gC]content (mol) per unit of volume
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
        self.phloemConductivity = 1.118e-12*1e4 #1.118e-12*1e4  #k in cm2 . see cpp PIA_munch file  1e-11#
        self.phloemViscosity = 1.7e-3 #Pa s-1 , cf Zhou 2020
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
        self.radiiCells = 0
        self.conductivityLeaves=10
        self.phiLim = 0
        self.lim = 0.1e-3/12 #mol Suc ml-3, see zhou et al. 2019
        self.VarMu = True
       
        
    def mp2mesh(self, TairK = 293, VariableMu = True): #TODO: adapt to make more time efficient. do it in c++?
        """ Converts a mappedPlant into a network of 1D grids            
            outputs: 
            creates mesh and stores it in self.mesh
            fills self.phi => allows for visualization of the mesh with vtk"""
        self.TairK = TairK
        nodes = self.get_nodes() 
        self.VarMu = VariableMu 
        #cells = [[bottom 1], [top 1], [bottom 2], [top 2]]
        #bottom 1 = x_node: link to segment bellow or parent seg bellow
        #top 1: y_node: link to segment above
        #bottom 2: x_node: link to child bellow or parent cell above
        #top 2: y_node: link to child above 
        #only bottom 1 and top 1 are used to define cellCenter
        
        orgs = self.plant.getOrgans()

        #cellFaceIDs
        bot1 = bot1bu = np.concatenate((list(map(lambda x: np.array(list(map(lambda x: np.array(x)[0], x.getSegments()))), orgs))))
        top1 = top1bu = np.concatenate((list(map(lambda x: np.array(list(map(lambda x: np.array(x)[1], x.getSegments()))), orgs))))
        unique, counts = np.unique(bot1, return_counts=True)
        unique=unique[np.logical_and(counts ==2 , unique!=0)]
        first_index = np.array([np.where(bot1 == ui)[0][0] for ui in unique])
        snd_index = [np.where(bot1 == ui)[0][1] for ui in unique]
        bot2 = np.full(len(bot1), -1)
        top2 = np.full(len(bot1), -1)
        new_indx = np.array(range(max(top1)+1, max(top1)+1 + len(unique)))
        if(len(orgs) > 2): #more than just main root and stem
            bot2[first_index] =bot2[snd_index] = new_indx
            bot1[snd_index] =top2[first_index -1]  = np.array(range(max(new_indx)+1, max(new_indx)+1+len(unique)))
        cells =  np.vstack((bot1,top1, bot2, top2))
        cells = MA.masked_values(cells, -1)
        print(cells)
        #node IDs maps and faceVertexIDs:
        oldNodesId = np.unique(self.get_segments().flatten())
        nodesBU = np.concatenate((bot1bu,bot1bu, top1bu,top1bu))
        newNodes =  np.concatenate((bot1,bot2, top1,top2))
        newNodesId = [list(np.unique(newNodes[np.where(np.logical_and(nodesBU == bi, newNodes !=-1))])) for bi in oldNodesId]
        self.old2newNodeID = dict(zip(oldNodesId,newNodesId ))      
        teilold = np.concatenate([np.full(len(self.old2newNodeID[bi]),bi) for bi in oldNodesId ])
        teilnew = np.concatenate([self.old2newNodeID[bi] for bi in oldNodesId ])
        self.new2oldNodeID = dict(zip(teilnew,teilold ))
                
        faces = np.array((np.arange(0,  max(np.maximum(top2, top1)) + 1),), dtype=np.int64) #index of verticies on each face
        print(faces)
        #vertices coords
        nodes_x_coord = [nodes[self.new2oldNodeID[ni]][0] for ni in faces[0]]     
        nodes_y_coord = [nodes[self.new2oldNodeID[ni]][1] for ni in faces[0]]     
        nodes_z_coord = [nodes[self.new2oldNodeID[ni]][2] for ni in faces[0]]     
        vertices = np.vstack((nodes_x_coord,nodes_y_coord,nodes_z_coord)) #change how coords are stored 
        
        #wat. psi. per cell
        if(not isinstance(self.Px, float)):#get mean rx of the segment (used by GrSink)
            self.rxCells = np.mean([[self.Px[self.new2oldNodeID[xi]] for xi in cells[0]], 
                                   [self.Px[self.new2oldNodeID[xi]] for xi in cells[1]]], axis=0)
            #rxCell >= rxThrMax -> rxFactor = 1
        else:
            self.rxCells = np.minimum(self.rxThrMax,self.Px)
        
        #cell length        
        length = np.array([norm(nodes[self.new2oldNodeID[cells[0][i]]]-nodes[self.new2oldNodeID[cells[1][i]]]) for i in range(len(cells[0]))])
        
        #radii per cell and face 
        radiiVertices = np.array([self.rs.radii[max((self.new2oldNodeID[xi]-1,0))] for xi in faces[0]]) #radius for each node, seed node take radius of 1st segment  of root  
        self.radiiCells = np.array([self.rs.radii[self.new2oldNodeID[xi]-1] for xi in cells[1]])
        
        #fipymesh
        #set radius of 0.0013cm to get same volume as with zhou2020
        #radiusPloemCells = (0.01**0.5)/3.14 #0.01=max size phloem cell https://github.com/granar/granar_examples/blob/master/modelparam/Wheat_F_Atkinson_2017.xml
        #6e-4 radius wheat phloem cell => Thompson and Holbrook (2003)
        self.mesh = Mesh1Dmod( radiiVertices = np.full(len(radiiVertices),0.0013), 
            radiiCells = np.full(len(self.radiiCells),0.0013),length = length, vertexCoords=vertices, 
            faceVertexIDs=faces, cellFaceIDs=cells) #1e-2 to go from radius of segment to radius of phloem
        
        #volume per cell    
        self.volumesSEG = np.array([(np.pi * self.radiiCells[xi]**2) *self.mesh.length[xi] for xi in range(self.mesh.numberOfCells) ])
        
        #map node2orgId
        self.orgID = np.concatenate([np.full(o.getNumberOfNodes()-1, o.getId()) for o in orgs])
        
        #store max growth rate and max organ length for computaiton of the CW-limited growth
        self.orgLength = np.array([sum(length[np.where(self.orgID == oid )[0]]) for oid in self.orgID])
        organTypes = self.get_organ_types()
        subTypes = self.get_subtypes()
        self.orgLengthThr = np.array( [self.rs.organParam[organTypes[self.new2oldNodeID[xi] - 1]][subTypes[self.new2oldNodeID[xi] -1] + 2*(organTypes[self.new2oldNodeID[xi] - 1]==4)].getParameter('lmax') for xi in cells[1]])
        self.maxGrowth = np.array( [self.rs.organParam[organTypes[self.new2oldNodeID[xi] - 1]][subTypes[self.new2oldNodeID[xi] -1]+ 2*(organTypes[self.new2oldNodeID[xi] - 1]==4)].getParameter('r')/(24*60*60) for xi in cells[1]]) #maxgrwth rate from cm/d to cm/s
         
        self.phi = CellVariablemod(name="solution variable", mesh=self.mesh,value = 0., hasOld = True)
        self.SucSink = CellVariablemod(name="total sucrose sink", mesh=self.mesh,value = 0., hasOld = True) #sucrose sink
        self.Mesophyll= CellVariablemod(name="leaf vacuole", mesh=self.mesh,value = 0., hasOld = True) #sucrose sink

    
    
    def _setOutflow(self): #only outflow when all segments have rm gr satisfaction?
        nodeYRoot = self.get_segments_index(2) + 1 #ynode_index of root
        nodeYRoot_newID = np.array([self.old2newNodeID[xi][0] for xi in nodeYRoot]).flatten() #new index
        self.rootCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in nodeYRoot_newID]) 
        
        #outFlowMax = self.rootCells*self.phi*self.radConductivity * ((self.mesh.radiiCells)**2)*self.mesh.length#* self.satisfaction #reset for security 
        #here rad of sefg or of phloem?
        #self.outFlow =(self.phi > 1e-10) *(self.phi - 1e-10)*self.rootCells*self.radConductivity * ((self.mesh.radiiCells)**2)*self.mesh.length#*self.phi* ((np.minimum(self.phi, self.phiThrMax)- self.phiThrMin)/(self.phiThrMax- self.phiThrMin))#*((self.phi -  np.min(self.phi, self.phiThrMin))/(0.003- self.phiThrMin))#(outFlowMax * (self.phi - outFlowMax > self.phiThrMin) + self.phi * (self.phi - outFlowMax <= self.phiThrMin))* (self.phi > self.phiThrMin)
        #self.outFlowMax = self.phi*self.rootCells*self.radConductivity * ((self.radiiCells)**2)*self.mesh.length#*self.phi* ((np.minimum(self.phi, self.phiThrMax)- self.phiThrMin)/(self.phiThrMax- self.phiThrMin))#*((self.phi -  np.min(self.phi, self.phiThrMin))/(0.003- self.phiThrMin))#(outFlowMax * (self.phi - outFlowMax > self.phiThrMin) + self.phi * (self.phi - outFlowMax <= self.phiThrMin))* (self.phi > self.phiThrMin)
        #self.outFlowMax = 2.82627e+95 * self.phi / (  1e+99 + self.phi)
        #self.outFlowMax =1e-3*0.12*0.2*sidesurface* self.phi /(60*60)#after arametrisation with piaf-munch. 60*600 hr -> s, *12-> sucrose to C (no need?)
        #==>
        self.outFlowMax =self.rootCells*(2.4e-04* (self.phi*12*1000) * np.pi*self.radiiCells*2*self.mesh.length)/(12*1000*3600)/self.mesh.cellVolumes
        print('out factor ',2.4e-05,  np.pi,self.radiiCells,2,self.mesh.length)
        #PiafMunchoutflowmax = (( 0.12*2*1e-3*np.pi* lengths[i] * surfacesOrg[i])* C_ST[i]);/

    def _setintCoeff(self): #use to get right value of k/mu* RT
        # R=8.314⋅m3⋅Pa⋅K−1mol−1 = 8.314⋅e6 mL⋅Pa⋅K−1mol−1
        #for fipy, use resistivity (conductivity) and not resistance (conductance). See: https://www.ctcms.nist.gov/fipy/documentation/numerical/discret.html#diffusion-term-nabla-cdot-left-gamma-1-nabla-phi-right
        waterViscosity = 2.414e-5 * (10**(247.8/(self.TairK - 140))) #in [Pa.s = N-s/m^2] https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
        sucroseDensity = 1.59 #g/cm³, https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose#section=Density
        sucroseMolarMass = 342.3 #g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose
        sucroseMolarVolume = sucroseMolarMass/sucroseDensity #cm³/mol like the one used by lacointe 2019
        if( self.VarMu):
            self.VolFractSucrose = self.phi.faceValue* sucroseMolarVolume #  #[mol/cm³] / [cm³/mol] = cm³/cm³ * sucroseMolarVolume
        else:
            self.VolFractSucrose = (0.5e-3)* sucroseMolarVolume
        #valid for VolFractSucrose <= 0.65 <-> [sucrose] <= 0.003019281 mol/cm3
        self.phloemViscosity = waterViscosity * np.exp((4.68 * 0.956 * self.VolFractSucrose)/(1 - 0.956 * self.VolFractSucrose)) #in [Pa.s = N-s/m^2]
        mu = self.phloemViscosity * 1e-6 /( 60*60 )#Pa s-1 to MPa hr-1
        k = self.phloemConductivity ## conductivity (cm2)
        #r_ST = mu/k * self.mesh.length / (self.radiiVertices**2)*np.pi #
        #self.phloemConductivity = (self.mesh.radiiVertices**2)/8#specific conductivity. assume that there are no sieve plate (see eq 5 from Thomphson 2006)
        R =  self.Param['R']*1e6  #Pa mL K-1 mol-1 
        #conductance = self.phloemConductivity/self.phloemViscosity  * np.pi*self.mesh.radiiVertices**2 /self.mesh._cellDistances ##ml-1 Pa s-1
        self.intCoeff = self.phloemConductivity/self.phloemViscosity * R * self.TairK # [cm/(Pa s)]* R[Pa mL K-1 mol-1] * T[K] = (cm mL)/(s mol-1)  => added in divergence: C (mol ml-1) * A (cm2) / L (cm) => cm3/s
        #self.dP_ST = self.phi.faceGrad * R * self.TairK #MPa mL
        #print('int coef ', 1/(conductance)*1e-6/3600,  r_ST , self.dP_ST)
        #print('S and L ,',np.pi*self.mesh.radiiVertices**2 ,self.mesh._cellDistances)
        #conductivity2conductance = #add surface/lenth to go from conductivity to conductance
        #logfilecondu = open('results/conductivity_10b.txt', "a")        
        #logfilecondu.write(','.join([num for num in map(str,self.phloemConductivity)])  +'\n')
        #logfilecondu.close()
        # flow = (k/mu) * RT * d(C)/dx = [cm2]/[Pa s] * [Pa * mL /(K mol)] * K * d[mol/mL] 
        # = [cm2 Pa-1 s-1] * d[Pa ] = cm2 s-1 (?)

    def _setAnSource(self, simDuration): #sourceterm in leaf segment
        leafCells = self.get_segments_index(4) + 1 #ynode_index of root
        leafCells_newID = np.array([self.old2newNodeID[xi][0] for xi in leafCells]).flatten() #new index
        leafCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in leafCells_newID]) 
        
        self.Agross = CellVariablemod(name = 'An', mesh = self.mesh, value = 0.)
        # find seg index with leaf, should be ordered like the An vector
        nodeYLeaf = self.get_segments_index(4) + 1 #ynode_index of leaves
        nodeYLeaf_newID = np.array([self.old2newNodeID[xi] for xi in nodeYLeaf]).flatten() #new index
        surfaceSide = (2*np.pi*((self.radiiCells)/100)*(self.mesh.length/100)) #m2. 
        Ag = np.add(self.An , self.Rd)
        #print('ag ', Ag, 'an ',self.An,'rd ', self.Rd) #in mol m-2 s-1
        for i, an in enumerate(Ag): #1.8e-11/self.mesh.cellVolumes in zhou2020
                
                if (an <= 0 or (simDuration >= 12 and simDuration < 24) ):
                    sucrose = 0 
                else:
                    v = self.mesh.cellVolumes[self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]]
                    s = surfaceSide[self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]] /2 #only assimilation on uper part of leaf?
                    #sucrose = 6.40314e-05*10/3600/1000#an*s/v / 12#mol sucrose s-1 ml-1,6.403140e-06#
                    #==>
                    sucrose = 6.40314e-004/3600/1000/v/12
                    #print('suc vol ',sucrose, v)
                self.Agross.constrain(sucrose, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #(sucrose* surfaceSide /self.mesh.cellVolumes, where= self.mesh.cellFaceIDs[1] == nodeYLeaf_newID[i]) #0 node = y node for stem
        self.JS_Mesophyll = np.maximum(1.33e-5*((self.Mesophyll-self.phi)*1000*12)*leafCells,self.phi*0)/self.mesh.cellVolumes/(1000*12*3600) #from Uys2020 but with lower k coef to avoid error, 1e-1*(self.Mesophyll-self.phi)/(1+self.Mesophyll+self.phi) 
        #print("js.meso", self.JS_Mesophyll)
        
    def _setRmSink(self): # Daudet 2002: Rm = (k1 + k2 * C)Sr. 
        self.RmMax = ((self.rhoC*1000) * self.volumesSEG *((self.phi*1000*12) * self.k1 + self.k2))/1000/12/3600/self.mesh.cellVolumes
        #RespMaint[i] =(5e-7* phi*1000*12+1e-7)*StructC[i]*1000/(60*60);//Rm
        
        #print('Rm, ',self.RmMax*1000*12*3600*self.mesh.cellVolumes)
       
    def _setGrSink(self): #Rg + C allocation = Grsink = growth * (1/GrowEff)
    
        #length and rx :  factor from 0 to 1 
        tipOrgans = np.concatenate((self.get_organ_nodes_tips()))
        tipOrgans_newID = np.array([self.old2newNodeID[xi] for xi in tipOrgans]).flatten() 
        tipOrgans = sum([self.mesh.cellFaceIDs[1] == xi for xi in tipOrgans_newID])
        
        leafCells = self.get_segments_index(4) + 1 #ynode_index of root
        leafCells_newID = np.array([self.old2newNodeID[xi][0] for xi in leafCells]).flatten() #new index
        leafCells = sum([self.mesh.cellFaceIDs[1] == xi for xi in leafCells_newID]) 
        
        growthCells = tipOrgans + leafCells > 0 #either tip of organ or leaf segment
        
        leafSegLength = (self.mesh.length/self.orgLength) * leafCells + 1 * (1-leafCells) #for leaf cells, weight of sink separated according to length of segment
        
        rxFactor =1#((self.rxCells - self.rxThrMin)/(self.rxThrMax - self.rxThrMin)) *(self.rxCells > self.rxThrMin)
        ####
        #Att!
        #
        ###
        
        self.lengthFactor =1# (self.orgLength*(-self.orgLength + self.orgLengthThr)/((self.orgLengthThr/2)**2)) * (self.orgLength < self.orgLengthThr)
        maxCNeeded = (((self.maxGrowth*np.pi * (self.radiiCells**2))*self.rhoC /self.GrEff))
        
        #maxCconcentationNeeded = (((self.maxGrowth*np.pi * (self.mesh.radiiCells*1e2)**2)*self.rhoC /self.GrEff))/self.mesh.cellVolumes
        #print("4grmax ", tipOrgans, leafCells,maxCconcentationNeeded, leafSegLength , self.maxGrowth, self.mesh.length,self.orgLength, rxFactor, self.lengthFactor)
        self.GrMax =  growthCells*maxCNeeded/12 *  rxFactor * self.lengthFactor * leafSegLength/self.mesh.cellVolumes# * self.CSat 
        #print('GrMax ,',self.GrMax*1000*12*3600*self.mesh.cellVolumes)
        
    def getCLimitedSinks(self, sink_dot, dt):
        #C-limited maintenance sink
        self.Rm = np.maximum(np.full(len(sink_dot),0.),np.minimum(sink_dot , self.RmMax.value))
        #C-limited growth sink
        self.GrSink =  np.maximum(np.full(len(sink_dot),0.),np.minimum(sink_dot - self.Rm, self.GrMax))
        
        #C-limited outflow
        self.outFlow =np.maximum(np.full(len(sink_dot),0.),np.minimum(sink_dot - self.Rm - self.GrSink, self.outFlowMax.value))
        #print('Gr ',self.totalSinkFactor,self.RmMax,  self.GrSink, self.orgID, np.array([max(0,xi) for xi in self.GrSink]))
        
        self.Growth0 = self.GrSink* self.mesh.cellVolumes *12/ self.rhoC * self.GrEff  #get increase in cm³ s-1
        self.Growth = self.Growth0 /(np.pi * (self.radiiCells**2))*dt#get length increase (cm)
        self.orgGr = np.array([sum(self.Growth[np.where(self.orgID == xi)[0]]) for xi in range(max(self.orgID )+1)]) #growth per organ
        Uerror = (self.Rm + self.GrSink + self.outFlow)/( self.RmMax + self.GrMax +self.outFlowMax)-self.Ufact
        print('lim climsink ',max(abs(Uerror)),
            self.Ufact[np.where([Uerror==max(abs(Uerror))][0])])
        
        self.sink_dot = sink_dot
        #print("self.sink_dot - self.Rm - self.GrSink - self.outFlow", self.sink_dot, self.sink_dot - self.Rm - self.GrSink - self.outFlow, self.outFlow)
        #print('minmax, ', self.RmMax.value*dt -         self.Rm,self.GrMax*dt- self.GrSink,  self.outFlowMax.value*dt-self.outFlow)
        #print(self.outFlowMax)
        #print('getCLimitedSinks, max error ', max(abs(self.sink_dot - self.Rm - self.GrSink - self.outFlow)))
        #print('ssdot ',self.sink_dot )
        if any(abs(( self.Rm + self.GrSink + self.outFlow)/( self.RmMax + self.GrMax +self.outFlowMax)-self.Ufact)> 0.01*self.Ufact):
            print((self.Rm + self.GrSink + self.outFlow)/( self.RmMax + self.GrMax +self.outFlowMax)-self.Ufact, self.Ufact,np.where([Uerror==max(abs(Uerror))][0]))
            raise Exception('getCLimitedSinks')

    @property        
    def organGr(self):#per cell not node
    #ATT: changed to dummy value for testing
        return np.array([0.01 for xi in range(max(self.orgID )+1)]) #self.orgGr

        
    def resetValues(self, simDuration):
        self._setAnSource(simDuration)
        self._setintCoeff()
        self._setRmSink()
        self._setGrSink()
        self._setOutflow()
        
        #self.totalSinkFactor = (self.RmMax + self.GrMax +self.outFlowMax)/( self.phi + (self.RmMax + self.GrMax + self.outFlowMax)*1.3)#0.4e-3*0.5 )
        
        
        #(self.phi+(self.RmMax + self.GrMax + self.outFlowMax)*3600)#)/self.phi )* ((self.RmMax + self.Gr )/self.phi <= 0.8) + 0.8 * ((self.RmMax + self.Gr )/self.phi >= 0.8))*(self.phi>0)
        #self.phiLim = (( self.phi - self.lim+ ((self.phi-self.lim)**2 )**0.5)/2)#max val between 0 and phi-lim
        #self.totalSinkFactor = (((self.RmMax + self.GrMax +self.outFlowMax)+self.phiLim)\
         #   -(((self.RmMax + self.GrMax +self.outFlowMax)+self.phiLim)**2 \
          #  - (4*0.9*(self.RmMax + self.GrMax +self.outFlowMax)*self.phiLim))**0.5)/(2*0.9) #min val between (self.RmMax + self.GrMax +self.outFlowMax) and self.phiLim
        #self.test = np.minimum(self.RmMax + self.GrMax +self.outFlowMax, self.phi -self.lim) #does update at each sweep
        #self.totalSinkFactor = np.maximum(np.minimum((self.RmMax + self.GrMax +self.outFlowMax)*12*1000*3600, (self.phi -self.lim)*12*1000), self.phi*0 )
        self.Ufact = np.maximum((np.minimum(self.phi, 0.3e-3/12) -self.lim)/(0.3e-3/12 - self.lim), self.phi*0 )
        self.totalSinkFactor = (self.RmMax + self.GrMax +self.outFlowMax)*self.Ufact
        #print((np.minimum(self.phi, 0.3e-3/12) -self.lim)/(0.3e-3/12 - self.lim))
        #print('totalSinkFactor ', min(self.totalSinkFactor.value), max(self.totalSinkFactor.value), sum(self.totalSinkFactor.value)/len(self.totalSinkFactor.value))
        #print('max(self.phi - self.RmMax ) ',max(self.phi - self.RmMax ),min(self.phi - self.totalSinkFactor ))
        #print("tot sink ", self.totalSinkFactor)
        
    def get_nodes_child_base(self):
        segments = self.get_segments() 	
        nodesx =[xi[0] for xi in segments] #take x node of all the segments
        unique, counts = np.unique(nodesx, return_counts=True)
        N3 = unique[np.logical_and(counts > 1 , unique != 0)] #take index of node which are at the base of 2 segment exept for seed node => all nodes linking 3 segments
        return(N3)
        
        
    def doCheck(self,cumulAn, cumulOut,cumulRm,cumulGr, initphi, timeSinceLastCheck , timeSpanCheck, numcheck,issue, issueRes,  issueLoop,phiConcentrationError, dt):
        phimin = self.phi * (self.phi < 0 )
        
        print('sucrose concentration < 0 at ',self.phi[np.where(self.phi < 0 )[0]],np.where(self.phi < 0 )[0])
        if(self.VarMu):
            print(' sucrose vol. fraction > 0.65 at ',
                np.where(self.VolFractSucrose > 0.65 )[0],self.VolFractSucrose[np.where(self.VolFractSucrose > 0.65 )[0]] )
        #print(self.sink_dot, self.Rm, self.GrSink, self.outFlow, self.mesh.cellVolumes)
        #print('Y nan? ',self.Rm, self.mesh.cellVolumes,self.sink_dot)
        if(any(self.sink_dot !=0)):
            print('current sink ratio of tot sink: Rm = ',int(np.around(sum(self.Rm* self.mesh.cellVolumes)/sum(self.sink_dot* self.mesh.cellVolumes)*100)),
                '%, Rg + Gr = ',int(np.around(sum(self.GrSink* self.mesh.cellVolumes)/sum(self.sink_dot* self.mesh.cellVolumes)*100)),'%, Out = ',
                int(np.around(sum(self.outFlow* self.mesh.cellVolumes)/sum(self.sink_dot* self.mesh.cellVolumes)*100)),'%' )
        print('cumulated loading: ', sum(cumulAn) + sum(initphi* self.mesh.cellVolumes), 'cumulated growth sink: ',sum(cumulGr), 'cumulated maintenance sink: ',sum(cumulRm), 
            'cumulated outflow: ',sum(cumulOut),'total suc. content in plant: ',sum(self.phi* self.mesh.cellVolumes), 
            'mean, max, min suc. concentration in plant (mol ml-1): ',sum(self.phi.value)/len(self.phi.value), ", ", max(self.phi.value),", ",min(self.phi.value))
        error1 =-sum(initphi* self.mesh.cellVolumes)-sum(cumulAn) + sum(self.phi* self.mesh.cellVolumes) \
            +sum(cumulOut+cumulRm+cumulGr)
        print('suc_plant + Rm + Rg + Gr + Out  - source :',error1)#+ sum(self.SucSink* self.mesh.cellVolumes))
        print('suc_plant +  TotSink -  source BIS :',-sum(initphi* self.mesh.cellVolumes)-sum(cumulAn) + sum(self.phi* self.mesh.cellVolumes) \
            + sum(self.SucSink))#+ sum(self.SucSink* self.mesh.cellVolumes))
        #print('totalSinkFactor ', min(self.totalSinkFactor.value), max(self.totalSinkFactor.value), 
            #sum(self.totalSinkFactor.value)/len(self.totalSinkFactor.value))
        error3=sum(self.SucSink) - sum(cumulGr)  - sum(cumulRm) - sum(cumulOut)
        print('C_TotSink - C_Gr - C_Rm - C_Out :', error3)
        if(abs(error1) < abs(error3)):
            raise Exception('issue with totsink values')
        diff = -sum(initphi* self.mesh.cellVolumes) -sum(cumulAn) + sum(self.phi* self.mesh.cellVolumes)  + sum(self.SucSink)
        if abs(diff/ sum(self.phi* self.mesh.cellVolumes)) > 0.01:
            print(sum(initphi* self.mesh.cellVolumes) ,sum(cumulAn) , sum(self.phi* self.mesh.cellVolumes)  , sum(self.SucSink* self.mesh.cellVolumes))
            raise Exception("diff/ sum(phl.phi* phl.mesh.cellVolumes) > 0.001")
        #print("check min phi ",np.where(self.totalSinkFactor > 0)[0],self.phi[np.where(self.totalSinkFactor > 0)[0]])
        valsink = np.where(self.sink_dot > 0)[0]
        #if len(valsink)>0 and  min(self.phi[valsink ] )<self.lim*0.999 :
         #   print(valsink)
          #  print(self.phi[valsink ] )
           # print(self.phi[valsink ] -self.lim )
            #raise Exception("(self.totalSinkFactor > 0) and (self.phi < 1e-5)")
        ###
        #add a check t see ratio of Rm, Rg + Gr, Out of the last time step.
        #
        ##
        print("check? ", timeSinceLastCheck , timeSpanCheck, numcheck)
        #print('no convergence at step(s) ',issue)#, ' res: ', issueRes, ' loops ', issueLoop)
        #print('phi value error at step(s) ',phiConcentrationError)

    def printParaview(self,cumulAn, cumulOut,cumulRm,cumulGr, fluxes ,simDuration):
        reorderedPhi = np.array([x for _,x in sorted(zip(self.cellsID,self.phi.value))])
        reorderedOutFlow = np.array([x for _,x in sorted(zip(self.cellsID,cumulOut))])
        reorderedRm = np.array([x for _,x in sorted(zip(self.cellsID,cumulRm))])
        reorderedGrSink = np.array([x for _,x in sorted(zip(self.cellsID,cumulGr))])
        reorderedVolume = np.array([x for _,x in sorted(zip(self.cellsID,self.mesh.cellVolumes))])
        ana = pb.SegmentAnalyser(self.rs)
        ana.addData("phi",  np.around(reorderedPhi*1e10, 5)) 
        #need to do around or too many digit for paraiew. and add 1e10 otherwise just yields 0
        ana.addData("outFlow",  np.around(reorderedOutFlow/reorderedVolume*1e10, 5))
        ana.addData("Rm",  np.around(reorderedRm/reorderedVolume*1e10, 5))
        ana.addData("GrSink",  np.around(reorderedGrSink/reorderedVolume*1e10, 5))
        ana.addData("rx", self.Px)
        ana.addData("fluxes", fluxes) 
        ana.write("results/%s_example10b.vtp" %(simDuration), ["radius", "surface", "phi", "outFlow", "Rm", "GrSink", "rx","fluxes"])
        
        #print('phi4ana ',reorderedPhi, np.around(reorderedPhi*1e10, 5))

    def doLog(self, duration,  organTypes, dt, issue, issueRes, issueLoop, alldts, simDuration):
        segZ = self.get_meanSegZ()
        subTypes = self.get_subtypes()
        
        logfileZ = open('results/xyz_10b.txt', "a")
        logfileOrgType = open('results/OrgType_10b.txt', "a")
        logfilesubType = open('results/subType_10b.txt', "a")
        logfileCSink = open('results/CSink_10b.txt', "a")
        logfileInput = open('results/Input_10b.txt', "a")
        logfileInput2 = open('results/Input2_10b.txt', "a")
        logfileMeso = open('results/meso_10b.txt', "a")
        logfilephi = open('results/phi_10b.txt', "a")
        logfilerg = open('results/rg_10b.txt', "a")
        logfilerm = open('results/rm_10b.txt', "a")
        logfileFlow = open('results/Flow_10b.txt', "a")
        logfileOut = open('results/outFlow_10b.txt', "a")
        logfileVol = open('results/vol_10b.txt', "a")
        logfilergMAX = open('results/rgMAX_10b.txt', "a")
        logfilermMAX = open('results/rmMAX_10b.txt', "a")
        logfileOutMAX = open('results/outFlowMAX_10b.txt', "a")
        logfileTimeR = open('results/TimeR_10b.txt', "a")
        logfileTime = open('results/Time_10b.txt', "a")
        logfileError = open('results/Error_10b.txt', "w")
        logfileucoef = open('results/uCoef_10b.txt', "a")
        logfileJW = open('results/JW_10b.txt', "a")
        
        #all in content and not concentration
        
        logfileZ.write(','.join([num for num in map(str,[segZ[self.new2oldNodeID[xi] - 1] for xi in self.mesh.cellFaceIDs[1]])])  +'\n')
        logfileOrgType.write(','.join([num for num in map(str,[organTypes[self.new2oldNodeID[xi] - 1] for xi in self.mesh.cellFaceIDs[1]])])  +'\n')
        logfileOrgType.write(','.join([num for num in map(str,[subTypes[self.new2oldNodeID[xi] - 1] for xi in self.mesh.cellFaceIDs[1]])])  +'\n')
        logfilephi.write(','.join([num for num in map(str,self.phi.value* self.mesh.cellVolumes)])  +'\n')
        logfileCSink.write(','.join([num for num in map(str,self.SucSink.value* self.mesh.cellVolumes)])  +'\n')#total C sink
        logfileInput.write(','.join([num for num in map(str, self.Agross.value*self.mesh.cellVolumes)])  +'\n')
        logfileInput2.write(','.join([num for num in map(str, self.JS_Mesophyll.value*self.mesh.cellVolumes)])  +'\n')
        logfileMeso.write(','.join([num for num in map(str, self.Mesophyll.value*self.mesh.cellVolumes)])  +'\n')
        logfileFlow.write(','.join([num for num in map(str, self.phi * (self.phi.faceGrad*self.intCoeff ).divergence*self.mesh.cellVolumes)])  +'\n')#write('\n;'+repr((phl.phi.faceValue * phl.phi.faceGrad*phl.intCoeff ).divergence)[7:-2])
        logfilerg.write(','.join([num for num in map(str, self.GrSink*self.mesh.cellVolumes)])  +'\n')
        logfilerm.write(','.join([num for num in map(str, self.Rm*self.mesh.cellVolumes)])  +'\n')
        logfileOut.write(','.join([num for num in map(str, self.outFlow*self.mesh.cellVolumes)])  +'\n')
        logfilergMAX.write(','.join([num for num in map(str, self.GrMax*self.mesh.cellVolumes)])  +'\n')
        logfilermMAX.write(','.join([num for num in map(str, self.RmMax.value*self.mesh.cellVolumes)])  +'\n')
        logfileOutMAX.write(','.join([num for num in map(str, self.outFlowMax*self.mesh.cellVolumes)])  +'\n')
        logfileVol.write(','.join([num for num in map(str,self.mesh.cellVolumes)])  +'\n')
        logfileTimeR.write(repr(duration)+'\n')
        logfileTime.write(repr(simDuration)+'\n')#;print('repr(simDuration) ', repr(simDuration))
        #logfiledt.write(','.join([num for num in map(str,alldts)]) )#write('\n;'+repr(phl.outFlow.value)[7:-2])
        logfileError.write('no convergence at step(s):\n'+repr(issue)+ '\nres:\n'+repr( issueRes)+ '\nloops:\n'+repr( issueLoop))
        logfileucoef.write(','.join([num for num in map(str,self.Ufact)])  +'\n')
        logfileJW.write(','.join([num for num in map(str, (self.phi.faceGrad*self.intCoeff ).divergence)])  +'\n')
        
        logfileJW.close()
        logfileMeso.close()
        logfileCSink.close()
        logfileInput2.close()
        logfileInput.close()
        logfilephi.close()
        logfilerg.close()
        logfilerm.close()
        logfileOut.close()
        logfileFlow.close()
        logfileVol.close()
        logfilergMAX.close()
        logfilermMAX.close()
        logfileOutMAX.close()    
        logfileTimeR.close()    
        logfileTime.close()
        logfileZ.close()
        logfileOrgType.close()
        logfilesubType.close()
        logfileError.close()
        logfileucoef.close()
        #logfiledt.close()
        
        
        if(self.VarMu):
            logfileintcoef = open('results/intcoef_10b.txt', "a")
            logfileintcoef.write(','.join([num for num in map(str,self.intCoeff)])  +'\n')
            logfileintcoef.close()
            logfilevisco = open('results/visco_10b.txt', "a")
            logfilevisco.write(','.join([num for num in map(str,self.phloemViscosity)])  +'\n')
            logfilevisco.close()
        
    def valueError(self):
        if(any(self.phi > 0.003)): #risk of fraction phi going over 0.65
            opt = np.get_printoptions()
            np.set_printoptions(threshold=np.inf)
            print(self.mesh.cellFaceIDs[0])
            print(self.mesh.cellFaceIDs[1])
            print(self.mesh.cellFaceIDs[2])
            print(self.mesh.cellFaceIDs[3])
            np.set_printoptions(**opt)
            raise Exception("phl.phi > 0.003")
        if(any(self.phi < 0)):
            raise Exception("phl.phi < 0")
        if(any([math.isnan(ci) for ci in self.phi])):
            raise Exception("phl.phi is nan")
        if(any(self.Rm < 0)):
            raise Exception("phl.Rm < 0")
        if(any(self.GrSink < 0)):
            raise Exception("phl.GrSink < 0")
        if(any(self.outFlow < 0)):
            raise Exception("phl.outFlow < 0") 
        
def GetTime(seconds):
    sec = timedelta(seconds=int(seconds))
    d = datetime(1,1,1) + sec
    return(f'{d.day-1}d {d.hour}hr {d.minute}mn {d.second}s')        
        
