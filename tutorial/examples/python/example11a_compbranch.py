from fipy import *
from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D as Grid1D
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D
from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython
import plantbox as pb
import vtk_plot as vp
import math
import os
from io import StringIO


class comp2mesh(PhloemFluxPython):
    def __init__(self):
        """ @param rs is either a pb.MappedPlant or a string containing a rsml filename                     
        """
        self.mesh = [] #to store grid
        self.phi = [] # to store solution
        self.outFlow = 1.
        self.satisfaction = 0
        self.Rm = 0.
        self.GrSink = 0. # C used for Rg and growth
        self.Gr = 0.05 #C-limited growth during time step
        self.Growth = 0.05 #C-limited growth during time step
        self.RmMax = 0.
        self.orgGr = np.array([])
        self.phiFactor2= -1
        self.intCoeff = 100.
        self.lengthMesh= 0
        self.Source = 0
        self.outFlow = 0
        
    def resetValmesh2D(self,phl):
        self.lengthMesh = self.mesh.cellVolumes/(np.pi * phl.mesh.radiiCells**2)
    
        waterViscosity = 2.414e-5 * (10**(247.8/(phl.TairK - 140))) #in [Pa.s = N-s/m^2] https://www.engineersedge.com/physics/water__density_viscosity_specific_weight_13146.htm
        sucroseDensity = 1.59 #g/cm³, https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose#section=Density
        sucroseMolarMass = 342.3 #g/mol https://pubchem.ncbi.nlm.nih.gov/compound/Sucrose
        sucroseMolarVolume = sucroseMolarMass/sucroseDensity #cm³/mol
        VolFractSucrose = self.phi.faceValue* sucroseMolarVolume#np.minimum(np.full(len(self.phi.faceValue),0.65),np.array(self.phi.faceValue* sucroseMolarVolume) ) #[mol/cm³] / [cm³/mol] = cm³/cm³ * sucroseMolarVolume
        phloemViscosity = waterViscosity * np.exp((4.68 * 0.956 * VolFractSucrose)/(1 - 0.956 * VolFractSucrose))
        R = phl.Param['R'] * 1e6 # Pa * mL /(K mol)
        self.intCoeff = 10.#phl.phloemConductivity/phloemViscosity * R * phl.TairK #(self.osmoCoeff)
                
        self.RmMax = phl.rhoC * self.mesh.cellVolumes * (self.phi * phl.k1 + phl.k2 )
        self.Rm = (self.phi > phl.phiThrMin)*( self.RmMax * (self.phi - self.RmMax >phl.phiThrMin )+ (self.phi )* (self.phi - self.RmMax <= phl.phiThrMin ) )#*(self.phi >0)#*(self.phi> self.phiThrGrowth))
        self.CSat = (self.Rm == self.RmMax) 

        maxCconcentationNeeded = (((phl.maxGrowth*np.pi * phl.mesh.radiiCells**2)*phl.rhoC /phl.GrEff))\
            * (self.lengthMesh/phl.orgLength)/self.mesh.cellVolumes
        self.Gr =  maxCconcentationNeeded/12 * 1 * phl.lengthFactor * self.CSat #* abs(self.phiFactor1)
        self.phiFactor2 = ((self.phi - self.Gr) > phl.phiThrMin)  
        self.GrSink =  (self.phi > phl.phiThrMin)*((self.Gr  * self.phiFactor2 \
        + (self.phi - phl.phiThrMin) * (~self.phiFactor2)))
        self.Growth0 = self.GrSink / phl.rhoC * phl.GrEff * self.mesh.cellVolumes #get increase in cm³
        self.Growth = self.Growth0 /(np.pi * phl.mesh.radiiCells**2)#get length increase (cm)
        self.satisfaction = self.CSat * (self.GrSink == self.Gr)
        
        
        outFlowMax = phl.rootCells*self.phi* phl.radConductivity \
        * (phl.mesh.radiiCells**2)*self.lengthMesh* self.satisfaction
        self.outFlow = outFlowMax * (self.phi - outFlowMax > 0) * (self.phi > 0) \
        + self.phi * (self.phi - outFlowMax <= 0)* (self.phi > 0)
        

######################
#
# plant
#
####################### 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "oneroot_mgiraud" #"manyleaves"#"smallPlant"# "oneroot" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")
start = 1.#5 #1.275

""" Parameters xylem"""
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004 #  cm3/day radial conductivity between xylem and guard cell
p_s = -200  # static water potential (saturation) 33kPa in cm
#p_g = -2000 # water potential of the guard cell
RH = 0.5 # relative humidity
TairC = 20
p_a =  -1000  #default outer water potential 
#simtime = start = 15 # [day] for task b
k_soil = []
Q = 900e-6 # mol quanta m-2 s-1 light, example from leuning1995
cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
TairK = TairC + 273.15


es = 0.61078 * math.exp(17.27 * TairC / (TairC + 237.3)) 
ea = es * RH 
VPD = es - ea 


""" Parameters phloem """
""" soil """
min_ = np.array([-50, -50, -150])
max_ = np.array([90, 40, 10])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

pl.initialize()
pl.simulate(start, False)
phl = PhloemFluxPython(pl)
segs = pl.getPolylines() #only first segment
phl.mp2mesh(segs) #creates grid

cellsIDOld = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]
comp = comp2mesh()
comp.mesh = phl.meshes
comp.phi = CellVariable(mesh = comp.mesh, value= 0., hasOld = True)

dt = 0.9*( 0.04 ** 2) / (2 *10)
phl.Source = CellVariable(mesh= phl.mesh, value = 0.)
comp.Source = CellVariable(mesh= phl.meshes, value = 0.)

phl.Source.constrain(100., where = phl.mesh.cellFaceIDs[1] == 5)
comp.Source.constrain(100., where = phl.meshes.cellFaceIDs[0] == 0)


cumulatedIn = 0.
cumulatedInPlant = 0.
steps =1
comp.phi.updateOld()
phl.phi.updateOld()
print(phl.mesh.cellFaceIDs,phl.meshes.cellFaceIDs)
phl.resetValues()
comp.resetValmesh2D(phl)
lim = 2
issueLoop=issueRes= []
cumulOut = 0.
cumulRm = 0.
cumulGr = 0.
for _ in range(steps):
          
    eq = (TransientTerm(var = comp.phi) == (comp.phi.faceValue * comp.phi.faceGrad).divergence+  comp.Source + comp.Rm + comp.GrSink + comp.outFlow) #
    #change for term rather than variable
    #see if, when fill mesh (no mask), affects diffusion
    eqplant = (TransientTerm(var = phl.phi) == (phl.phi.faceValue * phl.phi.faceGrad).divergence+ phl.Source + phl.Rm + phl.GrSink + phl.outFlow)
    print('steps ', _, dt)
    res2 =res=1e5
    loop =0 
    
    while ((res2 > 1e-10 or res > 1e-10) and loop < 100):# 1e-6* max(abs(phl.phi.value)) or res2 > 1e-6* max(abs(phl.Source.value))) and loop < 1000:#(res > 1e-6* max(abs(Source.value)) or res > 1e-6* max(abs(phi.value)) or res2 > 1e-6* max(abs(phl.phi.value)) or res2 > 1e-6* max(abs(phl.Source.value))) and loop < 1000:
            res =  eq.sweep(dt= dt)
            res2 =  eqplant.sweep( dt= dt)
            loop += 1
            print('res: ',res, res2)
            print('Rm ', comp.Rm,'\n', phl.Rm)
            print('Grsink ',comp.GrSink, '\n', phl.GrSink)
            print ('outflow ',comp.outFlow, '\n', phl.outFlow)
            print(phl.satisfaction ,phl.CSat , (phl.GrSink == phl.Gr))
            res = max(res,res2)
    if(loop > 100):
        print('no convergence! res: ', res, res2)
        issueRes = np.append(issueRes, [res2])
        issueLoop = np.append(issueLoop, [loop])
    comp.phi.updateOld()
    phl.phi.updateOld()
    cumulatedInPlant += sum(phl.Source * dt*phl.mesh.cellVolumes)
    cumulOut += sum(phl.outFlow * dt*phl.mesh.cellVolumes)
    cumulRm += sum(phl.Rm * dt*phl.mesh.cellVolumes)
    cumulGr +=  sum(phl.Gr * dt*phl.mesh.cellVolumes)
    print('mass balance ', cumulGr +cumulRm +cumulOut+sum(phl.phi*phl.mesh.cellVolumes)- cumulatedInPlant) #sum(phi*phl.meshes.cellVolumes) - cumulatedIn, 
    print(cumulGr ,cumulRm ,cumulOut,sum(phl.phi*phl.mesh.cellVolumes), cumulatedInPlant)
    print('phi 2',phl.phi,'\n',comp.phi)
    dt = 0.9*( min(phl.mesh.length) ** 2) / (2 * max((phl.phi.faceValue * phl.phi.faceGrad).divergence))
   
print('n conv at ', issueLoop, issueRes)
'''
print(phl.mesh.faceCellIDs[0],'\n',phl.mesh.faceCellIDs[1])
#print(phl.meshes.faceCellIDs[0][phl.meshes.interiorFaces],'\n',phl.meshes.faceCellIDs[1][phl.meshes.interiorFaces])
#print(phl.meshes.faceCellIDs)
print(phl.meshes._cellToFaceOrientations)
print(phl.meshes._orientedFaceNormals)
print(phl.mesh._cellToFaceOrientations)
print(phl.mesh._orientedFaceNormals)
#faceCellIDs ok?
#_cellToFaceOrientations?
#
print('faceValue\n',phl.phi.faceValue, '\nface value end\n')
print('faceGrad\n',phl.phi.faceGrad, '\nface grad end\n')
print('phi 2',phl.phi,'\n', phi, '\n', sum(phi*phl.meshes.cellVolumes),sum(phl.phi*phl.mesh.cellVolumes))
print(phl.meshes.cellCenters)
print(phl.mesh.cellCenters)
print('phi ', phl.phi, '\n',phi)
print('source ', phl.Source, '\n',Source)
print("facegradintdiff\n",(phl.phi.faceValue * phl.phi.faceGrad).divergence * dt, '\n',(phi.faceValue*phi.faceGrad).divergence*dt)
phi.setValue(0.)
phl.phi.setValue(0.)
phl.Source = CellVariable(mesh= phl.mesh, value = 0.)
Source = CellVariable(mesh= phl.meshes, value = 0.)
print('phi empty', phl.phi, '\n',phi)
print('source empty', phl.Source, '\n',Source)
phl.Source.constrain(100., where = phl.mesh.cellFaceIDs[0] == 16)
Source.constrain(100., where = phl.meshes.cellFaceIDs[0] == 37)
phi.setValue(phi.value + Source.value)
phl.phi.setValue(phl.phi.value + phl.Source.value)
print('phi ', phl.phi, '\n',phi)
print('source ', phl.Source, '\n',Source)
print("facegradintdiff\n",(phl.phi.faceValue * phl.phi.faceGrad).divergence * dt, '\n',(phi.faceValue*phi.faceGrad).divergence*dt)

#prblm: 
#lat cell connected to upper of parent cell but not bottom of other parent cell
#
#
#
#





'''

"""
print('\n\n\n\n')
steps =100
for _ in range(steps):
    dt = 0.9*( 0.04 ** 2) / (2 *10)
            
    eq = (TransientTerm(var = phi) == DiffusionTerm(var = phi,coeff= 10 *phi.faceValue )+  Source )
    eqplant = (TransientTerm(var = phl.phi) == DiffusionTerm(var = phl.phi,coeff=10*phl.phi.faceValue )+ phl.Source)
    #print('results plant ',phl.phi + ((phl.phi.faceValue)* (phl.phi.faceGrad )).divergence*dt + phl.Source * dt)
    #print('results plant ',phi + ((phi.faceValue )* ( phi.faceGrad )).divergence*dt +  Source * dt)
    print('steps ', _, dt)
    #res =res2 = 1e5
    #loop =0 
    eq.solve(dt= dt)
    eqplant.solve( dt= dt)
    '''
    while (res > 1e-3* max(abs(Source.value)) or res > 1e-3* max(abs(phi.value)) or res2 > 1e-3* max(abs(phl.phi.value)) or res2 > 1e-3* max(abs(phl.Source.value))) and loop < 1000:
            res =  eq.sweep(dt= dt)
            res2 =  eqplant.sweep( dt= dt)
            loop += 1
            print('res: ',res, res2)
            #res = max(res,res2)
            print('maxRes: ',1e-3* max(abs(phi.faceValue)), 1e-3* max(abs(phl.phi.faceValue)) )
    '''
    cumulatedIn += sum(Source * dt*phl.meshes.cellVolumes)
    cumulatedInPlant += sum(phl.Source * dt*phl.mesh.cellVolumes)
    #phi.updateOld
    #phl.phi.updateOld
    print('mass balance ', sum(phi*phl.meshes.cellVolumes) - cumulatedIn, sum(phl.phi*phl.mesh.cellVolumes)- cumulatedInPlant)
    print('phi 2',phl.phi,'\n', phi, '\n', sum(phi*phl.meshes.cellVolumes),sum(phl.phi*phl.mesh.cellVolumes))
"""
'''
viewer =  MatplotlibViewer(vars = phi, xmin = -1, xmax = 2., ymin = -1.5, ymax = 1)#fipy.viewers.Viewer(vars=(phi,))
viewer.plot("results/output1.png")

print('volume\n', phl.mesh.cellVolumes, '\n',phl.meshes.cellVolumes)

print("faceValueinterior\n",phl.phi.faceValue,'\n', phl.phi.faceValue[np.where(phl.mesh.interiorFaces)[0]], '\n',phi.faceValue[np.where(phl.meshes.interiorFaces)[0]])
print("faceGradinterior\n",phl.phi.faceGrad[2][np.where(phl.mesh.interiorFaces)[0]], '\n',phi.faceGrad[0][np.where(phl.meshes.interiorFaces)[0]], '\n',
    phi.faceGrad[1][np.where(phl.meshes.interiorFaces)[0]])
print("facegradint\n",phl.phi.faceGrad.divergence, '\n',phi.faceGrad.divergence)#[np.where(phl.mesh.exteriorFaces)[0]]
print("facegradintdiff\n",(phl.phi.faceValue * phl.phi.faceGrad).divergence * dt, '\n',(phi.faceValue*phi.faceGrad).divergence*dt)#[np.where(phl.mesh.exteriorFaces)[0]]

print("facegradint\n",phl.phi.faceValue.divergence, '\n',phi.faceValue.divergence)#[np.where(phl.mesh.exteriorFaces)[0]]
print(phl.Source, Source)
print(phl.meshes.cellCenters)
'''
"""

######################

#tip indexes 
#tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips()
#tiproots_newID = np.array([phl.old2newNodeID[xi] for xi in tiproots]).flatten()
#tiproots=  [np.any([tiproots_newID == yi]) for yi in phl.mesh.faceVertexIDs[0]]
#tiprootsfaces= sum([phl.mesh.faceVertexIDs[0] == xi for xi in tiproots_newID]) 
#tiprootsfacesID= np.where(tiprootsfaces)[0]
#tiprootscells= sum([phl.mesh.cellFaceIDs[1] == xi for xi in tiprootsfaces])
#tiprootscellsID= np.where(tiprootscells)[0]

####
#
# Equations
#
####


cumulOutOld = CellVariable(mesh = phl.mesh, value=0.)
cumulGrOld = CellVariable(mesh = phl.mesh, value=0.)
cumulRmOld = CellVariable(mesh = phl.mesh, value=0.)
cumulAnOld = CellVariable(mesh = phl.mesh, value=0.)
phiOld = CellVariable(mesh = phl.mesh, value=0., hasOld = True)
cellVolumeOld = phl.mesh.cellVolumes

logfilermM = open('results/rmMax_10a.txt', "w")
logfilerm = open('results/rm_10a.txt', "w")
logfilephi = open('results/phi_10a.txt', "w")
logfilerg = open('results/rg_10a.txt', "w")
logfilergS = open('results/rgSink_10a.txt', "w")


phl.phi.updateOld()
growthSteps = []
issue = []
issueRes = []
issueLoop = []
phiConcentrationError = []

#phl.phi = CellVariable(mesh = phl.mesh, value= phiOld.value * cellVolumeOld/phl.mesh.cellVolumes, hasOld =True)
cumulOut = cumulOutOld
cumulGr =cumulGrOld
cumulRm = cumulRmOld
cumulAn = cumulAnOld

numsegRoot = len(phl.get_segments_index(2))
numsegStem = len(phl.get_segments_index(3))
#print(numsegRoot,numsegStem, len(phl.mesh.length[:numsegRoot]))
#print(np.flip(phl.mesh.length[:numsegRoot]),phl.mesh.length[numsegRoot:] )

#meshroot = Grid2D(nx = numsegRoot,ny = 1, dx = np.flip(phl.mesh.length[:numsegRoot]), 
 #   dy = np.mean(phl.mesh.cellVolumes[:numsegRoot]/phl.mesh.length[:numsegRoot])) + [[sum(phl.mesh.length[numsegRoot:])], [0]] 
#meshstem = Grid2D(nx = numsegStem,ny = 1, dx = phl.mesh.length[numsegRoot:], 
#    dy = np.mean(phl.mesh.cellVolumes[numsegRoot:]/phl.mesh.length[numsegRoot:]))  #np.flip(
comp = comp2mesh()
comp.lengthMesh = np.concatenate((np.flip(phl.mesh.length[:numsegRoot]), phl.mesh.length[numsegRoot:]) )
comp.mesh = Grid2D(nx = len(comp.lengthMesh ),ny = 1, dx = comp.lengthMesh , 
    dy = np.mean(comp.lengthMesh))#meshroot + meshstem

print('numsegRoot ', numsegRoot, comp.lengthMesh)
comp.phi = CellVariable(mesh = comp.mesh, value=0., hasOld = True)
comp.resetValmesh2D(phl)

mesh2meshID = np.concatenate((np.flip(np.array(range(numsegRoot))), np.array(range(numsegRoot,len( phl.mesh.length)))))
######
#
#step and init konz
#
#######
#print(list(range(1/len(phl.mesh.cellFaceIDs[1]))))
initval = np.array(range(0,len(phl.mesh.cellFaceIDs[1])))*1e-2
initval2 = np.concatenate((np.flip(initval[:numsegRoot]),initval[numsegRoot:]) )
phl.phi.setValue(initval2)
comp.phi.setValue(initval)
#phl.phi.setValue(2e-4, where = phl.mesh.cellFaceIDs[1] == (phl.nodeYLeaf_newID -1))
#comp.phi.setValue(2e-4, where = phl.mesh.cellFaceIDs[1] == (phl.nodeYLeaf_newID -1))
phl.phi.updateOld()
comp.phi.updateOld()
steps = 0

phl.resetValues()    
comp.resetValmesh2D(phl)
#sourceInit = phl.Source + 1e-4
#print(sourceInit )
#dt =  0.9 * min(phl.mesh.length) ** 2 / (2 *phl.intCoeff) *2
#dtSave = dt
eq = (TransientTerm(var = phl.phi) == DiffusionTerm(var = phl.phi,coeff= phl.intCoeff * phl.phi.faceValue)+ phl.Source - phl.Rm - phl.GrSink  - phl.outFlow) #+ sourceInit

eqcomp = (TransientTerm(var = comp.phi) == DiffusionTerm(var =  comp.phi,coeff=  comp.intCoeff *  comp.phi.faceValue)+  phl.Source -  comp.Rm   -  comp.GrSink -  comp.outFlow) #+ sourceInit


print('volume\n', np.take(phl.mesh.cellVolumes,mesh2meshID ), '\n',comp.mesh.cellVolumes)

print("faceValueinterior\n",np.take(phl.phi.faceValue[np.where(phl.mesh.interiorFaces)[0]],mesh2meshID[:-1] ), '\n',comp.phi.faceValue[np.where(comp.mesh.interiorFaces)[0]])
print("facegradint\n",np.take(phl.phi.faceGrad.divergence,mesh2meshID ), '\n',comp.phi.faceGrad.divergence)#[np.where(phl.mesh.exteriorFaces)[0]]


print("phi + source\n",np.take(phl.phi + phl.Source ,mesh2meshID ) , '\n',
    comp.phi + phl.Source )

print("diff\n",np.take(  (phl.phi.faceValue*phl.intCoeff*phl.phi.faceGrad).divergence,mesh2meshID ), '\n',
    (comp.phi.faceValue* comp.intCoeff* comp.phi.faceGrad).divergence)
    
print("faceval\n",np.take(  (phl.phi.faceValue).divergence,mesh2meshID ), '\n',
    (comp.phi.faceValue).divergence)
print("facegrad\n",np.take(phl.phi.faceGrad[2][np.where(phl.mesh.interiorFaces)[0]],mesh2meshID[:-1])
    , '\n',    comp.phi.faceGrad[0][np.where(comp.mesh.interiorFaces)[0]])
print("facegrad\n",np.take(  (phl.phi.faceGrad).divergence,mesh2meshID ), '\n',
    ( comp.phi.faceGrad).divergence)

print('mesh1 ', phl.mesh._cellDistances[np.where(phl.mesh.interiorFaces)[0]] , '\n', comp.mesh._cellDistances[np.where(phl.mesh.interiorFaces)[0]], 
    '\n', 'mesh2\n', phl.mesh._faceAreas[np.where(phl.mesh.interiorFaces)[0]],  '\n',comp.mesh._faceAreas[np.where(phl.mesh.interiorFaces)[0]])

#print('mesh1a ', phl.mesh._cellVolumes  -comp.mesh._cellVolumes, 
 #   '\n', 'mesh2a\n',phl.mesh._cellAreas , '\n', comp.mesh._cellAreas[1], "\n",comp.mesh._cellAreas[3], )##mesh1a ok. 2D: area = length or width
print('mesh3 ', np.take(phl.mesh._cellDistances ,mesh2meshID ), '\n', comp.mesh._cellDistances, '\n', 'mesh4\n',np.take( phl.mesh._faceAreas,mesh2meshID ),  '\n',comp.mesh._faceAreas)        

for _ in range(steps):
    if(_ == 10):
        phl.Source.setValue(0.)
    temp = abs(np.take(phl.Gr,mesh2meshID)-comp.Gr)
    temp2 = abs(np.take(phl.Gr,mesh2meshID)-comp.Gr)/abs(np.take(phl.Gr,mesh2meshID))*100
    #print('\n\nmax diff ', max(temp)," Xmol/cm3 ", max(temp2), "%")
    if( max(temp2) > 1e-3):
        print("difference! gr")
        print(np.where(temp2 == max(temp2)))
        break
        
        
    temp = abs(np.take(phl.GrSink,mesh2meshID)-comp.GrSink)
    temp2 = abs(np.take(phl.GrSink,mesh2meshID)-comp.GrSink)/abs(np.take(phl.GrSink,mesh2meshID))*100
    print('\n\nmax diff ', max(temp)," Xmol/cm3 ", max(temp2), "%")
    if( max(temp2) > 1e-3):
        print("difference! grsink")
        print(np.where(temp2 == max(temp2)))
        break
        
    dt =  0.01 * min(phl.mesh.length) ** 2 / (2 * max(phl.intCoeff* phl.phi.faceValue) )
    print("\n\nstep ", _, start + _ * dt, dt)
    '''
    '''
    #sys.stdout = NullIO()  
    ####
    #
    # plant growth
    #
    ####
    pl.simulate(dt , False)
    phl = PhloemFluxPython(pl)
    segs = pl.getPolylines() #first segment
    
    
    ####
    #
    # xylem
    #
    ####
    
    nodes = phl.get_nodes()
    tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips() #end node of end segment of each organ
    node_tips = np.concatenate((tiproots, tipstem, tipleaf))
    tiproots, tipstem, tipleaf = phl.get_organ_segments_tips() #end segment of each organ
    seg_tips = np.concatenate((tiproots, tipstem, tipleaf))


    phl.setKr([[kr],[kr_stem],[gmax]]) 
    phl.setKx([[kz]])
    phl.airPressure = p_a

    phl.seg_ind = seg_tips # segment indices for Neumann b.c.
    phl.node_ind = node_tips
    phl.Px = phl.solve_leuning( sim_time = dt,sxx=[p_s], cells = True, Qlight = Q,VPD = VPD,
        Tl = TairK,p_linit = p_s,ci_init = cs,cs=cs, soil_k = [], log = False)
    fluxes = phl.radial_fluxes(dt, phl.Px, [p_s], k_soil, True)  # cm3/day
    
    ####
    #
    # phloem
    #
    ####
    
    phl.mp2mesh(segs) #creates grid
    
    if(len(phl.phi.value) > len(phiOld)): #place according to index ID
        print('len ',len(phl.phi.value) , len(phiOld))
        CphiOld = phiOld * cellVolumeOld
        orderedCPhi = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, CphiOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
        #print(orderedCPhi, phl.cellsID,phl.mesh.cellVolumes)
        orderedCPhiNew = np.take(orderedCPhi, phl.cellsID)/phl.mesh.cellVolumes#np.concatenate((CphiOld , np.full(len(phl.phi.value)-len(phiOld),0.)))/phl.mesh.cellVolumes#
        growthSteps = np.append(growthSteps, [_])
        phl.phi = CellVariable(mesh = phl.mesh, value= orderedCPhiNew , hasOld =True)
        
        phl.phi.updateOld()
        
        orderedcumulOut = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulOutOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
        orderedcumulOutNew = np.take(orderedcumulOut, phl.cellsID)        #np.concatenate((cumulOutOld , np.full(len(phl.phi.value)-len(phiOld),0.)))#
        cumulOut = CellVariable(mesh = phl.mesh, value= orderedcumulOutNew)
        
        orderedcumulGr = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulGrOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
        orderedcumulGrNew = np.take(orderedcumulGr, phl.cellsID)   #np.concatenate((cumulGrOld , np.full(len(phl.phi.value)-len(phiOld),0.)))#             
        cumulGr = CellVariable(mesh = phl.mesh, value= orderedcumulGrNew)
        
        orderedcumulRm = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulRmOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
        orderedcumulRmNew = np.take(orderedcumulRm, phl.cellsID)       #np.concatenate((cumulRmOld , np.full(len(phl.phi.value)-len(phiOld),0.)))#         
        cumulRm = CellVariable(mesh = phl.mesh, value= orderedcumulRmNew)
        
        orderedcumulAn = np.concatenate((np.array([x for _,x in sorted(zip(cellsIDOld, cumulAnOld))]), np.full(len(phl.phi.value)-len(phiOld),0.)))        
        orderedcumulAnNew = np.take(orderedcumulAn, phl.cellsID)   #np.concatenate((cumulAnOld , np.full(len(phl.phi.value)-len(phiOld),0.)))#             
        cumulAn = CellVariable(mesh = phl.mesh, value= orderedcumulAnNew)
    
    else: 
        phl.phi = CellVariable(mesh = phl.mesh, value= phiOld.value * cellVolumeOld/phl.mesh.cellVolumes, hasOld =True)
        cumulOut = cumulOutOld
        cumulGr =cumulGrOld
        cumulRm = cumulRmOld
        cumulAn = cumulAnOld
    phl.resetValues()    
    comp.resetValmesh2D(phl)
    '''
    res =res2 =  1e+10
    resOld = 2e+10
    loop = 0
    
    #sys.stdout = sys.__stdout__
    #print('rm1 ',comp.Rm   ,  comp.GrSink ,   comp.outFlow,comp.intCoeff )print('check if value updates ',sum(phl.CSat), sum(phl.RmMax), sum(phl.Rm), sum(phl.Gr), sum(phl.GrSink), sum(phl.outFlow)) #check if value updates
        
    #print('check if value updates ',phl.intCoeff * phl.phi.faceValue) #check if value updates
    while res > max(1e-5, 1e-5* max(phl.phi)) and loop < 1000:# and resOld != res:
        resOld = res
        res =  eq.sweep(dt= dt)
        res2 =  eqcomp.sweep(dt= dt)
        loop += 1
        print('res :',res,res2)
        res = max(res,res2)
    if res > max(1e-5, 1e-5 * max(phl.phi)) :
        print('no convergence! res: ', res)
        issue = np.append(issue, [_])
        issueRes = np.append(issueRes, [res])
        issueLoop = np.append(issueLoop, [loop])
        
        
    phl.phi.updateOld()
    comp.phi.updateOld()
    #print('rm2 ',comp.Rm   ,  comp.GrSink ,  comp.outFlow,comp.intCoeff )
    cumulRm.setValue(cumulRm.value + phl.Rm.value * dt* phl.mesh.cellVolumes)
    cumulAn.setValue(cumulAn.value + phl.Source.value * dt* phl.mesh.cellVolumes)
    cumulGr.setValue(cumulGr.value + phl.GrSink.value * dt* phl.mesh.cellVolumes)
    cumulOut.setValue(cumulOut.value + phl.outFlow.value * dt* phl.mesh.cellVolumes)
    
    if _ % 1 == 0: #choose how often get output       
        ####
        #
        # copy of phi
        #
        ####
        '''
        #print('sorted ',phl.cellsID,phl.phi.value)
        reorderedPhi = np.array([x for _,x in sorted(zip(phl.cellsID,phl.phi.value))])
        reorderedOutFlow = np.array([x for _,x in sorted(zip(phl.cellsID,cumulOut))])
        reorderedRm = np.array([x for _,x in sorted(zip(phl.cellsID,cumulRm))])
        reorderedGrSink = np.array([x for _,x in sorted(zip(phl.cellsID,cumulGr))])
        orgNum = np.array([phl.newCell2organID[xi] for xi in phl.newCell2organID])
        reorderedorgNum = np.array([x for _,x in sorted(zip(phl.cellsID,orgNum))])
        ana = pb.SegmentAnalyser(phl.rs)
        ana.addData("phi",  np.around(reorderedPhi, 5))
        ana.addData("outFlow",  np.around(reorderedOutFlow, 5))
        ana.addData("Rm",  np.around(reorderedRm, 5))
        ana.addData("GrSink",  np.around(reorderedGrSink, 5))
        #ana.addData("rx", phl.Px)
        #ana.addData("fluxes", fluxes) 
        ana.addData("orgNr", reorderedorgNum) 
        ana.write("results/%s_example10a.vtp" %(_), ["radius", "surface", "phi", "outFlow", "Rm", "GrSink","orgNr"])# "rx","fluxes",
        '''
        
        ####
        #
        # check balance
        #
        ####
        #print('coef ', phl.intCoeff * phl.phi.faceValue)
        phimin = phl.phi * (phl.phi < 0 )
        print('sink priority check ',sum(phl.GrSink * (~phl.CSat)),sum(phl.outFlow*~phl.satisfaction)) #tip of stem negative at the beginning
        print('sucrose concentration < 0 at ',phl.phi[np.where(phl.phi < 0 )[0]],np.where(phl.phi < 0 )[0], ' sucrose vol. fraction > 0.65 at ',np.where(phl.VolFractSucrose > 0.65 )[0])
        #,sum(sourceInit * dt * phl.mesh.cellVolumes)
        print('cumulated loading: ', sum(initval2* phl.mesh.cellVolumes)+sum(cumulAn) 
            ,'cumulated growth sink: ',sum(cumulGr), 'cumulated maintenance sink: ',sum(cumulRm), 'cumulated outflow: ',sum(cumulOut),'total C content in plant: ',sum(phl.phi* phl.mesh.cellVolumes))
        print('C_plant - C_Sink + C_source :',-sum(initval2* phl.mesh.cellVolumes) -sum(cumulAn) + sum(phl.phi* phl.mesh.cellVolumes) + sum(cumulGr)+ sum(cumulRm)+ sum(cumulOut))
        print('check if value updates ',sum(phl.CSat), sum(phl.RmMax), sum(phl.Rm), sum(phl.Gr), sum(phl.GrSink), sum(phl.outFlow)) #check if value updates
        #print( 'Growth sink in segments with C limitation for maintenance (should be 0): ',) #check if value updates
        
        logfilephi.write('\n'+repr(phl.phi.value)[7:-2])
        logfilerm.write('\n'+repr(phl.Rm.value)[7:-2])
        logfilermM.write('\n'+repr(phl.RmMax.value)[7:-2])
        logfilergS.write('\n'+repr(phl.GrSink.value)[7:-2])
        logfilerg.write('\n'+repr(phl.Gr.value)[7:-2])
        print('no convergence at step(s) ',issue, ' res: ', issueRes, ' loops ', issueLoop)
        print('phi value error at step(s) ',phiConcentrationError)
        #print('became bigger at step(s) ',growthSteps)
    
    sourceInit = phl.Source * 0.
    #eq = TransientTerm(var = phl.phi) == DiffusionTerm(var = phl.phi,coeff= phl.intCoeff * phl.phi.faceValue)\
     #   - phl.Rm   - phl.GrSink + phl.Source - phl.outFlow 

    #eqcomp = TransientTerm(var = comp.phi) == DiffusionTerm(var =  comp.phi,coeff=  comp.intCoeff *  comp.phi.faceValue)\
     #   -  comp.Rm   -  comp.GrSink +  phl.Source -  comp.outFlow 
    if(sum(abs(phl.GrSink * (~phl.CSat)))+sum(abs(phl.outFlow*~phl.satisfaction))\
        + sum(abs(phl.phi[np.where(phl.phi < 0 )[0]])) != 0):
        print("error priority")
        break
        
        
        
    
    #print("faceValueinterior\n",phl.phi.faceValue[np.where(phl.mesh.interiorFaces)[0]], '\n',comp.phi.faceValue[np.where(comp.mesh.interiorFaces)[0]])
    #print("faceValueextrior\n",np.take(phl.phi.faceGrad.divergence,mesh2meshID),
    #'\n',comp.phi.faceGrad.divergence)
    #print("phi\n",phl.phi,'\n', np.take(phl.phi,mesh2meshID), '\n',comp.phi)#,"\n", comp.phi)
    
    
    
    temp = abs(np.take(phl.phi,mesh2meshID)-comp.phi)
    temp2 = abs(np.take(phl.phi,mesh2meshID)-comp.phi)/abs(np.take(phl.phi,mesh2meshID))*100
    print('\n\nmax diff phi', max(temp)," Xmol/cm3 ", max(temp2), "%")
    if( max(temp2) > 1e-1):
        print("difference!")
        print(temp2)
        print(np.where(temp2 == max(temp2)))
        break
        
        
        
    temp = abs(np.take(phl.Gr,mesh2meshID)-comp.Gr)
    temp2 = abs(np.take(phl.Gr,mesh2meshID)-comp.Gr)/abs(np.take(phl.Gr,mesh2meshID))*100
    print('\n\nmax diff gr', max(temp),"\nXmol/cm3\n", temp, '\n', temp2)
    if( max(temp2) > 1e-3):
        print(np.take(phl.Gr,mesh2meshID),comp.Gr)
        print("difference! gr")
        print(np.where(temp2 == max(temp2)))
        break
        
        
    temp = abs(np.take(phl.GrSink,mesh2meshID)-comp.GrSink)
    temp2 = abs(np.take(phl.GrSink,mesh2meshID)-comp.GrSink)/abs(np.take(phl.GrSink,mesh2meshID))*100
    print('\n\nmax diff grrsink', max(temp),"\nXmol/cm3\n", temp, '\n', temp2)
    if( max(temp2) > 1e-3):
        print(np.take(phl.GrSink,mesh2meshID),comp.GrSink)
        print("difference! grsink")
        print(np.where(temp2 == max(temp2)))
        break
        
        
    temp = abs(np.take(phl.RmMax,mesh2meshID)-comp.RmMax)
    temp2 = abs(np.take(phl.RmMax,mesh2meshID)-comp.RmMax)/abs(np.take(phl.RmMax,mesh2meshID))*100
    print('\n\nmax diff RmMax', max(temp),"\nXmol/cm3\n", temp, '\n', temp2)
    if( max(temp2) > 1e-3):
        print(np.take(phl.RmMax,mesh2meshID),comp.RmMax)
        print("difference! RmMax")
        print(np.where(temp2 == max(temp2)))
        break
        
        
    temp = abs(np.take(phl.Rm,mesh2meshID)-comp.Rm)
    temp2 = abs(np.take(phl.Rm,mesh2meshID)-comp.Rm)/abs(np.take(phl.Rm,mesh2meshID))*100
    print('\n\nmax diff Rm', max(temp),"\nXmol/cm3\n", temp, '\n', temp2)
    if( max(temp2) > 1e-3):
        print(np.take(phl.Rm,mesh2meshID),comp.Rm)
        print("difference! Rm")
        print(np.where(temp2 == max(temp2)))
        break
        
    temp = abs(np.take(phl.outFlow,mesh2meshID)-comp.outFlow)
    temp2 = abs(np.take(phl.outFlow,mesh2meshID)-comp.outFlow)/abs(np.take(phl.outFlow,mesh2meshID))*100
    print('\n\nmax diff out', max(temp),"\nXmol/cm3\n", temp, '\n', temp2)

    phiOld = phl.phi.copy()
    cumulOutOld = cumulOut.copy()
    cumulGrOld = cumulGr.copy()
    cumulRmOld = cumulRm.copy()
    cumulAnOld = cumulAn.copy()
    cellVolumeOld = phl.mesh.cellVolumes
    cellsIDOld = phl.cellsID
    organGr = phl.organGr * dt
    #print('organGr ',organGr , dt)
    #pl.setCWLimGr(organGr , dt)
print('rm\n',np.take(phl.phi,mesh2meshID),'\n',np.take(phl.RmMax,mesh2meshID),'\n',np.take(phl.Rm,mesh2meshID))
print('rootcells\n',np.take(phl.rootCells,mesh2meshID),'\n',np.take(phl.CSat,mesh2meshID),'\n', np.take(phl.satisfaction,mesh2meshID))
logfilermM.close()
logfilerm.close()
logfilephi.close()
logfilerg.close()
logfilergS.close()
"""