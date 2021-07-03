from fipy import *
from fipy.meshes.nonUniformGrid1D import NonUniformGrid1D as Grid1D
from fipy.meshes.nonUniformGrid2D import NonUniformGrid2D as Grid2D
from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython, GetTime

from CellVariablemod import CellVariablemod
import plantbox as pb
import vtk_plot as vp
import math
import os
from io import StringIO
from datetime import datetime, timedelta

home_dir = os.getcwd()
dir_name = "/results"
dir_name2 = home_dir + dir_name
test = os.listdir(dir_name2)
for item in test:
    if item.endswith("10s.txt"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10s.vtk"):
        os.remove(os.path.join(dir_name2, item))
    if item.endswith("10s.vtp"):
        os.remove(os.path.join(dir_name2, item))
        
        
np.set_printoptions(threshold=sys.maxsize)

class NullIO(StringIO):
    def write(self, txt):
       pass
       

    
######################
#
# plant
#
####################### 
pl = pb.MappedPlant() 

path = "../../../modelparameter/plant/"  
name =  "oneroot_mgiraud"
pl.readParameters(path + name + ".xml")
start =2


""" soil """
min_ = np.array([-50, -50, -150])
max_ = np.array([90, 40, 10])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments
pl.initialize(True)

pl.simulate(start, False)



phl = PhloemFluxPython(pl)
phl.mp2mesh( VariableMu = False) #cst viscosity
organTypes = phl.get_organ_types()
ana = pb.SegmentAnalyser(phl.rs)
ana.write("results/example10s.vtp" )

cellsIDOld = [phl.new2oldNodeID[xi] - 1 for xi in phl.mesh.cellFaceIDs[1] ]


tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips()
print(tiproots, tipstem, tipleaf)
tiproots_newID = np.array([phl.old2newNodeID[xi] for xi in tiproots]).flatten()
tiproots=  [np.any([tiproots_newID == yi]) for yi in phl.mesh.faceVertexIDs[0]]
tiprootsfaces= sum([phl.mesh.faceVertexIDs[0] == xi for xi in tiproots_newID]) 
tiprootsfacesID= np.where(tiprootsfaces)[0]
tiprootscells = sum([phl.mesh.cellFaceIDs[1] == xi for xi in tiprootsfaces])
tiprootscellsID = np.where(tiprootscells)[0]


tipstem_newID = np.array([phl.old2newNodeID[xi] for xi in tipstem]).flatten()
tipstemfaces= sum([phl.mesh.faceVertexIDs[0] == xi for xi in tipstem_newID]) 
tipstemcells = sum([phl.mesh.cellFaceIDs[1] == xi for xi in tipstemfaces])

Input = 0.1 * tipstemcells #mol Suc ml-1 s-1 
Output = 0.9 * phl.phi * tiprootscells # mol Suc ml-1 s-1 
print('In and out ', Input, Output)      
simDuration=0
beginning = datetime.now()
step= 0
dtVal = dt= 1
increase = 1.5
growthSteps = []
issue = []
issueRes = []
issueLoop = []
phiConcentrationError = []

timeSpanCheck = 60*10
timeSinceLastCheck = np.inf
numcheck = 0

while simDuration < 36*60*60: #
    end = datetime.now()
    duration = end - beginning
    print('\n\nstep n°', step,', dt: ',dt,'s, ','tot sim time: ', GetTime(simDuration), ', computation time: ',duration)

    
    convergence = False
    
    if(step > 0 ): #convergence criteria according to phi value before loop (as might get odd values during loop)
        phi4res = phl.phi.copy() 
    else: #for first loop, phi == 0, so set phi4res as CellVariable => can update at each loop (otherwise convergence never reached)
        phi4res = phl.phi 
    while (not convergence):
        eq = (TransientTerm(var = phl.phi) ==   (phl.phi.faceValue * phl.phi.faceGrad*phl.intCoeff ).divergence + Input - Output )
        res = 1e+10
        resOld = 2e+10
        loop = 0
        
        while ( ( res > 1e-5*max(phi4res* phl.mesh.cellVolumes)or np.isnan(res) or any(phl.phi <0)) and loop <= 100 ) : 
            resOld = res
            res = eq.sweep(dt= dt) #quid solver?
            loop += 1
            
        if loop > 100:
            print('no convergence! res: ', res,' lim = ', 1e-5*max(phi4res* phl.mesh.cellVolumes))
            issue = np.append(issue, [step])
            issueRes = np.append(issueRes, [res])
            issueLoop = np.append(issueLoop, [loop])
            dtVal = dtVal/2
            dt = dtVal
            maxdt = 1
            increase = 1.1
            print("change dt from ", dtVal*2, " to ", dtVal)
        else:
            print('convergence reached')
            convergence = True
           
    phl.phi.updateOld()
    print(phl.phi)
    print(phl.intCoeff )
    simDuration += dt
    timeSinceLastCheck +=dt
    
    if timeSinceLastCheck >= timeSpanCheck: #choose how often get output 
        segZ = phl.get_meanSegZ()
        numcheck +=1
        end = datetime.now()
        duration = end - beginning
        
        logfileZ = open('results/xyz_10s.txt', "a")
        logfileInput = open('results/Input_10s.txt', "a")
        logfilephi = open('results/phi_10s.txt', "a")
        logfileFlow = open('results/Flow_10s.txt', "a")
        logfileOut = open('results/outFlow_10s.txt', "a")
        logfileVol = open('results/vol_10s.txt', "a")
        logfileTime = open('results/Time_10s.txt', "a")
        logfileJW = open('results/JW_10s.txt', "a")
        
        logfileZ.write(','.join([num for num in map(str,[segZ[phl.new2oldNodeID[xi] - 1] for xi in phl.mesh.cellFaceIDs[1]])])  +'\n')
        logfileInput.write(','.join([num for num in map(str, Input*phl.mesh.cellVolumes)])  +'\n')
        logfilephi.write(','.join([num for num in map(str,phl.phi.value* phl.mesh.cellVolumes)])  +'\n')
        logfileFlow.write(','.join([num for num in map(str, phl.phi * (phl.phi.faceGrad*phl.intCoeff ).divergence*phl.mesh.cellVolumes)])  +'\n')
        logfileOut.write(','.join([num for num in map(str, Output*phl.mesh.cellVolumes)])  +'\n')
        logfileVol.write(','.join([num for num in map(str,phl.mesh.cellVolumes)])  +'\n')
        logfileTime.write(repr(simDuration)+'\n')
        logfileJW.write(','.join([num for num in map(str, (phl.phi.faceGrad*phl.intCoeff ).divergence)])  +'\n')
        
        logfileZ.close()
        logfileInput.close()
        logfilephi.close()
        logfileFlow.close()
        logfileOut.close()
        logfileVol.close()
        logfileTime.close()
        logfileJW.close()
        
        
        timeSinceLastCheck = 0.
        
    
    step += 1
    
    if loop < 10:
        dtVal = min(dtVal*increase, maxdt)
    dt = dtVal
    
    
    
end = datetime.now()
duration = end - beginning
print('duration simulation ', duration.seconds, 's for a simulation of ', step ,' steps and ', GetTime(simDuration),  "\n", dt)

