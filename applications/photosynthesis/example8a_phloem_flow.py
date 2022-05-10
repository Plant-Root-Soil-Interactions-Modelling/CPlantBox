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


dir_name = "/home/rbtlm2004/DUMUX/CPlantBox/tutorial/examples/python/results"
test = os.listdir(dir_name)
for item in test:
    if item.endswith("8a.vtk"):
        os.remove(os.path.join(dir_name, item))
        
        
np.set_printoptions(threshold=sys.maxsize)

######################
#
# plant
#
####################### 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves"# "oneroot" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

pl.initialize()
pl.simulate(14, False)
phl = PhloemFluxPython(pl)
segs = pl.getPolylines() #segments regrouped per organ
phl.An = np.linspace(60000.,1000000., len(phl.get_segments_index(4)))
phl.mp2mesh(segs) #creates grid


######################
#mesh = Grid2D(nx = len(phl.mesh.length),ny = 1, dx = phl.mesh.length, dy = np.mean(phl.mesh.cellVolumes/phl.mesh.length))

#tip indexes 
tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips()
tiproots_newID = np.array([phl.old2newNodeID[xi] for xi in tiproots]).flatten()
tiproots=  [np.any([tiproots_newID == yi]) for yi in phl.mesh.faceVertexIDs[0]]
tiprootsfaces= sum([phl.mesh.faceVertexIDs[0] == xi for xi in tiproots_newID]) 
tiprootsfacesID= np.where(tiprootsfaces)[0]
tiprootscells= sum([phl.mesh.cellFaceIDs[1] == xi for xi in tiprootsfaces])
tiprootscellsID= np.where(tiprootscells)[0]

####
#
# Equations
#
####


eq1 = TransientTerm(var = phl.phi) == DiffusionTerm(var = phl.phi,coeff= phl.intCoeff * phl.phi.faceValue) - phl.Rm   - phl.GrSink + phl.Source#root
eq2 = TransientTerm(var = phl.phi) == ImplicitSourceTerm(var = phl.phi,coeff = (phl.extCoeff * phl.phi.faceGrad).divergence) 

cumulOut = CellVariable(mesh = phl.mesh, value=0.)
cumulGr = CellVariable(mesh = phl.mesh, value=0.)
cumulRm = CellVariable(mesh = phl.mesh, value=0.)
cumulAn = CellVariable(mesh = phl.mesh, value=0.)

logfilermM = open('8a_rmMax.txt', "w")
logfilerm = open('8a_rm.txt', "w")
logfilephi = open('8a_phi.txt', "w")
logfilerg = open('8a_rg.txt', "w")
logfilergS = open('8a_rgSink.txt', "w")


phl.phi.updateOld()
dt = min(0.1, 0.9 * min(phl.mesh.length) ** 2 / (2 * phl.osmoCoeff))
steps = 1000
issue = []
issueRes = []
issueLoop = []
for _ in range(steps):
    res = 1e+10
    resOld = 2e+10
    loop = 0
    print('step ', _)
    while res > max(1e-10, 1e-5 * max(phl.phi)) and loop < 1000 and resOld != res:
        resOld = res
        res =  eq1.sweep(dt= dt)
        loop += 1
    if res > max(1e-10, 1e-5 * max(phl.phi)) :
        print('no convergence! res: ', res)
        issue = np.append(issue, [_])
        issueRes = np.append(issueRes, [res])
        issueLoop = np.append(issueLoop, [loop])
        
    phl.phi.updateOld()
    cumulRm.setValue(cumulRm.value + phl.Rm.value * dt* phl.mesh.cellVolumes)
    cumulAn.setValue(cumulAn.value + phl.Source.value * dt* phl.mesh.cellVolumes)
    cumulGr.setValue(cumulGr.value + phl.GrSink.value * dt* phl.mesh.cellVolumes)
    
    loop= 0    
    resOld = 2e+10
    res = 1e+10
    while res > max(1e-10, 1e-5 * max(phl.phi))  and loop < 1000 and resOld != res:
        resOld = res
        res =  eq2.sweep(dt= dt)
        loop += 1
    if res > max(1e-10, 1e-5 * max(phl.phi)):
        print('no convergence! res: ', res)
        issue = np.append(issue, [_])
        issueRes = np.append(issueRes, [res])
        issueLoop = np.append(issueLoop, [loop])
    
    phl.phi.updateOld()
    cumulOut.setValue(cumulOut.value -  phl.phi * (phl.extCoeff * phl.phi.faceGrad).divergence * phl.mesh.cellVolumes * dt)
    
    vw = VTKCellViewer(vars=(phl.phi))
    vw.plot(filename="results/%s_example8a.vtk" %(_))
    
    phimin = phl.phi * (phl.phi < 0 )
    print(phl.phi[np.where(phl.phi < 0 )[0]],np.where(phl.phi < 0 )[0] ) #tip of stem negative at the beginning
    
    print('cumulated loading: ', sum(cumulAn), 'cumulated growth sink: ',sum(cumulGr), 'cumulated maintenance sink: ',sum(cumulRm), 'cumulated outflow: ',sum(cumulOut),'total C content in plant: ',sum(phl.phi* phl.mesh.cellVolumes))
    print('C_plant - C_Sink + C_source :', -sum(cumulAn) + sum(phl.phi* phl.mesh.cellVolumes) + sum(cumulGr)+ sum(cumulRm)+ sum(cumulOut))
    #print(sum(phl.CSat), sum(phl.RmMax), sum(phl.Rm), sum(phl.Gr), sum(phl.GrSink)) #check if value updates
    print( 'Growth sink in segments with C limitation for maintenance (should be 0): ',sum(phl.GrSink * (~phl.CSat))) #check if value updates
    
    logfilephi.write('\n'+repr(phl.phi.value)[7:-2])
    logfilerm.write('\n'+repr(phl.Rm.value)[7:-2])
    logfilermM.write('\n'+repr(phl.RmMax.value)[7:-2])
    logfilergS.write('\n'+repr(phl.GrSink.value)[7:-2])
    logfilerg.write('\n'+repr(phl.Gr.value)[7:-2])

print('no convergence at step(s) ',issue, ' res: ', issueRes, ' loops ', issueLoop)

logfilermM.close()
logfilerm.close()
logfilephi.close()
logfilerg.close()
logfilergS.close()
