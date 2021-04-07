from fipy import *

from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython
import plantbox as pb
import vtk_plot as vp
import math

dir_name = "/home/rbtlm2004/DUMUX/CPlantBox/tutorial/examples/python/results"
test = os.listdir(dir_name)
for item in test:
    if item.endswith("8a.vtk"):
        os.remove(os.path.join(dir_name, item))

np.set_printoptions(threshold=sys.maxsize)
# plant 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
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
phl.An = np.linspace(600.,1000., len(phl.get_segments_index(4)))
print(phl.An)
phl.mp2mesh(segs) #creates grid

R = T=  kz = mu= 1. #TODO redefine
phl.osmoCoeff = R * T* (kz/mu) 
phl.setD() #update D
tiproots, tipstem, tipleaf = phl.get_organ_nodes_tips()
tiproots_newID = np.array([phl.old2newNodeID[xi] for xi in tiproots]).flatten()
tiproots=  [np.any([tiproots_newID == yi]) for yi in phl.mesh.faceVertexIDs[0]]
tiprootscells= np.where(sum([phl.mesh.faceVertexIDs[0] == xi for xi in tiproots_newID]) ==1)[0]

logfilerm = open('8a_rm.txt', "w")
logfilephi = open('8a_phi.txt', "w")
logfilerg = open('8a_rg.txt', "w")

steps = 10000
timeStepDuration = 0.1
eqX =TransientTerm() ==  DiffusionTerm(coeff= phl.D* phl.phi.faceValue)   + phl.Source   - phl.GrSink - phl.Rm 
#negative fr implicit source term?
for _ in range(steps):
    print('step: ', _)
    eqX.solve(var=phl.phi, dt=timeStepDuration) 
    phl.phi.updateOld()
    res = 1e+10
    loop = 0
    while res > 1e-4 and loop < 1000:
        res =  eqX.sweep(var=phl.phi,dt= timeStepDuration)
        #print(res)
        loop += 1
    if res > 1e-4:
        print('no convergence! res: ', res)
        break
    if __name__ == '__main__' and _%1000 ==0:
        vw = VTKCellViewer(vars=(phl.phi))
        vw.plot(filename="results/%s_example8a.vtk" %(_))
        logfilerm.write('\n'+repr(phl.Rm.value)[7:-2])
        logfilephi.write('\n'+repr(phl.phi.value)[7:-2])
        logfilerg.write('\n'+repr(phl.GrSink.value)[7:-2])


logfilerm.close()
logfilerm.close()
logfilerg.close()

