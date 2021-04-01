from fipy import *

from fipy.tools import numerix as np
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from phloem_flux import PhloemFluxPython
import plantbox as pb
import vtk_plot as vp
import math


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

R = T= medPhi = kz = mu= 1. #TODO redefine
phl.osmoCoeff = R * T* medPhi* (kz/mu) 
phl.setD() #update D


steps = 10000
timeStepDuration = 0.1
eqX =TransientTerm() ==  DiffusionTerm(coeff= phl.D) + phl.Source - ImplicitSourceTerm( phl.exteriorCoeff.divergence ) #+((phl.mesh.faces&leafface ) * fluxBottom).divergence - ImplicitSourceTerm(exteriorCoeff.divergence)
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
        #print(_)
        #print(phl.phi)
        vw = VTKCellViewer(vars=(phl.phi))
        vw.plot(filename="results/%s_example8a.vtk" %(_))


