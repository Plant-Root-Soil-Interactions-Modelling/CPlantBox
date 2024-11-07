import os
onCloud = False # the notebook is not being run on a local machine but on a cloud-based service (e.g. mybinder, googlecolab)

if onCloud:
    import xvfbwrapper
    vdisplay = xvfbwrapper.Xvfb()
    vdisplay.start()
    os.environ["CUDA_VISIBLE_DEVICES"] = "1"
    doInteractiveImage = False
else:
    doInteractiveImage = True

    
sourcedir = os.getcwd()+"/../../.."
filedir = os.getcwd()
import sys; sys.path.append(sourcedir); sys.path.append(sourcedir+"/src")
import importlib
import plantbox as pb
import visualisation.vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt # for 2d plots
import numpy as np

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.photosynthesis_cpp import PhotosynthesisPython

""" root system """
simtime = 14  # [day] 
rs = pb.MappedRootSystem() # handles conductivity and mapping to the soil cells
path = "../../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
random_parameters = rs.getOrganRandomParameter(pb.root)
for p in random_parameters[1:]:
    p.dx = 0.25 
verbose = False

rs.initialize(verbose) # note that an artificial root with type = 0 is added 
rs.simulate(simtime,verbose)
_ = vp.plot_roots(pb.SegmentAnalyser(rs.mappedSegments()), "subType",interactiveImage = doInteractiveImage) 

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity [1/day]
p_top = -300  # top soil pressure [cm]
p0 = -500  # dirichlet bc at root collar [cm]
trans = -1.2  # neuman bc at root collar [cm3/day]

""" prepare soil matric potentials per segment"""
segs = rs.segments # MappedRootSystem has access to segments and nodes 
nodes = rs.nodes
p_s = np.zeros((len(segs),)) # soil total potentials around each root segment
for i, s in enumerate(segs):
    p_s[i] = p_top - 0.5 * (nodes[s.x].z + nodes[s.y].z)  # constant total potential (hydraulic equilibrium)

""" root problem """
r = XylemFluxPython(rs)
r.setKr([0., kr, kr , kr, kr, kr]) # no radial flux into the artificial root segment
r.setKx([1., kz, kz, kz, kz, kz])

""" Numerical solution """
#
# In the following 'cells = False' means total soil potentials 'p_s' is given for each segment
#

# rx = r.solve_neumann(simtime, trans, p_s, cells = False) # use Neumann bc 
# rx = r.solve_dirichlet(simtime, p0, 0, p_s, cells = False) # use Dirichlet bc
rx = r.solve(simtime, trans, 0, p_s, cells= False, wilting_point = -15000) # use Neumann, switch to Dirichlet if below wilting_point

fluxes1 = r.segFluxes(simtime, rx, p_s, cells = False)  # [cm3/day]
print("Transpiration", r.collar_flux(simtime, rx, [p_s], k_soil = [], cells = False), "cm3/day")