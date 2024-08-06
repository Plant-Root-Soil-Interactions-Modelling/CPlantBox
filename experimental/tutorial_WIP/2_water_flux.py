
import os
sourcedir = os.getcwd()+"/../.."

import sys; sys.path.append(sourcedir); sys.path.append(sourcedir+"/src")
import importlib
import plantbox as pb
import visualisation.vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt # for 2d plots
import numpy as np
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.photosynthesis_cpp import PhotosynthesisPython

isColab = False
if True:
    """ root system """
    simtime = 14  # [day] 
    rs = pb.MappedRootSystem() # handles conductivity and mapping to the soil cells
    path = "../../modelparameter/structural/rootsystem/"
    name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
    rs.readParameters(path + name + ".xml")
    random_parameters = rs.getRootRandomParameter()
    verbose = False
    for p in random_parameters[1:]:
        p.dx = 0.25 
    rs.initialize(verbose) # note that an artificial root with type =0 is added in MappedRootSystem
    rs.simulate(simtime,verbose)
    _ = vp.plot_roots(pb.SegmentAnalyser(rs.mappedSegments()), "subType",interactiveImage = not isColab) 

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


    surfs = np.multiply(np.array(r.rs.segLength()), 2*np.array(r.rs.radii)*np.pi)  # root segment side surface [cm2]
    fluxes = np.divide(fluxes1, surfs)  # we convert to [cm3/(cm2 day)]

    """ plot results """
    ana = pb.SegmentAnalyser(r.rs.mappedSegments())
    ana.addData("rx", rx) # add simulation result
    _ = vp.plot_roots(ana, "rx", "Xylem matric potential (cm)",interactiveImage = not isColab)  
    ana.addData("fluxes", fluxes) # add simulation result 
    _ = vp.plot_roots(ana, "fluxes", "Segment flux (cm/day)",interactiveImage = not isColab) 


    fig, ax = plt.subplots()
    name = ["root", "stem", "leaf"]
    color = ['tab:blue', 'tab:orange', 'tab:green']
    for ndType in [2, 3, 4]:
        #nodes = self.get_nodes()
        y = r.get_nodes_organ_type(ndType)#coordinates
        x = np.array(rx)[r.get_nodes_index(ndType)]
        ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
                   alpha=0.3, edgecolors='none')

    ax.legend()
    ax.grid(True)
    plt.xlabel("Xylem pressure (cm)")
    plt.ylabel("Depth (cm)")
    plt.title("Xylem matric potential (cm)")
    plt.show()


    fig, ax = plt.subplots()
    name = ["root", "stem", "leaf"]
    color = ['tab:blue', 'tab:orange', 'tab:green']
    for ndType in [2, 3, 4]:
        segIdx = r.get_segments_index(ndType)
        nodesy = segIdx + np.ones(segIdx.shape, dtype = np.int64)
        nodes = r.get_nodes()
        y = np.array(nodes)[nodesy]#coordinates
        x = fluxes[segIdx]
        ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
                   alpha=0.3, edgecolors='none')

    ax.legend()
    ax.grid(True)
    plt.xlabel("Fluxes (cm/day)")
    plt.ylabel("Depth (cm)")
    plt.title("water fluxes")
    plt.show()
        

simtime = 14  # [day] 
pl = pb.MappedPlant() #for plant objects
path = "../../modelparameter/structural/plant/" 
name = "fspm2023" 
pl.readParameters(path + name + ".xml")

soilSpace = pb.SDF_PlantContainer(500,500, 500, True) #to avoid root growing aboveground
pl.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil

verbose = False
pl.initialize(verbose)
pl.simulate(simtime,verbose)

kz = 4.32e-2  # axial conductivity [cm3/day]
kr_root = 1.728e-4  # radial conductivity [1/day]
kr_stem = 0  # radial conductivity of stem  [1/day], set to 0
kr_leaf =  0.004 #  radial conductivity of leaf xylem membrane [1/day] different from stomatal conductance
depth = 60 #cm
p_top = -1000  # top soil total potential [cm]
p_bot = p_top +depth # bot total potential [cm]

#envirenmental variables
RH = 0.5 # relative humidity [-]
TairC = 20 #air temperature [°C]
TairK = TairC + 273.15
Q = 900e-6 #irradiance [mol quanta m-2 s-1]
cs = 350e-6 #co2 partial pressure at leaf surface [mol mol-1]
es = 0.61078 * np.exp(17.27 * TairC / (TairC + 237.3))/10 # [hPa] #from FAO56
ea = es * RH


Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
rho_h2o = 1 #g/cm3
Mh2o = 18.05 #g/mol
MPa2hPa = 10000
hPa2cm = 1.0197                     
p_a = np.log(RH)*Rgaz*rho_h2o*TairK/Mh2o *MPa2hPa*hPa2cm #air water potential [cm]


#we give as input the initial guess of water potential and internal CO2 concentration
#for the fixed-point iteration
r = PhotosynthesisPython(pl, p_top, cs*0.7) #Equivalent to XylemFluxPython() for whole plant


" Coupling to soil "
z_ = np.linspace(0,-depth,depth) # 0, -1, -2, ... -35

p_s = np.linspace(p_top, p_bot, depth)
soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1
pl.setSoilGrid(soil_index)


""" Set up organ conductivities and atmospheric conditions"""

r.setKr([[kr_root],[kr_stem],[kr_leaf]]) 
r.setKx([[kz]])
r.airPressure = p_a 
r.Qlight = Q;


#
# In the following 'cells = True' means total soil potentials 'p_s' is given for soil voxel
#
r.solve_photosynthesis(ea_=ea, es_ = es, sim_time_ = simtime, sxx_=p_s, cells_ = True, TairC_ = TairC)
# Make sure to let the solve_photosynthesis() function finish before starting the block of code below.


plantWatPot = r.psiXyl
fluxes = r.radial_fluxes(simtime, plantWatPot, p_s,[], True)  # cm3/day
surfs = np.multiply(np.array(r.rs.segLength()), 2*np.array(r.rs.radii)*np.pi)  # root segment side surface [cm2]
fluxes = np.divide(fluxes, surfs)  # we convert to [cm3/(cm2 day)]

""" plot results """

fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    #nodes = self.get_nodes()
    y = r.get_nodes_organ_type(ndType)#coordinates
    x = np.array(plantWatPot)[r.get_nodes_index(ndType)]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (cm)")
plt.title("Xylem matric potential (cm)")
plt.show()


fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    segIdx = r.get_segments_index(ndType)
    nodesy = segIdx + np.ones(segIdx.shape, dtype = np.int64)
    nodes = r.get_nodes()
    y = np.array(nodes)[nodesy]#coordinates
    x = fluxes[segIdx]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Fluxes (cm/day)")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()


print("Transpiration", sum(fluxes[r.get_segments_index(4)]), "cm3/day")
print("Net assimilation", np.sum(r.An)*1e3, "mmol CO2 m-2 s-1")