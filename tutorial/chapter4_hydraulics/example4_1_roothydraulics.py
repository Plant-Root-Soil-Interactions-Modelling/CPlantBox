""" water movement within the root (static soil) """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l41:imports}|
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier  # |\label{l41:imports_end}|

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """  # |\label{l41:parameters}|
initial_age = 14  # root system age [day]
kx = 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]
Hs = -300  # soil total potential [cm]
h0 = -1000  # dirichlet bc at top [cm] |\label{l41:h0}|
t_pot = -1  # potential plant transpiration [cm3/day] |\label{l41:t_pot}|

""" root system """  # |\label{l41:rootsystem}|
plant = pb.MappedPlant()  # |\label{l41:mappedplant}|
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + name + ".xml")
plant.initialize()
plant.simulate(initial_age)  # |\label{l41:rootsystem_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l41:hydraulicparams}|
params.set_kr_const(kr)
params.set_kx_const(kx)
hm = HydraulicModel_Doussan(plant, params)  # |\label{l41:model}|
# hm = HydraulicModel_Meunier(plant, params)

""" Numerical solution """  # |\label{l41:numerical}|
ns = plant.getNumberOfMappedSegments()
hsr = plant.total2matric(Hs * np.ones((ns,)))  # |\label{l41:hsr}|
hx = hm.solve_dirichlet(initial_age, h0, hsr, cells = False)  # |\label{l41:dirichlet}|
print("Root collar potential {:g} [cm], transpiration {:g} [cm3/day]".format(hx[0], hm.get_transpiration(initial_age, hx, hsr)))
hx = hm.solve_neumann(initial_age, t_pot, hsr, cells = False)  # |\label{l41:neumann}|
print("Root collar potential {:g} [cm], transpiration {:g} [cm3/day]".format(hx[0], hm.get_transpiration(initial_age, hx, hsr)))  # |\label{l41:numerical_end}|

""" Additional vtk plot """
ana = pb.SegmentAnalyser(hm.ms.mappedSegments())  # |\label{l41:sa}|
ana.addData("hx", hx)  # xylem potentials [cm]
ana.addData("SUF", hm.get_suf(initial_age))  # standard uptake fraction [1]
ana.addAge(initial_age)  # age [day] |\label{l41:age}|
ana.addHydraulicConductivities(params, initial_age)  # kr [1/day], kx [cm3/day] |\label{l41:conductivities}|
ana.addFluxes(hm, hx, hsr, initial_age)  # "axial_flux" [cm3/day], "radial_flux" [ (cm3/cm2) / day] |\label{l41:fluxes}|
vp.plot_plant(ana, "radial_flux")  # |\label{l41:sa_end}|

""" output for paraview """
ana.write("restuls/example4_1_roothydraulics.vtp",  # |\label{l41:paraview}|
          types = ["radius", "subType", "age", "hx", "SUF", "kr", "kx", "axial_flux", "radial_flux"])

