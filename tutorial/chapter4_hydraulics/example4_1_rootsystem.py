""" water movement within the root (static soil) """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
initial_age = 14  # root system age [day]
kx = 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]
Hs = -300  # soil total potential [cm]
p0 = -1000  # dirichlet bc at top [cm]
t_pot = -1  # potential plant transpiration [cm3/day]

""" root system """
plant = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
plant.readParameters(path + name + ".xml")
plant.initialize()
plant.simulate(initial_age, False)

""" root hydraulic properties """
params = PlantHydraulicParameters()
params.set_radial_conductivity(kr)
params.set_axial_conductivity(kx)
r = HydraulicModel_Doussan(plant, params)  # hydraulic model

""" Numerical solution """
ns = plant.getNumberOfMappedSegments()
hsr = plant.total2matric(Hs * np.ones((ns,)))

hx = r.solve_dirichlet(initial_age, p0, hsr, cells = False)
print("Root collar potential {:g} [cm], transpiration {:g} [cm3/day]".format(hx[0], r.get_transpiration(initial_age, hx, hsr)))
# hx = r.solve_neumann(initial_age, t_pot, hsr, cells = False)
# print("Root collar potential {:g} [cm], transpiration {:g} [cm3/day]".format(hx[0], r.get_transpiration(initial_age, hx, hsr)))

fluxes = r.radial_fluxes(initial_age, hx, hsr, False)  # cm3/day

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.ms.mappedSegments())
ana.addData("hx", hx)  # xylem potentials [cm]
ana.addData("SUF", r.get_suf(initial_age))  # standard uptake fraction [1]
ana.addAge(initial_age)  # age [day]
ana.addHydraulicConductivities(params, initial_age)  # kr [1/day], kx [cm3/day]
ana.addFluxes(r, hx, hsr, initial_age)  # "axial_flux" [cm3/day], "radial_flux" [ (cm3/cm2) / day]
vp.plot_roots(ana, "radial_flux")

""" output for paraview """
ana.write("example4_1_root_hydraulics.vtp",
          types = ["radius", "subType", "age", "hx", "SUF", "kr", "kx", "axial_flux", "radial_flux"])
