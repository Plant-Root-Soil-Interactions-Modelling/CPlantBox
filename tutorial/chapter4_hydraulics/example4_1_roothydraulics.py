"""water movement within the root (static soil)"""

import numpy as np

import plantbox as pb
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan  # |\label{l41:imports_end}|
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l41:imports}|
import plantbox.visualisation.vtk_plot as vp

# Parameters |\label{l41:parameters}|
initial_age = 14  # root system age (day)
kx = 4.32e-2  # axial conductivity (cm3 day-1)
kr = 1.728e-4  # radial conductivity (day-1)
h_s_initial = -300  # soil total potential (cm)
h_x_collar = -1000  # dirichlet bc at top (cm) |\label{l41:h0}|
t_pot = -1  # potential plant transpiration (cm3 day-1) |\label{l41:t_pot}|

# Root system |\label{l41:rootsystem}|
plant = pb.MappedPlant()  # |\label{l41:mappedplant}|
path = "../../modelparameter/structural/rootsystem/"
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")
plant.initialize()
plant.simulate(initial_age)  # |\label{l41:rootsystem_end}|

# Root hydraulic properties
params = PlantHydraulicParameters()  # |\label{l41:hydraulicparams}|
params.set_kr_const(kr)  # day-1
params.set_kx_const(kx)  # cm3 day-1
hm = HydraulicModel_Doussan(plant, params)  # |\label{l41:model}|
# hm = HydraulicModel_Meunier(plant, params)

# Numerical solution  |\label{l41:numerical}|
ns = plant.getNumberOfMappedSegments()
h_sr = plant.total2matric(h_s_initial * np.ones((ns,)))  # |\label{l41:h_sr}|
h_x = hm.solve_dirichlet(initial_age, h_x_collar, h_sr, cells=False)  # |\label{l41:dirichlet}|
print(f"Root collar potential {h_x[0]:g} [cm], transpiration {hm.get_transpiration(initial_age, h_x, h_sr):g} (cm3 day-1)")
h_x = hm.solve_neumann(initial_age, t_pot, h_sr, cells=False)  # |\label{l41:neumann}|
print(f"Root collar potential {h_x[0]:g} [cm], transpiration {hm.get_transpiration(initial_age, h_x, h_sr):g} (cm3 day-1)")  # |\label{l41:numerical_end}|

# Additional vtk plot
ana = pb.SegmentAnalyser(hm.ms.mappedSegments())  # |\label{l41:sa}|
ana.addData("h_x", h_x)  # xylem potentials (cm)
ana.addData("SUF", hm.get_suf(initial_age))  # standard uptake fraction
ana.addAge(initial_age)  # age (day) |\label{l41:age}|
ana.addHydraulicConductivities(params, initial_age)  # kr (day-1), kx (cm3 day-1) |\label{l41:conductivities}|
ana.addFluxes(hm, h_x, h_sr, initial_age)  # "axial_flux" (cm3 day-1), "radial_flux" (cm3 cm-3 day-1) |\label{l41:fluxes}|
vp.plot_plant(ana, "radial_flux")  # |\label{l41:sa_end}|

# Output for Paraview
ana.write(
    "results/example4_1_roothydraulics.vtp",  # |\label{l41:paraview}|
    types=["radius", "subType", "age", "h_x", "SUF", "kr", "kx", "axial_flux", "radial_flux"],
)
