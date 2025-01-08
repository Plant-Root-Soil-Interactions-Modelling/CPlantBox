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
kx = 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]
p_s = -200  # constant soil potential [cm]
p0 = -500  # dirichlet bc at top [cm]
initial_age = 14  # root system age [day]

""" root system """
rs = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(initial_age, False)

""" root problem """
params = PlantHydraulicParameters()
params.setKx([kx, kx, kx, kx, kx, kx])  # axial conductivities per root order
params.setKr([kr * 0, kr, kr , kr, kr, kr])  # radial conductivities per root order

r = HydraulicModel_Doussan(rs, params)  # hydraulic model

""" Numerical solution """
soil = [p_s]  # soil with a single soil cell
rx = r.solve_dirichlet(initial_age, p0, p_s, soil, True)
fluxes = r.segFluxes(initial_age, rx, -200 * np.ones(rx.shape), False)  # cm3/day
print("Transpiration", r.collar_flux(simtime, rx, [p_s]), "cm3/day")

""" Macroscopic root system parameter """
suf = r.get_suf(simtime)
krs, _ = r.get_krs(simtime)
print("Krs: ", krs, "cm2/day")

""" plot results """
nodes = r.get_nodes()
plt.plot(rx, nodes[:, 2] , "r*")
plt.xlabel("Xylem potentials (cm)")
plt.ylabel("Depth (m)")
plt.show()

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.rs.mappedSegments())
ana.addData("rx", rx)  # xylem potentials [cm]
ana.addData("SUF", suf)  # standard uptake fraction [1]
ana.addAge(simtime)  # age [day]
ana.addConductivities(r, simtime)  # kr [1/day], kx [cm3/day]
ana.addFluxes(r, rx, p_s * np.ones(rx.shape), simtime)  # "axial_flux" [cm3/day], "radial_flux" [ (cm3/cm2) / day]
vp.plot_roots(ana, "subType")  # "rx", "SUF", "age", kr, "axial_flux", "radial_flux"
