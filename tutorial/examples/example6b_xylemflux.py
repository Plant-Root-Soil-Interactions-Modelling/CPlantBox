""" water movement within the root (static soil) """

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Meunier
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]
p_s = -200  # constant soil potential [cm]
p0 = -500  # dirichlet bc at top [cm]
simtime = 14  # [day]

""" root system """
rs = pb.MappedPlant()
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(simtime, False)

""" root problem """
params = PlantHydraulicParameters(rs)
for sub_type in range(6):
    params.set_kx_const(kz, subType=sub_type)
    params.set_kr_const(0.0 if sub_type == 0 else kr, subType=sub_type)

r = HydraulicModel_Meunier(rs, params, cached=False)
soil_index = lambda x, y, z: 0  # maps every coordinate to soil cell 0
r.ms.setSoilGrid(soil_index)

""" Numerical solution """
soil = [p_s]  # soil with a single soil cell
rx = r.solve_dirichlet(simtime, p0, soil, cells=True)
fluxes = r.radial_fluxes(simtime, rx, soil, cells=True)  # cm3/day
print("Transpiration", r.get_transpiration(simtime, rx, soil, cells=True), "cm3/day")

""" Macroscopic root system parameter """
suf = r.get_suf(simtime)
krs, _ = r.get_krs(simtime)
print("Krs: ", krs, "cm2/day")

""" plot results """
nodes = r.get_nodes()
plt.plot(rx, nodes[:, 2], "r*")
plt.xlabel("Xylem potentials (cm)")
plt.ylabel("Depth (cm)")
plt.show()

""" Additional vtk plot """
ana = pb.SegmentAnalyser(r.ms)
ana.addData("rx", rx)  # xylem potentials [cm]
ana.addData("SUF", suf)  # standard uptake fraction [1]
ana.addAge(simtime)  # age [day]
ana.addData("kr", r.get_kr(simtime))  # [1/day]
ana.addData("kx", r.get_kx(simtime))  # [cm3/day]
ana.addData("axial_flux", r.axial_fluxes(simtime, rx, soil, cells=True))  # [cm3/day]
ana.addData("radial_flux", fluxes)  # [cm3/day]
vp.plot_roots(ana, "subType")  # "rx", "SUF", "age", "kr", "kx", "axial_flux", "radial_flux"
