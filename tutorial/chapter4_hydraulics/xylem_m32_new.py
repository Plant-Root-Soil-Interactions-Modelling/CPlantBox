import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.PlantHydraulicParameters import PlantHydraulicParameters

import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import numpy as np

""" 
Benchmark M3.2 Root system: steady state small root system solved with the Python/cpp Hybrid solver
(does not work in parallel)
"""
fig, (ax1, ax2) = plt.subplots(1, 2)

""" Parameters """
kz = 4.32e-2  # [cm^3/day]
kr = 1.728e-4  # [1/day]
p_s = -200  # static soil pressure [cm]
p0 = -500  # dircichlet bc at top
simtime = 14  # [day] for task b

""" root problem """
param = PlantHydraulicParameters()
param.set_kr_const(kr)
param.set_kx_const(kz)

r = HydraulicModel_Meunier("../../grids/RootSystem.rsml", param, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
r.test()

nodes = r.get_nodes()
soil_index = lambda x, y, z: 0
r.ms.setSoilGrid(soil_index)

""" Numerical solution (a) """
rx_a = r.solve_dirichlet(0., p0, [p_s], cells = True)
print("Transpiration", r.get_transpiration(0., rx_a, [p_s], cells = True), "cm3/day")
np.savetxt("results_m32a", np.vstack((nodes[:, 2], rx_a)), delimiter = ',')

ax1.plot(rx_a, nodes[:, 2] , "r*")
ax1.set_xlabel("Xylem pressure (cm)")
ax1.set_ylabel("Depth (m)")
ax1.set_title("Constant conductivities")

print()

""" Numerical solution (b) """
kr0 = np.array([[0, 1.14e-03], [2, 1.09e-03], [4, 1.03e-03], [6, 9.83e-04], [8, 9.35e-04], [10, 8.90e-04], [12, 8.47e-04], [14, 8.06e-04], [16, 7.67e-04], [18, 7.30e-04], [20, 6.95e-04], [22, 6.62e-04], [24, 6.30e-04], [26, 5.99e-04], [28, 5.70e-04], [30, 5.43e-04], [32, 5.17e-04]])
kr1 = np.array([[0, 4.11e-03], [1, 3.89e-03], [2, 3.67e-03], [3, 3.47e-03], [4, 3.28e-03], [5, 3.10e-03], [6, 2.93e-03], [7, 2.77e-03], [8, 2.62e-03], [9, 2.48e-03], [10, 2.34e-03], [11, 2.21e-03], [12, 2.09e-03], [13, 1.98e-03], [14, 1.87e-03], [15, 1.77e-03], [16, 1.67e-03], [17, 1.58e-03]])
param.set_kr_age(kr0[:, 0], kr0[:, 1], subType = 1)
param.set_kr_age(kr1[:, 0], kr1[:, 1], subType = [2, 3])

kx0 = np.array([[0, 6.74e-02], [2, 7.48e-02], [4, 8.30e-02], [6, 9.21e-02], [8, 1.02e-01], [10, 1.13e-01], [12, 1.26e-01], [14, 1.40e-01], [16, 1.55e-01], [18, 1.72e-01], [20, 1.91e-01], [22, 2.12e-01], [24, 2.35e-01], [26, 2.61e-01], [28, 2.90e-01], [30, 3.21e-01], [32, 3.57e-01]])
kx1 = np.array([[0, 4.07e-04], [1, 5.00e-04], [2, 6.15e-04], [3, 7.56e-04], [4, 9.30e-04], [5, 1.14e-03], [6, 1.41e-03], [7, 1.73e-03], [8, 2.12e-03], [9, 2.61e-03], [10, 3.21e-03], [11, 3.95e-03], [12, 4.86e-03], [13, 5.97e-03], [14, 7.34e-03], [15, 9.03e-03], [16, 1.11e-02], [17, 1.36e-02]])
param.set_kx_age(kx0[:, 0], kx0[:, 1], subType = 1)
param.set_kx_age(kx1[:, 0], kx1[:, 1], subType = [2, 3])

rx_b = r.solve_dirichlet(simtime, p0, [p_s], cells = True)
print("Transpiration", r.get_transpiration(simtime, rx_b, [p_s], cells = True), "cm3/day")
np.savetxt("results_m32b", np.vstack((nodes[:, 2], rx_b)), delimiter = ',')

ax2.plot(rx_b, nodes[:, 2] , "r*")
ax2.set_xlabel("Xylem pressure (cm)")
ax2.set_ylabel("Depth (m)")
ax2.set_title("Age dependent conductivities")

plt.show()

# """ Additional vtk plot """
# radial_fluxes = r.radial_fluxes(simtime, rx_b, [p_s])
# axial_fluxes = r.axial_fluxes(simtime, rx_b, [p_s]) # = axial_i
#
# axial_i = np.array([r.axial_flux(i, simtime, rx_b, [p_s], [], True, True) for i in range(0, len(axial_fluxes))]) # same as axial_fluxes, but in node j
# axial_j = np.array([r.axial_flux(i, simtime, rx_b, [p_s], [], True, False) for i in range(0, len(axial_fluxes))]) # same as axial_fluxes, but in node j
#
# ana = pb.SegmentAnalyser(r.rs)
# ana.addData("rx", rx_b)  # node data are converted to segment data
# ana.addData("radial", radial_fluxes)
# ana.addData("axial", axial_fluxes)
# ana.addData("net", axial_i-axial_j-radial_fluxes) # np.maximum(1.e-2,)
# pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "rx", "axial", "radial", "net"])
# vp.plot_roots(pd, "net") # axial, radial, rx

