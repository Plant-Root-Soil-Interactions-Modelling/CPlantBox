""" 
Benchmark M3.1 Single root: steady state vertical root solved with the NEW PlantHydraulicModel and PlantHydraulicParameters classes
"""

import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import plantbox as pb

from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.PlantHydraulicParameters import PlantHydraulicParameters
import visualisation.vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
g = 9.8065 * 100.*24.*3600.*24.*3600.  # gravitational acceleration [cm day-2]
rho = 1.  # density of water, [g/cm^3]
L = 50  # length of single straight root [cm]
a = 0.2  # radius [cm] <--------------------------------------------------------- ???
kz0 = 4.32e-2  # [cm^3/day]
kz = kz0 / (rho * g)  # axial conductivity [cm^5 s / g]
kr0 = 1.728e-4  # [1/day]
kr = kr0 / (rho * g)  # radial conductivity per root type [cm^2 s / g]
p_s = -200  # static soil pressure [cm]
p0 = -1000  # dircichlet bc at top

""" Analytical solution """
c = 2 * a * np.pi * kr / kz
p_r = lambda z: p_s + d[0] * np.exp(np.sqrt(c) * z) + d[1] * np.exp(-np.sqrt(c) * z)  #

AA = np.array([[1, 1], [np.sqrt(c) * np.exp(-np.sqrt(c) * L), -np.sqrt(c) * np.exp(np.sqrt(c) * L)] ])  # # Boundary conditions dirichlet top, neumann bot
bb = np.array([p0 - p_s, -1])  # -rho * g
d = np.linalg.solve(AA, bb)  # compute constants d_1 and d_2 from bc

za_ = np.linspace(0, -L, 100)  # Evaluate function
pr = list(map(p_r, za_))

plt.plot(pr, za_)

""" Numeric solution """
N = 100  # resolutionlinearSystem
z_ = np.linspace(0., -L, N)

nodes, segs, radii = [], [], []
for z in z_:
    nodes.append(pb.Vector3d(0, 0, z))
for s in range(0, N - 1):
    segs.append(pb.Vector2i(s, s + 1))
    radii.append(a)

# print("dx", np.diff(z_))

rs = pb.MappedSegments(nodes, segs, radii)
soil_index = lambda x, y, z: 0
rs.setSoilGrid(soil_index)

params = PlantHydraulicParameters(rs)
params.set_kr_const(kr0)
params.set_kx_const(kz0)

# r = HydraulicModel_Doussan(rs, params, cached = False)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
r = HydraulicModel_Meunier(rs, params, cached = False)
kx = params.getKx(0.)
# print(kx)
rx = r.solve_dirichlet(0., p0, [p_s], cells = True)
print("rx[0]", rx[0], "rx[1]", rx[1])
print("rx", rx.shape)

trans = r.get_transpiration(0., rx, [p_s], cells = True)
print("Transpiration", trans, "cm3/day")
plt.plot(rx, z_, "r*")

# get_collar_potential(t_act, rsx)

# #
# # check net fluxes
# #
# simtime = 0.
# radial_fluxes = r.radial_fluxes(simtime, rx, [p_s])
# axial_fluxes = r.axial_fluxes(simtime, rx, [p_s])
# # axial_i = np.array([r.axial_flux(i, simtime, rx, [p_s], True, True) for i in range(0, len(axial_fluxes))])  # same as axial_fluxes, but in node j
# # axial_j = np.array([r.axial_flux(i, simtime, rx, [p_s], True, False) for i in range(0, len(axial_fluxes))])  # same as axial_fluxes, but in node j
# ana = pb.SegmentAnalyser(r.ms)
# ana.addData("rx", rx)  # node data are converted to segment data
# ana.addData("radial", radial_fluxes)
# ana.addData("axial", axial_fluxes)
# # ana.addData("net", axial_i - axial_j - radial_fluxes)
# pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "length", "rx", "axial", "radial", "net"])
# vp.plot_roots(pd, "radial")  # axial, radial, rx

rx = r.solve_neumann(0., -2., [p_s], cells = True)  # or solve ...
# rx = r.solve_dirichlet(0., -869.6363009384586, [p_s], cells = True)  # or solve .

trans = r.get_transpiration(0., rx, [p_s], cells = True)
trans2 = r.axial_fluxes(0., rx, [p_s], cells = True)

print("Transpiration", trans, trans2[0], trans - trans2[0], "cm3/day")
plt.plot(rx, z_, "g*")

plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (m)")
plt.legend(["analytic solution", "numeric solution", "predescribed flux -2 cm$^3$ day$^{-1}$"])
plt.show()

