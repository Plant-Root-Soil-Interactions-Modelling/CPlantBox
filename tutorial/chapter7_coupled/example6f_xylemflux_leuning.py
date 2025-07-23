""" water movement within the root (static soil) """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython
from functional.Leuning import Leuning

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004  #  cm3/day radial conductivity between xylem and guard cell
p_s = -200  # static water potential (saturation) 33kPa in cm
# p_g = -2000 # water potential of the guard cell
RH = 0.5  # relative humidity
TairC = 20
p_a = -1000  # default outer water potential
simtime = 14.0  # [day] for task b
k_soil = []
Q = 900e-6  # mol quanta m-2 s-1 light, example from leuning1995
cs = 350e-6  # co2 paartial pressure at leaf surface (mol mol-1)
TairK = TairC + 273.15

es = 0.61078 * np.exp(17.27 * TairC / (TairC + 237.3))
ea = es * RH
VPD = es - ea

# root system
pl = pb.MappedPlant()
path = "../../modelparameter/structural/plant/"
name = "manyleaves"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

pl.initialize()
pl.simulate(simtime, False)
# rs.simulate(simtime, False) #test to see if works in case of several simulate

r = Leuning(pl)
nodes = r.get_nodes()

r.setKr([[kr], [kr_stem], [gmax]])
r.setKx([[kz]])
r.airPressure = p_a

# Numerical solution
rx = r.solve_leuning(sim_time = simtime, sxx = [p_s], cells = True, Qlight = Q, VPD = VPD,
Tl = TairK, p_linit = p_s, ci_init = cs, cs = cs, soil_k = [], log = True, verbose = True)

fluxes = r.radial_fluxes(simtime, rx, [p_s], k_soil, True)  # cm3/day
# r.summarize_fluxes(fluxes, simtime, rx, [p_s], k_soil, True, show_matrices = False)

# plot results
fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    y = r.get_nodes_organ_type(ndType)  # coordinates
    x = rx[r.get_nodes_index(ndType)]
    ax.scatter(x, y[:, 2], c = color[ndType - 2], label = name[ndType - 2],
               alpha = 0.3, edgecolors = 'none')

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
    y = nodes[nodesy]  # coordinates
    x = fluxes[segIdx]
    ax.scatter(x, y[:, 2], c = color[ndType - 2], label = name[ndType - 2],
               alpha = 0.3, edgecolors = 'none')

ax.legend()
ax.grid(True)
plt.xlabel("Fluxes (cm3/day)")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()

# Additional vtk plot
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes", fluxes)  # cut off for vizualisation
ana.write("results/example_6f.vtp", ["radius", "surface", "rx", "fluxes"])  #
# vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
