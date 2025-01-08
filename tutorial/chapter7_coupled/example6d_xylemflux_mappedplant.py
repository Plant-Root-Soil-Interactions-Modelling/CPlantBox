""" water movement within the root (static soil) """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gs = 0.03  # radial conductivity of leaves = stomatal conductivity [1/day]
p_s = -200  # static soil pressure [cm]
p_a = -100000  # static air pressure
# p0 = -500  # dircichlet bc at top
simtime = 14.0  # [day] for task b
k_soil = []

""" root system """
rs = pb.MappedPlant()
path = "../../modelparameter/structural/plant/"
name = "manyleaves"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
soil_index = lambda x, y, z: 0
rs.setSoilGrid(soil_index)
rs.initialize()
rs.simulate(3, False)
# rs.simulate(simtime, False) #test to see if works in case of several simulate

r = XylemFluxPython(rs)
nodes = r.get_nodes()
tiproots, tipstem, tipleaf = r.get_organ_nodes_tips()  # end node of end segment of each organ
node_tips = np.concatenate((tiproots, tipstem, tipleaf))
tiproots, tipstem, tipleaf = r.get_organ_segments_tips()  # end segment of each organ
seg_tips = np.concatenate((tiproots, tipstem, tipleaf))
"""
    we can give a kr/kx:
        constant across type: r.setKx([[kz]])
        by type: r.setKx([[kz], [kz2], [kz3]])
        by type and subtype: r.setKx([[kz, kz2], [kza, kzb], [kzd, kzf]])
    att: the bud (from which the leaf grows) has to have a kz
        one of the above + time dependant:
            constant across type: r.setKx([[kz1, kz2, kr3]], [[age1, age2, age3]])
            by type: r.setKxTables([[[kz1, kz2, kr3],[kza, kzb, krc]]], [[[age1, age2, age3], [agea, ageb, agec]]])
            by type and subtype: r.setKxTables([[[kz1, kz2, kr3],[kza, kzb, krc]],[[kz1, kz2, kr3],[kza, kzb, krc]]],
                    [[[age1, age2, age3], [agea, ageb, agec]],[[age1, age2, age3], [agea, ageb, agec]]])
"""

r.setKr([[kr], [kr_stem], [gs]])
r.setKx([[kz]])
r.airPressure = p_a

# Numerical solution
r.neumann_ind = seg_tips  # segment indices for Neumann b.c.
r.node_ind = node_tips
rx = r.solve_neumann(sim_time = simtime, value = 0, sxx = [p_s], cells = True)  # water matric pot given per segment
fluxes = r.radial_fluxes(simtime, rx, [p_s], k_soil, True)  # cm3/day
# too slow
# r.summarize_fluxes(fluxes, simtime, rx, [p_s], k_soil, True)#r.segFluxes(simtime, rx, p_out, False)  # cm3/day

# plot results
plt.plot(rx, nodes[:, 2] , "r*")
plt.xlabel("Xylem pressure (cm)")
plt.ylabel("Depth (cm)")
plt.title("Xylem matric potential (cm)")
plt.show()

plt.plot(fluxes, nodes[1:, 2] , "r*")
plt.xlabel("Fluxes")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()

# Additional vtk plot
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes", np.maximum(fluxes, -1.e-3))  # cut off for vizualisation
ana.write("results/example_6d.vtp", ["radius", "surface", "rx", "fluxes"])  #
# vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
