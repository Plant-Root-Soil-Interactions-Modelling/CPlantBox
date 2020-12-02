""" water movement within the root (static soil) """
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gs = 0.03  # radial conductivity of leaves = stomatal conductivity [1/day]
p_s = -200  # static soil pressure [cm]
p_a = -1000 #static air pressure
p0 = -500  # dircichlet bc at top
simtime = 14  # [day] for task b

""" root system """
rs = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
soil_index = lambda x, y, z : 0
rs.setSoilGrid(soil_index)
rs.initialize()
rs.simulate(3, False)
#rs.simulate(simtime, False) #test to see if works in case of several simulate


r = XylemFluxPython(rs) 
nodes = r.get_nodes()
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

r.setKr([[kr],[kr_stem],[gs]]) 
r.setKx([[kz]])
p_out = r.get_outer_matpot_matix(p_s, p_a) #create matrix to have p_a outer pressure for stem/leaves and p_s outer pressure for roots.
leavenodes = r.get_nodes_index(4) #only takes for nodes of leaves
# Numerical solution 
r.node_ind = r.get_nodes_index(4)
r.seg_ind = r.get_segments_index(4)
rx = r.solve_dirichlet(sim_time= 0., value=p0,sxc= p_s, sxx=p_out, cells=False) #water matric pot given per segment
print("Transpiration", r.collar_flux(simtime, rx, p_out),"cm3/day")

fluxes = r.segFluxes(simtime, rx, p_out, False)  # cm3/day

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

#Additional vtk plot
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes", np.maximum(fluxes, -1.e-3))  # cut off for vizualisation
ana.write("results/example_6d.vtp", ["radius", "surface", "rx", "fluxes"]) #
#vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"