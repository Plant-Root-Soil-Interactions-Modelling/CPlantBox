""" water movement within the root (static soil) """
import sys; sys.path.append("../../..")
from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-2  # axial conductivity [cm^3/day] => isn't it the conductance? (according to Meunier2017)
kr = 1.728e-4  # radial conductivity [1/day]
p_s = -200  # static soil pressure [cm]
p0init = -500  # dircichlet bc at top
simtime = 1  # [day] for task b

""" root system """
rs = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
rs.readParameters(path + name + ".xml")
soil_index = lambda x, y, z : 0
rs.setSoilGrid(soil_index)
rs.initialize()
rs.simulate(2, False)
rs.simulate(2, False) #test to see if works in case of several simulate


r = XylemFluxPython(rs) #just for root
nodes = r.get_nodes() #change: only takes for nodes of roots

r = XylemFluxPython(rs) #just for root
r.setKr([kr]) 
r.setKx([kz])


leavenodes = r.get_nodes_index(4) #only takes for nodes of leaves
p0 = np.full(len(leavenodes),p0init)
# Numerical solution 
rx = r.solve_dirichlet(sim_time= 0., value=p0,n0 = leavenodes,sxc= p_s, sxx=[p_s], cells=True)
print("Transpiration", r.collar_flux(simtime, rx, [p_s]),"cm3/day") #transpiration>0

fluxes = r.segFluxes(simtime, rx, -200 * np.ones(rx.shape), False)  # cm3/day

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
ana.write("results/example_6e.vtp", ["radius", "surface", "rx", "fluxes"]) #
#vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"