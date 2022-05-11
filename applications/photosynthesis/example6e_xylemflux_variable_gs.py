""" water movement within the root (static soil) """
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from StomataModel import StomataModel
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.0864 #  cm3/day radial conductivity of leaves = stomatal conductivity [1/day]
p_s = -200  # static water potential (saturation) 33kPa in cm
p_a =  -1000  #static air water potential 
simtime = 14.0  # [day] for task b
k_soil = []

# root system 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")


""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments


pl.initialize()
pl.simulate(simtime, False)
#rs.simulate(simtime, False) #test to see if works in case of several simulate


r = StomataModel(pl,PAR = 502, VPD= 0.03,TH=50, TL=10,  Topt=28, psi1=1.1, psi2=5, gmax =gmax) 
nodes = r.get_nodes()
tiproots, tipstem, tipleaf = r.get_organ_nodes_tips() #end node of end segment of each organ
node_tips = np.concatenate((tiproots, tipstem, tipleaf))
tiproots, tipstem, tipleaf = r.get_organ_segments_tips() #end segment of each organ
seg_tips = np.concatenate((tiproots, tipstem, tipleaf))


r.setKr([[kr],[kr_stem],[gmax]]) 
r.setKx([[kz]])
r.airPressure = p_a

# Numerical solution 
r.seg_ind = seg_tips # segment indices for Neumann b.c.
r.node_ind = node_tips
rx = r.solve_neumann_gs( sim_time = simtime,sxx=[p_s], cells = True, PAR=350,VPD = 10,Tair=20,p_linit = p_s,  soil_k = [])

fluxes = r.radial_fluxes(simtime, rx, [p_s], k_soil, True)  # cm3/day
r.summarize_fluxes(fluxes, simtime, rx, [p_s], k_soil, True, show_matrices = False)


# plot results 
fig, ax = plt.subplots()
name = ["root", "stem", "leaf"]
color = ['tab:blue', 'tab:orange', 'tab:green']
for ndType in [2, 3, 4]:
    y = r.get_nodes_organ_type(ndType)#coordinates
    x = rx[r.get_nodes_index(ndType)]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

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
    y = nodes[nodesy]#coordinates
    x = fluxes[segIdx]
    ax.scatter(x, y[:,2], c=color[ndType-2],  label=name[ndType-2],
               alpha=0.3, edgecolors='none')

ax.legend()
ax.grid(True)
plt.xlabel("Fluxes")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()

#Additional vtk plot
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes",fluxes)  # cut off for vizualisation
ana.write("results/example_6e.vtp", ["radius", "surface", "rx", "fluxes"]) #
#vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
