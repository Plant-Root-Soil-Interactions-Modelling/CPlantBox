""" water movement within the root (static soil) """
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
from xylem_flux import XylemFluxPython  # Python hybrid solver
from Leuning import Leuning
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day] 
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.0864 #  cm3/day radial conductivity of leaves = stomatal conductivity [1/day]
#p_s = -200  # static water potential (saturation) 33kPa in cm
p_a =  -1000  #static air water potential 
simtime = 14.0  # [day] for task b
k_soil = []

# root system 
pl = pb.MappedPlant() #pb.MappedRootSystem() #pb.MappedPlant()
path = "../../../modelparameter/plant/" #"../../../modelparameter/rootsystem/" 
name = "manyleaves" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")
RH = np.arange(0.1, 0.9, 0.01)
Q = np.arange(0, 1000e-6, 10e-6)#[250e-6,300e-6, 350e-6, 400e-6,450e-6, 500e-6,750e-6, 1000e-6] # mol quanta m-2 s-1 light, example from leuning1995
TairC = np.arange(-30,30, 5)
p_s = np.arange(-200, -45000, -100) #cm
N = np.arange(0.1, 6, 0.5) #%
es = 0.61078 * np.exp(17.27 * 20 / (20+ 237.3))
VPDvar = [es - es*rh for rh in RH]

cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments


pl.initialize()
pl.simulate(simtime, False)
#rs.simulate(simtime, False) #test to see if works in case of several simulate


r = Leuning(pl) 
nodes = r.get_nodes()
tiproots, tipstem, tipleaf = r.get_organ_nodes_tips() #end node of end segment of each organ
node_tips = np.concatenate((tiproots, tipstem, tipleaf))
tiproots, tipstem, tipleaf = r.get_organ_segments_tips() #end segment of each organ
seg_tips = np.concatenate((tiproots, tipstem, tipleaf))


r.setKr([[kr],[kr_stem],[gmax]]) #gmax will be changed by the leuning function 
r.setKx([[kz]])
r.airPressure = p_a

# Numerical solution 
r.seg_ind = seg_tips # segment indices for Neumann b.c.
r.node_ind = node_tips
leaf_nodes = r.get_nodes_index(4)
variables = [ Q,RH, TairC, p_s, N]
results=[[],[],[], [], []]
resultsAn=[[],[],[], [], []]
resultsgco2=[[],[],[], [], []]
resultsVc=[[],[],[], [], []]
resultsVj=[[],[],[], [], []]
resultscics=[[],[],[], [], []]
resultsfw=[[],[],[], [], []]
resultspl=[[],[],[], [], []]
for i in range(len(variables)):
    Q_input = 900e-6
    RH_input = 0.5
    Tair_input = 20 
    p_s_input = -200
    N_input = 4.4
    var = [Q_input,RH_input, Tair_input, p_s_input, N_input]
    trans = []
    An = []
    Vc = []
    Vj = []
    gco2 = []
    cics = []
    fw = []
    pl = []
    #VPDvar  = []
    for j in range(len(variables[i])):
        print(i, j)
        var[i] = variables[i][j]
        es = 0.61078 * np.exp(17.27 * var[2] / (var[2] + 237.3)) #2.338205 kPa
        ea = es * var[1] #1.169102 kPa
        VPD = es - ea # 1.169103 kPa
        #if(j == 1 or j ==2):
         #   VPDvar.append(VPD)
        rx = r.solve_leuning( sim_time = simtime,sxx=[var[3]], cells = True, Qlight = var[0],VPD = VPD,Tl = var[2] + 273.15,p_linit = var[3],
        ci_init = cs*0.7, cs =cs, soil_k = [], N = var[4], log = False)
        fluxes = r.radial_fluxes(simtime, rx, [var[3]], k_soil, True)  # cm3/day
        organTypes = np.array(r.rs.organTypes)
        trans.append(sum(np.where(organTypes == 4, fluxes,0)))
        An.append(np.mean(r.An)*1e6)
        Vc.append(np.mean(r.Vc)*1e6)
        Vj.append(np.mean(r.Vj)*1e6)
        gco2.append(np.mean(r.gco2))
        cics.append(np.mean(r.ci)/cs)
        fw.append(np.mean(r.fw))
        pl.append(np.mean(r.x[leaf_nodes]))
    results[i] = trans
    resultsAn[i] = An
    resultsVc[i] = Vc
    resultsVj[i] = Vj
    resultsgco2[i] = gco2
    resultscics[i] = cics
    resultsfw[i] = fw
    resultspl[i] = pl

logfile = open('leuning6finput.txt', "w")
logfile.write(repr(variables))
logfile.close()
logfile = open('leuning6fE.txt', "w")
logfile.write(repr(results))
logfile.close()
logfile = open('leuning6fAn.txt', "w")
logfile.write(repr(resultsAn))
logfile.close()
logfile = open('leuning6fVc.txt', "w")
logfile.write(repr(resultsVc))
logfile.close()
logfile = open('leuning6fVj.txt', "w")
logfile.write(repr(resultsVj))
logfile.close()
logfile = open('leuning6fco2.txt', "w")
logfile.write(repr(resultsgco2))
logfile.close()
logfile = open('leuning6fcics.txt', "w")
logfile.write(repr(resultscics))
logfile.close()
logfile = open('leuning6fpl.txt', "w")
logfile.write(repr(resultspl))
logfile.close()
logfile = open('leuning6ffw.txt', "w")
logfile.write(repr(resultsfw))
logfile.close()

# plot results 
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], results[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='total E (cm3/d)')
axs[0, 1].plot(variables[1], results[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, results[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], results[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='total E (cm3/d)')
axs[1, 1].plot(variables[3], results[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], results[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsAn[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean An (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variables[1], resultsAn[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultsAn[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultsAn[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean An (μmol CO2 m-2 s-1)')
axs[1, 1].plot(variables[3], resultsAn[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultsAn[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsVc[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean Vc (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variables[1], resultsVc[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultsVc[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultsVc[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean Vc (μmol CO2 m-2 s-1)')
axs[1, 1].plot(variables[3], resultsVc[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultsVc[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsVj[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean Vj (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variables[1], resultsVj[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultsVj[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultsVj[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean Vj (μmol CO2 m-2 s-1)')
axs[1, 1].plot(variables[3], resultsVj[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultsVj[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsgco2[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean gco2 (mol CO2 m-2 s-1)')
axs[0, 1].plot(variables[1], resultsgco2[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultsgco2[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultsgco2[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean gco2 (mol CO2 m-2 s-1)')
axs[1, 1].plot(variables[3], resultsgco2[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultsgco2[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultscics[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean ci/cs (-)')
axs[0, 1].plot(variables[1], resultscics[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultscics[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultscics[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean ci/cs (-)')
axs[1, 1].plot(variables[3], resultscics[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultscics[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsfw[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean fw (-)')
axs[0, 1].plot(variables[1], resultsfw[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultsfw[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultsfw[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean fw (-)')
axs[1, 1].plot(variables[3], resultsfw[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultsfw[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultspl[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean Ψleaf (cm)')
axs[0, 1].plot(variables[1], resultspl[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(VPDvar, resultspl[1], 'tab:pink')
axs[0, 2].set(xlabel='VPD (kPa)')
axs[1, 0].plot(variables[2], resultspl[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean Ψleaf (cm)')
axs[1, 1].plot(variables[3], resultspl[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultspl[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
"""
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
plt.xlabel("Fluxes (cm3/day)")
plt.ylabel("Depth (cm)")
plt.title("water fluxes")
plt.show()

#Additional vtk plot
ana = pb.SegmentAnalyser(r.rs)
ana.addData("rx", rx)
ana.addData("fluxes",fluxes)  # cut off for vizualisation
ana.write("results/example_6e.vtp", ["radius", "surface", "rx", "fluxes"]) #
#vp.plot_roots(ana, "rx", "Xylem matric potential (cm)")  # "fluxes"
"""