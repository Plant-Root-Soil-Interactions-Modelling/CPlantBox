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
results=[[],[]]
resultsAn=[[],[]]
resultsgco2=[[],[]]
resultsVc=[[],[]]
resultsVj=[[],[]]
resultscics=[[],[]]
resultsfw=[[],[]]
resultspl=[[],[]]
Q_input = 900e-6
RH_input = 0.5
Tair_input = 20 
p_s_input = -200
N_input = 4.4
cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)

cics_ = np.arange(0, 1, 0.1) 
pl_ = np.arange(0, -20000, -100) 
variablesinit = [ pl_, cics_]
for i in range(len(variablesinit)):
    plinit = -200
    cicsinit = 0.7
    varinit = [plinit, cicsinit]
    trans = []
    An = []
    Vc = []
    Vj = []
    gco2 = []
    cics = []
    fw = []
    pl = []
    for j in range(len(variablesinit[i])):
        print(i, j)
        varinit[i] = variablesinit[i][j]
        print(varinit[0],cs*varinit[1])
        es = 0.61078 * np.exp(17.27 * Tair_input / (Tair_input + 237.3)) #2.338205 kPa
        ea = es * RH_input 
        VPD = es - ea 
        rx = r.solve_leuning( sim_time = simtime,sxx=[p_s_input], cells = True, Qlight = Q_input,VPD = VPD,Tl = Tair_input + 273.15,p_linit = varinit[0],
        ci_init = cs*varinit[1], cs =cs, soil_k = [], N = N_input, log = False)
        fluxes = r.radial_fluxes(simtime, rx, [p_s_input], k_soil, True)  # cm3/day
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
logfile.write(repr(variablesinit))
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

fig, axs = plt.subplots(2,2)
axs[0, 0].plot(variablesinit[0], resultsAn[0])
axs[0, 0].set(xlabel='p_linit', ylabel='mean An (μmol CO2 m-2 s-1)')
axs[1, 0].plot(variablesinit[0], resultsVc[0], 'tab:red')
axs[1, 0].set(xlabel='p_linit', ylabel='mean Vc (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variablesinit[0], resultsVj[0], 'tab:brown')
axs[0, 1].set(xlabel='p_linit', ylabel='mean Vj (μmol CO2 m-2 s-1)')
axs[1, 1].plot(variablesinit[0], resultsgco2[0], 'tab:brown')
axs[1, 1].set(xlabel='p_linit', ylabel='mean gco2 mean Vc (μmol CO2 m-2 s-1)')
plt.show()
fig, axs = plt.subplots(2,2)
axs[0, 0].plot(variablesinit[0], results[0])
axs[0, 0].set(xlabel='p_linit', ylabel='E')
axs[0, 1].plot(variablesinit[0], resultsfw[0], 'tab:brown')
axs[0, 1].set(xlabel='p_linit', ylabel='fw')
axs[1, 0].plot(variablesinit[0], resultspl[0], 'tab:brown')
axs[1, 0].set(xlabel='p_linit', ylabel='pl')
axs[1, 1].plot(variablesinit[0], resultscics[0], 'tab:brown')
axs[1, 1].set(xlabel='p_linit', ylabel='ci/cs (-)')
plt.show()
fig, axs = plt.subplots(2,2)
axs[0, 0].plot(variablesinit[1], resultsAn[1])
axs[0, 0].set(xlabel='ci/cs init', ylabel='mean An (μmol CO2 m-2 s-1)')
axs[1, 0].plot(variablesinit[1], resultsVc[1], 'tab:red')
axs[1, 0].set(xlabel='ci/cs init', ylabel='mean Vc (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variablesinit[1], resultsVj[1], 'tab:brown')
axs[0, 1].set(xlabel='ci/cs init', ylabel='mean Vj (μmol CO2 m-2 s-1)')
axs[1, 1].plot(variablesinit[1], resultsgco2[1], 'tab:brown')
axs[1, 1].set(xlabel='ci/cs init', ylabel='mean gco2 mean Vc (μmol CO2 m-2 s-1)')
plt.show()
fig, axs = plt.subplots(2,2)
axs[0, 0].plot(variablesinit[1], results[1])
axs[0, 0].set(xlabel='ci/cs init', ylabel='E')
axs[0, 1].plot(variablesinit[1], resultsfw[1], 'tab:brown')
axs[0, 1].set(xlabel='ci/cs init', ylabel='fw')
axs[1, 0].plot(variablesinit[1], resultspl[1], 'tab:brown')
axs[1, 0].set(xlabel='ci/cs init', ylabel='pl')
axs[1, 1].plot(variablesinit[1], resultscics[1], 'tab:brown')
axs[1, 1].set(xlabel='ci/cs init', ylabel='ci/cs (-)')
plt.show()