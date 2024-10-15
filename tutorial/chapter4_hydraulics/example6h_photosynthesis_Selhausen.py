""" Example of the photosynthesis module, using real data from the Selhausen lysimeter setup """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython
from functional.Leuning import Leuning

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np

""" Parameters """
kz = 4.32e-1  # axial conductivity [cm^3/day]
kr = 1.728e-4  # radial conductivity of roots [1/day]
kr_stem = 1.e-20  # radial conductivity of stem  [1/day], set to almost 0
gmax = 0.004  #  cm3/day radial conductivity of leaves = stomatal conductivity [1/day]
p_a = -1000  # static air water potential
simtime = 14.0  # [day] for task b
k_soil = []
plotResults = True
saveResults = False

t_init = 70
t_end = 90

# root system
pl = pb.MappedPlant()
path = "../../modelparameter/structural/plant/"
name = "manyleaves"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
pl.readParameters(path + name + ".xml")

# load data
df = pd.read_csv('../../modelparameter/functional/climate/Selhausen_weather_data.txt', delimiter = "\t")

""" soil """
min_ = np.array([-5, -5, -15])
max_ = np.array([9, 4, 0])
res_ = np.array([5, 5, 5])
pl.setRectangularGrid(pb.Vector3d(min_), pb.Vector3d(max_), pb.Vector3d(res_), True)  # cut and map segments

pl.initialize()
pl.simulate(simtime, False)

r = Leuning(pl)
r.setKr([[kr], [kr_stem], [gmax]])  # gmax will be changed by the leuning function
r.setKx([[kz]])
r.airPressure = p_a
leaf_nodes = r.get_nodes_index(4)

# Numerical solution
results = []
resultsAn = []
resultsgco2 = []
resultsVc = []
resultsVj = []
resultscics = []
resultsfw = []
resultspl = []

for i in range(t_init, t_end):

    Q_input = df['PAR'][i]
    RH_input = df['RH'][i]
    Tair_input = df['Tair'][i]
    p_s_input = df['ps'][i]
    N_input = 4.4  # nitrogen satisfaction for small wheat plants
    cs_input = df['co2'][i]
    var = [Q_input, RH_input, Tair_input, p_s_input, N_input, cs_input]
    es = 0.61078 * np.exp(17.27 * var[2] / (var[2] + 237.3))  # FAO56
    ea = es * var[1]
    VPD = es - ea
    r.Param['Patm'] = df['Pair'][i]
    rx = r.solve_leuning(sim_time = simtime, sxx = [var[3]], cells = True, Qlight = var[0], VPD = VPD, Tl = var[2] + 273.15, p_linit = var[3],
    ci_init = var[5] * 0.7, cs = var[5], soil_k = [], N = var[4], log = False, verbose = False)
    fluxes = r.radial_fluxes(simtime, rx, [var[3]], k_soil, True)  # cm3/day
    organTypes = np.array(r.rs.organTypes)
    results.append(sum(np.where(organTypes == 4, fluxes, 0)))
    resultsAn.append(np.mean(r.An) * 1e6)
    resultsVc.append(np.mean(r.Vc) * 1e6)
    resultsVj.append(np.mean(r.Vj) * 1e6)
    resultsgco2.append(np.mean(r.gco2))
    resultscics.append(np.mean(r.ci) / var[5])
    resultsfw.append(np.mean(r.fw))
    resultspl.append(np.mean(r.x[leaf_nodes]))

if saveResults:
    logfile = open('results/leuning6hE.txt', "w")
    logfile.write(repr(results))
    logfile.close()
    logfile = open('results/leuning6hAn.txt', "w")
    logfile.write(repr(resultsAn))
    logfile.close()
    logfile = open('results/leuning6hVc.txt', "w")
    logfile.write(repr(resultsVc))
    logfile.close()
    logfile = open('results/leuning6hVj.txt', "w")
    logfile.write(repr(resultsVj))
    logfile.close()
    logfile = open('results/leuning6hco2.txt', "w")
    logfile.write(repr(resultsgco2))
    logfile.close()
    logfile = open('results/leuning6hcics.txt', "w")
    logfile.write(repr(resultscics))
    logfile.close()
    logfile = open('results/leuning6hpl.txt', "w")
    logfile.write(repr(resultspl))
    logfile.close()
    logfile = open('results/leuning6hfw.txt', "w")
    logfile.write(repr(resultsfw))
    logfile.close()

timePlot = df[t_init:t_end]['time']
if plotResults:
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(timePlot, resultsAn)
    axs[0, 0].set(xlabel = '', ylabel = 'mean An (μmol CO2 m-2 s-1)')
    axs[0, 0].xaxis.set_major_locator(MaxNLocator(5))
    axs[1, 0].plot(timePlot, resultsVc, 'tab:red')
    axs[1, 0].set(xlabel = 'time', ylabel = 'mean Vc (μmol CO2 m-2 s-1)')
    axs[1, 0].xaxis.set_major_locator(MaxNLocator(5))
    axs[0, 1].plot(timePlot, resultsVj, 'tab:brown')
    axs[0, 1].set(xlabel = '', ylabel = 'mean Vj (μmol CO2 m-2 s-1)')
    axs[0, 1].xaxis.set_major_locator(MaxNLocator(5))
    axs[1, 1].plot(timePlot, resultsgco2, 'tab:brown')
    axs[1, 1].set(xlabel = 'time', ylabel = 'mean gco2 (mol CO2 m-2 s-1)')
    axs[1, 1].xaxis.set_major_locator(MaxNLocator(5))
    plt.show()
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(timePlot, results)
    axs[0, 0].set(xlabel = '', ylabel = 'E')
    axs[0, 0].xaxis.set_major_locator(MaxNLocator(5))
    axs[0, 1].plot(timePlot, resultsfw, 'tab:brown')
    axs[0, 1].set(xlabel = '', ylabel = 'fw')
    axs[0, 1].xaxis.set_major_locator(MaxNLocator(5))
    axs[1, 0].plot(timePlot, resultspl, 'tab:brown')
    axs[1, 0].set(xlabel = 'time', ylabel = 'pl')
    axs[1, 0].xaxis.set_major_locator(MaxNLocator(5))
    axs[1, 1].plot(timePlot, resultscics, 'tab:brown')
    axs[1, 1].set(xlabel = 'time', ylabel = 'ci/cs (-)')
    axs[1, 1].xaxis.set_major_locator(MaxNLocator(5))
    plt.show()
'''
fig, axs = plt.subplots(2,3, sharey=True)
axs[0, 0].plot(variables[0], resultsAn[0])
axs[0, 0].set(xlabel='Q mol quanta m-2 s-1', ylabel='mean An (μmol CO2 m-2 s-1)')
axs[0, 1].plot(variables[1], resultsAn[1], 'tab:blue')
axs[0, 1].set(xlabel='RH (-)')
axs[0, 2].plot(variables[5], resultsAn[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 2].plot(variables[5], resultsVc[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 2].plot(variables[5], resultsVj[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 2].plot(variables[5], resultsgco2[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 2].plot(variables[5], resultscics[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 1].set(xlabel='RH (-))')
axs[0, 2].plot(variables[5], resultsfw[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
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
axs[0, 2].plot(variables[5], resultspl[5], 'tab:pink')
axs[0, 2].set(xlabel='cs (mol mol-1)')
axs[1, 0].plot(variables[2], resultspl[2], 'tab:red')
axs[1, 0].set(xlabel='Tair (°C)', ylabel='mean Ψleaf (cm)')
axs[1, 1].plot(variables[3], resultspl[3], 'tab:brown')
axs[1, 1].set(xlabel='Ψ at root collar (cm)')
axs[1, 2].plot(variables[4], resultspl[4], 'tab:green')
axs[1, 2].set(xlabel='N (%)')
plt.show()
'''
