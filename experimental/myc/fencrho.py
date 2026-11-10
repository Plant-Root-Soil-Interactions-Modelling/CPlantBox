import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import time
import matplotlib as mpl
import AMFAnalysis as amf
import AMFDensities as amfd
import AMFPlots as amfp
import AMFSimulationSetup as amfs
from sklearn.linear_model import LinearRegression

# I changed the c++ code, so set might not be valid anymore
# good_seeds = [10, 15, 25, 30, 45, 50, 70, 105, 115, 125]
# good_seeds = [25+ i for i in range(10)]
good_seeds = [25,10]
path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

## Setting up petri dish (geometry of the experiment, height from the filled medium volumes, see amfs.defaultPetridish)
hyphae_layer = 0.1 # hyphae grow on the cellophane in the fungal compartment, i.e. in a thin layer at the top
nRings = 15
start_dishes = time.perf_counter()
dish = amfs.defaultPetridish(nRings=nRings, hyphae_layer=hyphae_layer)
radius, height, barrier_thickness = dish["radius"], dish["height"], dish["barrier_thickness"]
petri_dish, small_hyphae_dish, half_dish, rings = dish["petri_dish"], dish["small_hyphae_dish"], dish["half_dish"], dish["rings"]

finished_dishes = time.perf_counter()
dt, nSteps = amfs.setUpSimulationTime(5,24)
# seed = good_seeds[0]
simulation = []
for seed in good_seeds:
    start_sim = time.perf_counter()
    mycp, time_star, times = amfs.makeSimulation(seed, path, name, height,petri_dish, small_hyphae_dish, half_dish,dt,nSteps, 60, "compare_COCO",animation=False,verbose=False)
    # mycp, times = amf.makeFreeSimulation(seed, path, name, dt, nSteps,"colonizationvisforfspm",True,False)
    finished_sim = time.perf_counter()
    # vp.plot_roots_and_container(mycp,petri_dish)
    # ana = amf.getMycSegmentAnalyser(mycp,filename="colonizationvisforfspm",stdwrite=True)
    # ana.crop(small_hyphae_dish)
    received_ana_data = time.perf_counter()
    # fewer_rings = rings[0:int(nRings/2)]
    # Ringe -> Zeile // Zeitabschnitt -> Spalte
    mats = amfd.getParasPerRing(["nodeTips","anastomosis","length"],times,mycp,rings)
    tipMat, anaMat, lenMat = mats["nodeTips"], mats["anastomosis"], mats["length"]
    lenMat_time = time.perf_counter()
    vol = amfs.ringVolume(radius, nRings, hyphae_layer, barrier_thickness) # one volume per ring, rings are in the layer
    simulation.append({"label": "seed " + str(seed), "mycp": mycp, "time_star": time_star, "times": times, "tipMat": tipMat, "anaMat": anaMat, "lenMat": lenMat, "vol": vol})

    print("Time to make the dishes ", finished_dishes - start_dishes)
    print("Time to run the simulation ", finished_sim - start_sim)
    print("Time to get all Paras per Ring", lenMat_time -received_ana_data)

all_results = amf.analyze_simulations(simulation, window = 0.05, step=0.01,min_points =20)

amfp.plot_data(all_results)
amfp.plot_linear_analysis(all_results)


# rho = lenMat / vol 
# f_enc = np.divide(anaMat, tipMat, out=np.zeros_like(anaMat, dtype=float), where=(np.array(anaMat) != 0) & (np.array(tipMat) != 0)) #
# rho = rho.reshape(-1,1)
# f_enc = f_enc.reshape(-1)
# # f_encdivrho = np.divide(f_enc, rho, out=np.zeros_like(f_enc, dtype=float), where=(np.array(f_enc) != 0) & (np.array(rho) != 0))
# model = LinearRegression(fit_intercept=False).fit(rho,f_enc)
# C = model.coef_[0]
# R2 = model.score(rho,f_enc)
# fig, ax = plt.subplots(figsize=(7, 5))
# ax.scatter(rho,f_enc,alpha=0.35,s=20)
# rho_fit = np.linspace(rho.min(),rho.max(),200)
# ax.plot(rho_fit,C * rho_fit,label=(rf"$C={C:.3e}$, " rf"$R^2={R2:.3f}$"))
# ax.set_xlabel(r"realized hyphal length density $\rho_h$")
# ax.set_ylabel(r"$f_{\mathrm{enc}}$")
# ax.set_title(rf"Encounter fraction vs. $\rho_h$ ")
# ax.legend()
# ax.grid(alpha=0.2)
# plt.tight_layout()
# plt.show()