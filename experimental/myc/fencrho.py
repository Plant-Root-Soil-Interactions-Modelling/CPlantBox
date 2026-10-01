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
from sklearn.linear_model import LinearRegression

# I changed the c++ code, so set might not be valid anymore
good_seeds = [10, 15, 25, 30, 45, 50, 70, 105, 115, 125]

path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

## Setting up petri dish
diameter = 9.4
radius = diameter /2
height = 1.6
    # introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height 
opening_length = 5.0
opening_height = 0.2

nRings = 15
start_dishes = time.perf_counter()
petri_dish, small_hyphae_dish, half_dish, rings = amf.makedishes(diameter, height, barrier_thickness, barrier_height, opening_length, opening_height, 0.1, nRings)

finished_dishes = time.perf_counter()
dt, nSteps = amf.setUpSimulationTime(5,24)
seed = good_seeds[0]
start_sim = time.perf_counter()
mycp, crossed_time, times = amf.makeSimulation(seed, path, name, height,petri_dish, small_hyphae_dish, half_dish,dt,nSteps, 60, "f_enc_vs_rho_h",animation=False,verbose=False)
finished_sim = time.perf_counter()
vp.plot_roots_and_container(mycp,petri_dish)
ana = pb.SegmentAnalyser()
ana = amf.getMycSegmentAnalyser(mycp)
ana.crop(small_hyphae_dish)
received_ana_data = time.perf_counter()
# fewer_rings = rings[0:int(nRings/2)]
## Ringe -> Zeile // Zeitabschnitt -> Spalte
tipMat = amf.getParaSumperRing("nodeTips",times,ana,rings[0:1])
tipMat_time = time.perf_counter()
anaMat = amf.getParaSumperRing("anastomosis",times,ana,rings[0:1])
anaMat_time = time.perf_counter()
lenMat = amf.getParaSumperRing("length",times,ana, rings[0:1])
lenMat_time = time.perf_counter()

print("Time to make the dishes ", finished_dishes - start_dishes)
print("Time to run the simulation ", finished_sim - start_sim)
print("Time to get all Paras per Ring", lenMat_time -received_ana_data)


area =np.pi /2 * radius*radius/nRings
vol =  area * height

rho = lenMat / vol 
f_enc = np.divide(anaMat, tipMat, out=np.zeros_like(anaMat, dtype=float), where=(np.array(anaMat) != 0) & (np.array(tipMat) != 0)) #
rho = rho.reshape(-1,1)
f_enc = f_enc.reshape(-1)
# f_encdivrho = np.divide(f_enc, rho, out=np.zeros_like(f_enc, dtype=float), where=(np.array(f_enc) != 0) & (np.array(rho) != 0))
model = LinearRegression(fit_intercept=False).fit(rho,f_enc)
C = model.coef_[0]
R2 = model.score(rho,f_enc)
fig, ax = plt.subplots(figsize=(7, 5))
ax.scatter(rho,f_enc,alpha=0.35,s=20)
rho_fit = np.linspace(rho.min(),rho.max(),200)
ax.plot(rho_fit,C * rho_fit,label=(rf"$C={C:.3e}$, " rf"$R^2={R2:.3f}$"))
ax.set_xlabel(r"realized hyphal length density $\rho_h$")
ax.set_ylabel(r"$f_{\mathrm{enc}}$")
ax.set_title(rf"Encounter fraction vs. $\rho_h$ ")
ax.legend()
ax.grid(alpha=0.2)
plt.tight_layout()
plt.show()