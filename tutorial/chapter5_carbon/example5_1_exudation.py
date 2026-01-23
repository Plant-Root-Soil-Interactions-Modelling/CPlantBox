"""analytic computation of citrate exudation around a growing root system"""

import time  # |\label{l5_2_exudation:libend}|

import matplotlib as mpl  # |\label{l5_2_exudation:libstart}|
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
from pyevtk.hl import gridToVTK

import plantbox as pb

font = {"size": 20}  # TODO move to figure_style.py
plt.rc("font", **font)
mpl.rcParams["mathtext.default"] = "regular"

# Root system
rs = pb.RootSystem()  # |\label{l5_2_exudation:rsstart}|
path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI"
rs.readParameters(path + name + ".xml")
rs.setSeed(1)   # |\label{l5_2_exudation:setRandomSeed}|
rs.initialize()
simtime = 10
rs.simulate(simtime, True)
rs.write("results/example5_2_exudation_rootsystem.vtp")  # |\label{l5_2_exudation:rsend}|

# Grid parameter
nodes = np.array([np.array(n) for n in rs.getNodes()])  # |\label{l5_2_exudation:gridstart}|
boxmin = nodes.min(axis=0)
boxmax = nodes.max(axis=0)
width = abs(max(boxmax[0], boxmax[1]) - min(boxmin[0], boxmin[1])) + 6  # cm
depth = abs(boxmin[2]) + 3
xres = 0.3
yres = 0.3
zres = 0.3
nx = int(width / xres)
ny = int(width / yres)
nz = int(depth / zres)
print(f"Width {width}, Depth {depth} at a Resolution {nx}*{ny}*{nz}")  # |\label{l5_2_exudation:gridend}|

# Model parameter
model = pb.ExudationModel(width, width, depth, nx, ny, nz, rs)  # |\label{l5_2_exudation:paramstart}|
model.Q = 18.4  # Citrate exudation rate (mu g/d/cm root)
model.Dl = 0.171  # Citrate liquid diffusion coefficient (cm2/d)
model.theta = 0.3  # Soil water content (cm^3/cm^3)
model.R = 16.7  # Retardation factor, (-)
model.k = 1.42  # Citrate decomposition rate  (d^-1)
model.l = 5  # Citrate depositon length behind the root tip (cm)   # |\label{l5_2_exudation:paramend}|

# Numerical parameter
model.type = pb.IntegrationType.mls  # mps, mps_straight, mls # |\label{l5_2_exudation:numparamstart}|
model.n0 = 10  # integration points per cm
model.thresh13 = 1.0e-15  # threshold to neglect diffusing g (eqn 13)
model.calc13 = True  # turns Eqn 13  on (True) and off (False)
model.observationRadius = 0.8  # limits computational domain around roots [cm]  # |\label{l5_2_exudation:numparamend}|

t = time.time()  # |\label{l5_2_exudation:runstart}|
C = model.calculate(simtime)
elapsed = time.time() - t
print(f"Computation took {elapsed} s")  # |\label{l5_2_exudation:runend}|

# Post processing...
C = np.reshape(C, (nx, ny, nz))  # |\label{l5_2_exudation:reshapestart}|
X = np.linspace(-width / 2, width / 2, nx)
Y = np.linspace(-width / 2, width / 2, ny)
Z = np.linspace(-depth, 0, nz)
X_, Y_, Z_ = np.meshgrid(X, Y, Z, indexing="ij")  # |\label{l5_2_exudation:reshapeend}|

gridToVTK("results/./example5_2_exudation_citrate", X, Y, Z, pointData={"Citrate concentration": C})  # |\label{l5_2_exudation:save}|

fig1 = plt.figure()  # |\label{l5_2_exudation:plotstart}|
ax = plt.axes()
C_ = C[:, int(ny / 2), :]
C_ = C_ + 10**-5
levels = np.linspace(np.min(np.log10(C_[:])), np.max(np.log10(C_[:])))
cs = ax.contourf(X_[:, int(ny / 2), :], Z_[:, int(ny / 2), :], np.log10(C_), levels=levels, cmap="jet")
plt.axis("equal")
cbar = fig1.colorbar(cs, pad=-0.4)
cbar.set_label(r"log10 Citrate concentration $(\mu g \, cm^{-3})$")
cbar.locator = MaxNLocator(nbins=5)
ax.set_xticks([])
ax.set_yticks([])
for pos in ["right", "top", "bottom", "left"]:
    plt.gca().spines[pos].set_visible(False)
plt.show()  # |\label{l5_2_exudation:plotend}|
