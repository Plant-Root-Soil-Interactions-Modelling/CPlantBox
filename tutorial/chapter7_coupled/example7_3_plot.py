""" postprocessing... plots results from 7.2 & 7.3
a) Transpiration over time

"""
import sys; sys.path.append(".."); sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import figure_style as st

import matplotlib.pyplot as plt
import numpy as np


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


name = "Zeamays_synMRI_modified"
trans = 250  # cm3 /day (sinusoidal) = mL/day
sim_time = 7.5  # [day]
dt = 360. / (24 * 3600)

N = round(sim_time / dt)
ind = int((2.5 / sim_time * N) // 10)

# vp.plot_roots_and_soil_files("example72_{:06d}".format(ind), "pressure head", "subType")  # "water content", "pressure head"

data72 = np.load("results/" + name + ".npy")
data73 = np.load("results/" + name + "_fp.npy")

""" Transpiration over time """
x1_ = data72[0,:]
y1_ = data72[1,:]
fig, ax1 = st.subplots21()
ax1.plot(x1_, trans * sinusoidal(x1_), 'k', label = "Potential")
ax1.plot(x1_, np.array(y1_), 'g', label = "Actual")
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend(loc = 'upper left')
ax2 = ax1.twinx()
ax2.plot(x1_, np.cumsum(np.array(y1_) * dt), 'c--', label = "Cumulative")
ax2.set_ylabel("Cumulative transpirtation $[mL]$")
ax2.legend(loc = 'upper right')
plt.tight_layout()
plt.savefig("results/figure7_2_trans.png")
plt.show()

""" Transpiration over time """
x2_ = data73[0,:]
y2_ = data73[1,:]
z2_ = data73[2,:]
fig, ax1 = st.subplots21()
ax1.plot(x2_, trans * sinusoidal(x2_), 'k', label = "Potential")
ax1.plot(x2_, np.array(y2_), 'g', label = "Actual nonlinear")
ax1.plot(x1_, np.array(y1_), 'g:', label = "Actual")
ax1.set_xlabel("Time $[d]$")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend(loc = 'upper left')
ax2 = ax1.twinx()
ax2.plot(x2_, np.cumsum(np.array(y2_) * dt), 'c', label = "Cumulative nonlinear")
ax2.plot(x1_, np.cumsum(np.array(y1_) * dt), 'c:', label = "Cumulative")
ax2.set_ylabel("Cumulative transpirtation $[mL]$")
ax2.legend(loc = 'upper right')
plt.tight_layout()
plt.savefig("results/figure7_3_trans.png")
plt.show()

