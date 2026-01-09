"""postprocessing... plots results from 7.2 & 7.3
a) Transpiration over time
TODO
"""

import figure_style as st
import matplotlib.pyplot as plt
import numpy as np


def sinusoidal(t):
    """sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0


filename = "Zeamays_synMRI_modified"
trans = 250  # cm3 day-1 (sinusoidal) = mL day-1
sim_time = 7.5  # day
dt = 360.0 / (24 * 3600)

n_steps = round(sim_time / dt)
ind = int((2.5 / sim_time * n_steps) // 10)

# vp.plot_roots_and_soil_files("example72_{:06d}".format(ind), "pressure head", "subType")  # "water content", "pressure head"

data72 = np.load("results/" + filename + ".npy")
data73 = np.load("results/" + filename + "_fp.npy")

# Transpiration over time
x1_ = data72[0,:]
y1_ = data72[1,:]
fig, ax = st.subplots21()
ax1 = ax[0]
ax1.plot(x1_, trans * sinusoidal(x1_), "k", label = "Potential")
ax1.plot(x1_, np.array(y1_), "g", label = "Actual")
ax1.set_xlabel("Time (day)")
ax1.set_ylabel("Transpiration (mL day$^{-1}$) per plant")
ax1.legend(loc = "upper left")
ax2 = ax1.twinx()
ax2.plot(x1_, np.cumsum(np.array(y1_) * dt), "c--", label = "Cumulative")
ax2.set_ylabel("Cumulative transpiration (mL)")
ax2.legend(loc = "upper right")
plt.tight_layout()
plt.savefig("results/figure7_2_trans.png")

# Transpiration over time
x2_ = data73[0,:]
y2_ = data73[1,:]
z2_ = data73[2,:]
ax1 = ax[1]
ax1.plot(x2_, trans * sinusoidal(x2_), "k", label = "Potential")
ax1.plot(x2_, np.array(y2_), "g", label = "Actual nonlinear")
ax1.plot(x1_, np.array(y1_), "g:", label = "Actual")
ax1.set_xlabel("Time (day)")
ax1.set_ylabel("Transpiration (mL day$^{-1}$) per plant")
ax1.legend(loc = "upper left")
ax2 = ax1.twinx()
ax2.plot(x2_, np.cumsum(np.array(y2_) * dt), "c", label = "Cumulative nonlinear")
ax2.plot(x1_, np.cumsum(np.array(y1_) * dt), "c:", label = "Cumulative")
ax2.set_ylabel("Cumulative transpiration (mL)")
ax2.legend(loc = "upper right")
plt.tight_layout()
plt.savefig("results/figure7_3_trans.png")
plt.show()
