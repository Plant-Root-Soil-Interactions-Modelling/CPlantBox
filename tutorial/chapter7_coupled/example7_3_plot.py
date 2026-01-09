"""plots results from 7.2_coupling & 7.3_coupling_fp into one plot"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from plantbox.visualisation import figure_style
import plantbox.visualisation.vtk_plot as vp


def sinusoidal(t):
    """sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0


filename = "Zeamays_synMRI_modified"
t_pot = 250  # cm3 day-1 (sinusoidal) = mL day-1

data72 = np.load("results/" + filename + ".npy")
sim_times72_ = data72[0, :]
t_act72_ = data72[1, :]

data73 = np.load("results/" + filename + "_fp.npy")
sim_times73_ = data73[0, :]
t_act73_ = data73[1, :]
# q_soil73_ = data73[2, :]

# Transpiration over time
fig, ax = figure_style.subplots12(1, 1)
ax.plot(sim_times73_, t_pot * sinusoidal(sim_times73_), "k", label="Potential")
ax.plot(sim_times73_, np.array(t_act73_), "g", label="Actual nonlinear")
# ax.plot(sim_times73_, np.array(-q_soil73_), "b:", label="Actual nonlinear")
ax.plot(sim_times72_, np.array(t_act72_), "g:", label="Actual")
ax.set_xlabel("Time (day)")
ax.set_ylabel("Transpiration (mL day$^{-1}$) per plant")
ax.legend(loc="upper left")
ax2 = ax.twinx()
ax2.plot(sim_times73_, cumulative_trapezoid(t_act73_, sim_times73_, initial=0), "c", label="Cumulative nonlinear")
ax2.plot(sim_times72_, cumulative_trapezoid(t_act72_, sim_times72_, initial=0), "c:", label="Cumulative")
ax2.set_ylabel("Cumulative transpiration (mL)")
ax2.legend(loc="upper right")
plt.tight_layout()
plt.savefig("results/figure7_3_trans.png")
plt.show()

# sim_time = 7.5  # day
# dt = 360.0 / (24 * 3600)
# n_steps = round(sim_time / dt)
# ind = int((2.5 / sim_time * n_steps) // 10) # only every 10th time step is written
# vp.plot_roots_and_soil_files("example72_{:06d}".format(ind), "pressure head", "subType")  # "water content", "pressure head"
# vp.plot_roots_and_soil_files("example73_{:06d}".format(ind), "pressure head", "subType")  # "water content", "pressure head"
