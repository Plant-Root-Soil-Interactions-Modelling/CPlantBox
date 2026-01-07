"""
produces the analytical solution Figure 4abc
from Vanderborght et al. (2005)

D. Leitner, 2018
"""
<<<<<<< HEAD
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

from functional.van_genuchten import *

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from math import *

sand = Parameters([0.045, 0.43, 0.15, 3, 1.1574e-04 * 100 * 3600 * 24])
loam = Parameters([0.08, 0.43, 0.04, 1.6, 5.7870e-06 * 100 * 3600 * 24])
clay = Parameters([0.1, 0.4, 0.01, 1.1, 1.1574e-06 * 100 * 3600 * 24])
=======

import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

import plantbox.functional.van_genuchten as vg

sand = vg.Parameters([0.045, 0.43, 0.15, 3, 1.1574e-04 * 100 * 3600 * 24])
loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 5.7870e-06 * 100 * 3600 * 24])
clay = vg.Parameters([0.1, 0.4, 0.01, 1.1, 1.1574e-06 * 100 * 3600 * 24])
>>>>>>> origin/master

jwpot_ = [-0.1, -0.1, -0.3, -0.3]
head_i_ = [-40, -200, -200, -200]

N = 1000
y = np.zeros((N, 4))
t = np.linspace(0, 10, N)  # days

for i, soil in enumerate([sand, loam, loam, clay]):
    head_i = head_i_[i]
<<<<<<< HEAD
    theta_i = water_content(head_i, soil)  # initial theta
    theta_sur = water_content(-10000, soil)  # critical vaule
    jwpot = jwpot_[i]

    # TH = (theta-theta_sur)/(theta_i-theta_sur)
    dw = lambda TH: water_diffusivity(TH, theta_i, theta_sur, soil)
    int_dw, err = integrate.quad(dw, 0, 1)

    theta_dw = lambda TH: TH * water_diffusivity(TH, theta_i, theta_sur, soil)
    int_theta_dw, err = integrate.quad(theta_dw, 0, 1)
    beta = pow(int_theta_dw / int_dw, 2)  # 43

    fun_dw = lambda TH: pow(1 - TH * beta, 2) * dw(TH)
=======
    theta_i = vg.water_content(head_i, soil)  # initial theta
    theta_sur = vg.water_content(-10000, soil)  # critical vaule
    jwpot = jwpot_[i]

    # TH = (theta-theta_sur)/(theta_i-theta_sur)
    def dw(TH):
        return vg.water_diffusivity(TH, theta_i, theta_sur, soil)

    int_dw, err = integrate.quad(dw, 0, 1)

    def theta_dw(TH):
        return TH * vg.water_diffusivity(TH, theta_i, theta_sur, soil)

    int_theta_dw, err = integrate.quad(theta_dw, 0, 1)
    beta = pow(int_theta_dw / int_dw, 2)  # 43

    def fun_dw(TH):
        return pow(1 - TH * beta, 2) * dw(TH)

>>>>>>> origin/master
    alpha, err = integrate.quad(fun_dw, 0, 1)
    alpha /= int_dw  # 42

    mu = (3 * beta * (1 + np.sqrt(1 - (14 / 9) * (1 - alpha / pow(1 - beta, 2))))) / (2 * (1 - beta) * (alpha / pow(1 - beta, 2) - 1))  # eq 41

<<<<<<< HEAD
    sw = lambda theta_sur, theta_i: (theta_i - theta_sur) * sqrt((4 / mu) * int_dw)  # eq 39
=======
    def sw(theta_sur, theta_i):
        return (theta_i - theta_sur) * np.sqrt((4 / mu) * int_dw)  # eq 39
>>>>>>> origin/master

    tdash = (sw(theta_sur, theta_i) * sw(theta_sur, theta_i)) / (4 * jwpot * jwpot)  # eq 44
    tpot = (sw(theta_sur, theta_i) * sw(theta_sur, theta_i)) / (2 * jwpot * jwpot)  # eq 45
    print("Scenario ", i, " tpot", tpot)

<<<<<<< HEAD
    jw = lambda t: (t < tpot) * jwpot + (t >= tpot) * sw(theta_sur, theta_i) / (2 * sqrt(abs(tdash + t - tpot)))  # eq 46 & 47
=======
    def jw(t):
        return (t < tpot) * jwpot + (t >= tpot) * sw(theta_sur, theta_i) / (2 * np.sqrt(abs(tdash + t - tpot)))  # eq 46 & 47
>>>>>>> origin/master

    y[:, i] = list(map(jw, t))  # evaluate

#
# prepare plot
#
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

<<<<<<< HEAD
ax1.plot(t, abs(y[:, 0]), 'b')
ax1.set_ylabel('$E_{act}$ (cm day$^{-1}$)',fontsize=20)
ax1.set_xlim(0, 1)
ax1.set_title("Sand",fontsize=20)

ax2.plot(t, abs(y[:, 1]), 'b')
ax2.set_xlim(0, 10)
ax2.set_title("Loam",fontsize=20)

ax3.plot(t, abs(y[:, 2]), 'b')
ax3.set_xlabel('$t$ (days)',fontsize=20)
ax3.set_ylabel('$E_{act}$ (cm day$^{-1}$)',fontsize=20)
ax3.set_xlim(0, 2)
ax3.set_title("Loam",fontsize=20)

ax4.plot(t, abs(y[:, 3]), 'b')
ax4.set_xlabel('$t$ (days)',fontsize=20)
ax4.set_xlim(0, 6)
ax4.set_title("Clay",fontsize=20)

ax1.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax3.tick_params(axis='both', which='major', labelsize=16)
ax4.tick_params(axis='both', which='major', labelsize=16)
=======
ax1.plot(t, abs(y[:, 0]), "b")
ax1.set_ylabel("$E_{act}$ (cm day$^{-1}$)", fontsize=20)
ax1.set_xlim(0, 1)
ax1.set_title("Sand", fontsize=20)

ax2.plot(t, abs(y[:, 1]), "b")
ax2.set_xlim(0, 10)
ax2.set_title("Loam", fontsize=20)

ax3.plot(t, abs(y[:, 2]), "b")
ax3.set_xlabel("$t$ (days)", fontsize=20)
ax3.set_ylabel("$E_{act}$ (cm day$^{-1}$)", fontsize=20)
ax3.set_xlim(0, 2)
ax3.set_title("Loam", fontsize=20)

ax4.plot(t, abs(y[:, 3]), "b")
ax4.set_xlabel("$t$ (days)", fontsize=20)
ax4.set_xlim(0, 6)
ax4.set_title("Clay", fontsize=20)

ax1.tick_params(axis="both", which="major", labelsize=16)
ax2.tick_params(axis="both", which="major", labelsize=16)
ax3.tick_params(axis="both", which="major", labelsize=16)
ax4.tick_params(axis="both", which="major", labelsize=16)

>>>>>>> origin/master
if __name__ == "__main__":
    plt.show()
