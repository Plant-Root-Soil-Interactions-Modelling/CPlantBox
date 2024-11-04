"""
produces the analytical solution Figure 4abc
from Vanderborght et al. (2005)

D. Leitner, 2018
"""
import sys; sys.path.append("../modules"); sys.path.append("../../../CPlantBox");  sys.path.append("../../../CPlantBox/src")

import functional.van_genuchten as vg

import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

sand = vg.Parameters([0.045, 0.43, 0.15, 3, 1000])
loam = vg.Parameters([0.08, 0.43, 0.04, 1.6, 50])
clay = vg.Parameters([0.1, 0.4, 0.01, 1.1, 10])

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax = [ax1, ax2, ax3]

for i, soil in enumerate([sand, loam, clay]):  # make three subplots

    if soil == sand:
        theta_sur = 0.2824
    else:
        theta_sur = soil.theta_S

    theta_i = vg.water_content(-400, soil);

    K_sur = vg.hydraulic_conductivity(vg.pressure_head(theta_sur, soil), soil);
    K_i = vg.hydraulic_conductivity(-400, soil)
    psi = lambda theta: vg.pressure_head(theta, soil)
    K = lambda psi: vg.hydraulic_conductivity(psi, soil)
    Dw = lambda psi: K(psi) / (vg.specific_moisture_storage(psi, soil))

    F = lambda theta: Dw(psi(theta)) / ((K_sur - K_i) * (theta - theta_i) - (K(psi(theta)) - K_i) * (theta_sur - theta_i))

    theta_a = (theta_sur + theta_i) / 2

    if soil == clay:  # todo: same same?
        theta_ = np.linspace (theta_i + 1e-3, theta_sur - 1e-3, 300)
    else:
        theta_ = np.linspace (theta_i + 1e-3, theta_sur - 1e-3, 300)

    delta_eta = np.zeros(len(theta_),)
    for j in range(0, len(theta_)):
        ans, err = integrate.quad(F, theta_[j], theta_a)
        delta_eta[j] = ans

    delta_eta = delta_eta * (theta_sur - theta_i)

    tv = [ [0.1, 0.2, 0.3],
         [0.2, 0.5, 1.0],
         [0.1, 0.2, 0.5]]

    x_aa = [43, 41, 27.5]  # [42.14103, 35.21381052, 23.0052]; %50;  #  how to choose reference water content and its position ????
    x_a = x_aa[i]

    t_a2 = [0.1, 0.2, 0.1]
    t_a = t_a2[i]

    eta_a = x_a - (K_sur - K_i) / (theta_sur - theta_i) * t_a
    eta = delta_eta + eta_a

    # finally, plot the thing
    lineStyle = ['b-', 'b-', 'b-']
    for j in range(0, len(tv[0])):
        t = tv[i][j]
        x = eta + (K_sur - K_i) * t / (theta_sur - theta_i);
        ax[i].plot(theta_, -x, lineStyle[i])

    ax[i].set_xlabel(r'$\theta$')
    ax[i].set_xlim(0, 0.5)

ax1.set_ylabel('Depth (cm)')
ax1.set_ylim(-150, 0)
ax2.set_ylim(-200, 0)
ax3.set_ylim(-120, 0)

if __name__ == "__main__":
    plt.show()
