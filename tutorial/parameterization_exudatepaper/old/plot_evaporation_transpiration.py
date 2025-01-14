""" 
    Maize using rhizosphere models  
"""
import sys;
sys.path.append("../");
sys.path.append("../../../build-cmake/cpp/python_binding/");
sys.path.append("../data_magda/");
sys.path.append("../../modules/");
sys.path.append("../../../../CPlantBox/");
sys.path.append("../../../../CPlantBox/src")
sys.path.append("../../../../CPlantBox/src/functional/");
sys.path.append("../../../../CPlantBox/src/rsml/");
sys.path.append("../../../../CPlantBox/src/visualisation/")
sys.path.append("../../../../CPlantBox/src/structural/")
sys.path.append("../../../../CPlantBox/src/external/")
import plantbox as pb  # CPlantBox
import visualisation.vtk_plot as vp
import numpy as np
import scenario_setup as scenario
#import evapotranspiration as evap
import matplotlib.pyplot as plt
import pandas as pd
from xylem_flux import sinusoidal2
from scipy import interpolate
import os

#colors
a = [0, 0, 1]
b = [106/255, 166/255, 1]
c = [1, 0, 0]
d = [1, 128/255, 0]
cols = [a,b,c,d]
label_root = ["L_WT", "L_rth3", "S_WT", "S_rth3"]
soil_type_ = ["loam", "loam", "sand", "sand"]
genotype_ = ["WT", "RTH3","WT", "RTH3"]
year = 2019

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title


fig, (ax1, ax2) = plt.subplots(1, 2)
for i in range(0, len(soil_type_)):
    soil_type = soil_type_[i]
    genotype = genotype_[i]

    soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
    sim_time = 154   #  [day]

    df = pd.read_csv("../data_magda/Inf.csv")  # precipitation data
    yd1  = df["Inf_"+str(year)].loc[0: sim_time].values*0.1 #nmm/d - cm/d
    t_ = np.linspace(0, sim_time, yd1.shape[0] * 24)  # relative time in hours
    precip = np.array([ yd1[int(t)] * sinusoidal2(t, 0.) for t in t_ ])

    """ 1. load ET, ET0 -> ETc """
    df2 = pd.read_csv("../data_magda/ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = np.array([ yd2[int(t)] * sinusoidal2(t, 0.) for t in t_ ])

    """2. interpolate Kc """
    f_Kc = interpolate.interp1d(np.linspace(0,sim_time, sim_time+1), Kc[0:sim_time+1])
    Kc_hours = f_Kc(t_)
    etc = et0 * Kc_hours

    """3. load LAI """
    if soil_type == 'loam':
        st = 'L'
    else:
        st = 'S'

    df3 = pd.read_csv("../data_magda/LAI_"+str(year)+".csv")  # LAI
    LAI =  df3[st+"_"+genotype].loc[0:sim_time].values
    lai = interpolate.interp1d(np.linspace(0,sim_time, sim_time+1), LAI)

    """ 4. ETc -> Tpot, Evap """
    k = 0.45
    tpot = np.multiply(etc, [1. - np.exp(-k *lai(t)) for t in t_])
    evap = etc - tpot
    evap = -evap
    net_inf = precip + evap

    y_ = []
    x_ = []
    for p in range(0, t_.shape[0] - 1):
        x_.extend([t_[p], t_[p + 1]])  # hour -> day
        y_.extend([net_inf[p], net_inf[p]])
    x_ = np.array(x_)
    y_ = np.array(y_) 


    """ 0. load ET, ET0 -> ETc """
    df2 = pd.read_csv("../data_magda/ET0.csv")  # evapotranspiration  in cm / day 
    yd2 = df2["ET0_"+str(year)].loc[0:sim_time].values # cm / d
    et0 = yd2
    etc = et0 * Kc[0:sim_time+1]

    """ 1. ETc -> Tpot, Evap """
    k = 0.45
    t_ = np.linspace(0, len(etc) - 1, len(etc))
    tpot = np.multiply(etc*10, [(1. - np.exp(-k * lai(t_[p]))) for p in range(0, len(etc))])
    evap = (etc*10 - tpot)

    ax1.plot(t_, tpot, color = cols[i],  label = label_root[i])
    ax2.plot(t_, evap, color = cols[i])

ax1.set_xlabel("Time (d)")
ax2.set_xlabel("Time (d)")
ax2.set_ylabel(r'Potential evaporation ($mm\ d^{-1}$)')
ax1.set_ylabel(r'Potential transpiration ($mm\ d^{-1}$)')
ax1.set_ylim([0, 7])
ax2.set_ylim([0, 7])
ax1.legend()

# ask matplotlib for the plotted objects and their labels
#lines, labels = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax2.legend(lines + lines2, labels + labels2, loc='best')
plt.show()
