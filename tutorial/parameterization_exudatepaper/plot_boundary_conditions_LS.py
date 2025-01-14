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

font = {'size'   : 24}
plt.rc('font', **font)


"""scenario"""
year = 2019
soil_type_ = ["loam", "sand"]
soil_type1 = ["Loam", "Sand"]
genotype_ = ["WT"]
linst = ['-','--']
hatch = [' ', '--']

fig, axs = plt.subplots(1)

for i in range(0, len(soil_type_)):
    soil_type = soil_type_[i]
    for j in range(0, len(genotype_)):
        genotype = genotype_[j]
    
        soil_, min_b, max_b, cell_number, area, Kc = scenario.maize_SPP(soil_type)
        sim_time = 154   #  [day]

        #x_, y_, lai = evap.net_infiltration(year, soil_type, genotype, sim_time, Kc)
        #sim_time += 1

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

        #sim_time += 1  # root system starts to grow one day earlier (+max dat)

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
        trans = lambda t, dt:-tpot[int((t + dt / 2))] * area * sinusoidal2(t, dt)


        """ 2. load infiltration (precipitation + irrigation) """
        df = pd.read_csv("../data_magda/Inf.csv")  # precipitation data
        precip  = df["Inf_"+str(year)].loc[0: sim_time].values  #cm/d 
        time = np.linspace(0, sim_time, precip.shape[0])  # relative time in hours


        t2_ = np.linspace(0, sim_time + 0.4, len(etc) * 24)  # more points more fun...
        tpot2 = [-trans(t, 0.) / area for t in t2_]

        if i == 0: 
            axs.plot(t_, tpot, 'b',  linestyle = linst[i], label = "Potential transpiration")
            axs.plot(t_, evap,  'r', linestyle = linst[i], label = "Potential evaporation")
            axs.plot(np.nan, np.nan, 'k-', label = 'Loam')
            axs.plot(np.nan, np.nan, 'k--', label = 'Sand')
        else: 
            axs.plot(t_, tpot, 'b',  linestyle = linst[i])
            axs.plot(t_, evap,  'r', linestyle = linst[i])
        axs.fill_between(x=t_, y1= lai(t_), color= "grey",hatch = hatch[i], alpha= 0.3, label = 'LAI, '+soil_type1[i])
        axs.set_xlabel("Time (d)")
        axs.set_ylabel(r'LAI (-) / Potential evaporation ($mm\ d^{-1}$) /''\n'' Potential transpiration ($mm\ d^{-1}$)')
        axs.set_ylim([0, 10])

        ax2 = axs.twinx()
        ax2.bar(time, precip, label = 'Precipitation and irrigation') 
        ax2.set_ylabel('Precipitation + irrigation ($mm\ d^{-1}$)')
        ax2.set_ylim([0, 100])
        ax2.invert_yaxis()

# ask matplotlib for the plotted objects and their labels
lines, labels = axs.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc='best')
plt.show()
