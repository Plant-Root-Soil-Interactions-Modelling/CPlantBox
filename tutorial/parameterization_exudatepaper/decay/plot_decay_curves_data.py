import sys
sys.path.append("../../../../../CPlantBox")
sys.path.append("../../../../../CPlantBox/src")
import plantbox as pb
import visualisation.vtk_plot as vp

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import pandas as pd
import numpy as np
import matplotlib as mpl

plt.rcParams.update({'font.size': 16})
mpl.rcParams['mathtext.default'] = 'regular'

#data Daniela raw
df1 = pd.read_csv("data/early_exudates.csv")
df2 = pd.read_csv("data/late_exudates.csv")

#data Daniela
DAS = [42, 154] #d
Vmax = [2.54e-6, 1.98e-6] #g C / g soil / h 
Km = [0.12e-3, 0.06e-3]#g C / g soil
bd = 1.4 #g/cm³
print(np.mean(Vmax), np.mean(Km))
print(np.mean(Vmax)*1e6*bd*24, np.mean(Km)*1e6*bd)
conc = np.linspace(0,.002,1000) #g C/ g soil 

v_42= np.zeros((len(conc)))
v_154 = np.zeros((len(conc)))
v_mean = np.zeros((len(conc)))
for i in range(0, len(conc)):
    v_42[i] = (Vmax[0]*conc[i])/(Km[0]+conc[i]) #g C / g soil / h
    v_154[i] = (Vmax[1]*conc[i])/(Km[1]+conc[i])
    v_mean[i] = (np.mean(Vmax)*conc[i])/(np.mean(Km)+conc[i])

#conc: [microg C / cm3 soil ]
#v: [microg C / cm3 soil / d]
#plt.plot(conc*1e6*bd, v_42*1e6*bd*24, label = 'early: 42 DAS') 
#plt.plot(conc*1e6*bd, v_154*1e6*bd*24, label = 'late: 154 DAS')
plt.scatter(df1["concentration"]*1e3*bd, df1["decay rate"]*bd*24, label = 'DAS42')
plt.scatter(df2["concentration"]*1e3*bd, df2["decay rate"]*bd*24, label = 'DAS154')
plt.plot(conc*1e6*bd, v_mean*1e6*bd*24,  color = 'k')
plt.text(400, 30, "$V_{max}$ =76 $\mu$gC $cm^{-3}$ $d^{-1}$")
plt.text(400, 20, "$K_{m}$ =126 $\mu$gC $cm^{-3}$")
plt.xlabel("Concentration ($\mu$g C $cm^{-3}$ soil)")
plt.ylabel("Decay rate ($\mu$g C $cm^{-3}$ soil  $d^{-1}$)")
plt.legend(loc = 'upper right')
plt.show()

#test
#plt.plot(conc*1e3*bd, v_42*1e6, label = '42 days') 
#plt.plot(conc*1e3*bd, v_154*1e6, label = '154 days')
#plt.xlabel('concentration ( mg C/ g soil)')
#plt.ylabel('decay ($\mu$ g C/ g soil / h)')
#plt.show()
