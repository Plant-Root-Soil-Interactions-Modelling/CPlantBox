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

font = {'size'   : 18}
plt.rc('font', **font)



#data Daniela
DAS = [42, 154] #d
Vmax = [2.54e-6, 1.98e-6] #g C / g soil / h 
Km = [0.12e-3, 0.06e-3]#g C / g soil
print(np.mean(Vmax), np.mean(Km))
bd = 1.4 #g/cm³
conc = np.linspace(0,.001,1000) #g C/ g soil 

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
plt.plot(conc*1e6*bd, v_mean*1e6*bd*24, label = 'mean')
plt.xlabel('concentration ($\mu$g C/ $cm^{3}$ soil)')
plt.ylabel('decay rate ($\mu$g C/ $cm^{3}$ soil / d)')
#plt.legend()
plt.tight_layout()
plt.show()

print('Vmax', np.mean(Vmax)*1e6*bd*24)
print('Km', np.mean(Km)*1e6*bd)

#test
#plt.plot(conc*1e3*bd, v_42*1e6, label = '42 days') 
#plt.plot(conc*1e3*bd, v_154*1e6, label = '154 days')
#plt.xlabel('concentration ( mg C/ g soil)')
#plt.ylabel('decay ($\mu$ g C/ g soil / h)')
#plt.show()
