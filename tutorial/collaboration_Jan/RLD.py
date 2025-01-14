"""root system surface density"""
import sys; sys.path.append("../..");
sys.path.append("../../src/")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

font = {'size'   : 16}
plt.rc('font', **font)

#simulation 
path = "rootsystem/"
names = ["RS_optimized_field_L_WT_new",  "RS_optimized_field_S_WT_new"]
title = ["Loam", "Sand"]
cols = ['r','b','g','c']
depth = 80
layers = 9
times = [42, 63, 98, 154]
BBCHs = [14,19,59,83]

fig, ax = plt.subplots(1, 2)
xx = []
for i in range(0,len(names)):
    rld_ = []
    dummy = 0
    name = names[i]
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    rs.initializeLB(5,4)
    rl_ = []
    for j in range(0,times[-1]+1):
        if i ==0:
            if j<=98:
                rs.simulate(1, True)
        elif i == 1:
            rs.simulate(1, True)

        if j in times:
            ana = pb.SegmentAnalyser(rs)
            rl = ana.distribution("length", 0., -depth, layers, True)
            soilvolume = (depth / layers) * 45 * 20
            RLD = np.array(rl) / soilvolume  # convert to density
            rld_.append(RLD)

            z_ = np.linspace(0, -depth, layers)  # z - axis
            ax[i].plot(RLD, z_, color = cols[dummy], linestyle = "--")
            ax[i].title.set_text(title[i])
            ax[i].set_xlabel("RLD $(cm$ $ cm^{-3})$")
            ax[i].set_ylabel("Depth $(cm)$")
            dummy = dummy+1

    #write computed RLD to csv file
    xx.extend(np.transpose(np.vstack((np.repeat(i,len(z_)),z_,rld_))))
np.savetxt('RLD.csv', xx, delimiter = ",", header="treatment,depth,BBCH14,BBCH19,BBCH59,BBCH83",fmt="%0.4f")


            
#Helena's data
df = pd.read_csv("data/RLD_2019_Helena.csv")
depth = df['depth'].values
BBCH = df["BBCH"].values
treatment = df["treatment"].values
RLD = df["RLD"].values

treatment_ = np.unique(treatment)
BBCH_ = np.unique(BBCH)
depth_ = np.unique(depth)

data = np.zeros((len(treatment_), len(BBCH_),len(depth_)))
data_SE = np.zeros((len(treatment_), len(BBCH_),len(depth_)))
for i in range(0, len(treatment_)):
    for j in range(0, len(BBCH_)):
        for k in range(0, len(depth_)):
            data_ = df[(df['treatment']==treatment_[i])&(df['BBCH']==BBCH_[j])&(df['depth']==depth_[k])]
            if data_["RLD"].loc[:].values.size > 0:
                data[i,j,k] = data_["RLD"].loc[:].values
                data_SE[i,j,k] = data_["RLD_SE"].loc[:].values

#PLOT
dummy = [0,2]
for i in range(0,2):
    for j in range(0, len(BBCH_)):
        x = data[dummy[i],j,:] 
        error = data_SE[dummy[i],j,:]
        ax[i].plot(x,depth_*-1+10, color = cols[j], label = "BBCH "+str(BBCHs[j]))
        ax[i].fill_betweenx(depth_*-1+10, x-error, x+error, color = cols[j],alpha = 0.3)
        ax[i].set_xlim([-0.1,4])
        ax[i].set_ylim([-75,0])


plt.plot(np.nan,np.nan, color = 'k', linestyle = '-', label = 'Helena')
plt.plot(np.nan,np.nan, color = 'k', linestyle = '--', label = 'Simulation')
        
plt.legend()
plt.show()
