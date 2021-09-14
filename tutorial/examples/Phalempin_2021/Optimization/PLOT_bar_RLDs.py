import sys; 
sys.path.append("../../../..")
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors
import math
import plantbox as rb
from scipy.stats import pearsonr
#############################################################

phase = ['p1', 'p2']
param = ['RLD']
corenum = 3
pearsonR = np.zeros((2,))
rRMSE = np.zeros((2,))

for ii in range(0,len(phase)):
    for jj in range(0,len(param)):    

        with open("Pearson/"+phase[ii]+"_"+param[jj]+"_meas.txt") as f: 
            meas = [[float(xx) for xx in line.split()] for line in f]
        with open("Pearson/"+phase[ii]+"_"+param[jj]+"_virt.txt") as f: 
            virt = [[float(xx) for xx in line.split()] for line in f]

            #print(meas)
            #print(virt) 
            meas = np.reshape(np.array(meas), (int(len(meas)),corenum))
            virt = np.reshape(np.array(virt), (int(len(virt)),corenum))
            print(np.nanmean(meas,axis= 0))
            print(np.nanmean(virt,axis= 0))

            # calculate Pearson's correlation
            x = np.nanmean(meas,axis=0)
            y = np.nanmean(virt,axis=0)
            print(x, y) 
            corr_coeff, _ = pearsonr(x,y)
            print(phase[ii], param[jj])
            print('Pearsons correlation: %.3f' % corr_coeff**2)
            pearsonR[ii] = round(corr_coeff,3)

            # calculate rRMSE
            dummy = (y-x)**2
            n = len(x)
            #print((np.sum(dummy/n))**0.5)
            rRMSE[ii] = (np.sum(dummy/n))**0.5*100/np.mean(x)
            print('rRMSE: %.3f' % rRMSE[ii])
            rRMSE[ii] = round(rRMSE[ii],0)

###############################################
vars = ["RLDs"]
vars1 = ["RLD"]
mult = [1]
unit = [' ($cm / cm^3$)']

plt.rcParams.update({'font.size': 16})
plt.rcParams['figure.constrained_layout.use'] = True
fig, axs = plt.subplots(2,3)
axs = axs.ravel()

for j in range(0,len(vars)): 

    A = np.zeros((6,2))
    B = np.zeros((6,2))

    with open(vars[j]+"/p1_measured_"+vars1[j]+"_std_loam_low.txt") as f: 
        meas1 = [[float(xx) for xx in line.split()] for line in f]
    with open(vars[j]+"/p2_measured_"+vars1[j]+"_std_loam_low.txt") as f: 
        meas2 = [[float(xx) for xx in line.split()] for line in f]
    A[:,0] = np.reshape(np.array(meas1), (6,))
    A[:,1] = np.reshape(np.array(meas2), (6,))

    with open(vars[j]+"/p1_virt_rm_std_loam_low.txt") as f: 
        virt1 = [[float(xx) for xx in line.split()] for line in f]
    with open(vars[j]+"/p2_virt_rm_std_loam_low.txt") as f: 
        virt2 = [[float(xx) for xx in line.split()] for line in f]
    B[:,0] = np.reshape(np.array(virt1), (6,))
    B[:,1] = np.reshape(np.array(virt2), (6,))



    X = np.arange(3)
    ym_p1 = A[0:3,0]
    yv_p1 = B[0:3,0]
    ym_p2 = A[0:3,1]
    yv_p2 = B[0:3,1]

    y_p1_m = A[3:6,0]
    y_p1_v = B[3:6,0]
    y_p2_m = A[3:6,1]
    y_p2_v = B[3:6,1]


    axs[j].barh(X - 0.125, ym_p1*mult[j], xerr = y_p1_m*mult[j], capsize=4, color = 'b', height = 0.25)
    axs[j].barh(X + 0.125, yv_p1*mult[j], xerr = y_p1_v*mult[j], capsize=4, color = 'g', height = 0.25)

    axs[j+1].barh(X - 0.125, ym_p2*mult[j], xerr = y_p2_m*mult[j], capsize=4, color = 'b', height = 0.25)
    axs[j+1].barh(X + 0.125, yv_p2*mult[j], xerr = y_p2_v*mult[j], capsize=4, color = 'g', height = 0.25)
    
    axs[j+1].barh(np.nan, np.nan, color = 'b', label = 'measurement')
    axs[j+1].barh(np.nan, np.nan, color = 'g', label = 'simulation')

    axs[j].set_yticks(range(0,3)) #,['core 3','core 2','core 1'])
    axs[j+1].set_yticks(range(0,3))
    axs[j].set_xlabel(vars1[j]+unit[j])
    axs[j+1].set_xlabel(vars1[j]+unit[j])

    axs[j].invert_yaxis()
    axs[j+1].invert_yaxis()

    labels = ['core 1','core 2','core 3']
    labels1 = [' ', ' ', ' ']
    axs[j].set_yticklabels(labels)
    axs[j+1].set_yticklabels(labels1)

    axs[j].set_xlim([0,50])
    axs[j+1].set_xlim([0,5])
    
    axs[j].set_title('phase 1')
    axs[j+1].set_title('phase 2')
    

    axs[j].text(33,0,"R = "+ str(pearsonR[0])+"\nrRMSE = "+str(int(rRMSE[0]))+"%")
    axs[j+1].text(3.2,0,"R = "+ str(pearsonR[1])+"\nrRMSE = "+str(int(rRMSE[1]))+"%")


    if j==1:
        axs[j+1].legend()

axs[2].axis('off')
axs[3].axis('off')
axs[4].axis('off')
axs[5].axis('off')
#fig.tight_layout()
plt.show()



