import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pandas as pd
import sys
import scipy

font = {'size'   : 20}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'

left  = 0.08  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.4   # the amount of width reserved for blank space between subplots
hspace = 0.4   # the amount of height reserved for white space between subplots

def rsquared(x, y):
    """ Return R^2 where x and y are array-like."""
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2
    
def func(x, a):
    return a * x

df1 = pd.read_csv("../results/Zea_mays_3_Postma_2011_day120_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df2 = pd.read_csv("../results/wheat_Morandage_day215_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df_ = [df1, df2]
depth_ = df1["depth"].loc[:].values
depth = np.unique(depth_)[::-1] 

plant = ['Maize', 'Winter wheat']
num = ['(a)', '(b)']
xtext = [0.4, 0.3, 0.07]
xlim = [1,4]
ylim = [3,12.5]
c = ['r', 'b', 'g', 'c', 'y','m']

fig, axs = plt.subplots(2,3)
for i in range(0,len(df_)): 
    df = df_[i]

    for j in range(0, len(depth)): 
        axs[0,i].set_ylim(0, ylim[i]) 
        axs[0,i].set_xlim(0, xlim[i]) 
        data = df[(df['depth']==depth[j])]
        x = data['vrld_sc_ip'].loc[:].values
        y = data['prld_stand'].loc[:].values
        axs[0,i].scatter(y,x,marker = '.', color = c[j], alpha = 0.1, edgecolor = 'None')
        axs[0,i].scatter(np.nan, np.nan, color = c[j], label = str(int(depth[j]))+' cm')
    axs[0,i].errorbar(np.nan, np.nan, np.nan, np.nan, color = 'k',marker = 'o',  label = 'Experiment')
    axs[0,i].scatter(np.nan, np.nan, color = 'k', marker = '.', label = 'Simulation')
    axs[0,i].set_xlabel('pRLD $(cm$ $cm^{-2})$')
    axs[0,i].set_ylabel('vRLD $(cm$ $cm^{-3})$')
    axs[0,i].text(0.1*xlim[i], (0.85)*ylim[i], num[i]+' '+plant[i])
    if i == 0: 
        axs[0,i].legend(bbox_to_anchor=(3.25,1))
    
for i in range(0,3): 
    axs[1,i].axis('off') 
axs[0,2].axis('off') 

#experimental data 
df1 = pd.read_csv("experimental_data/rhizoimgs_2020_planar.csv") 
df2 = pd.read_csv("experimental_data/rhizoimgs_2021_planar.csv") 
df3 = pd.read_csv("experimental_data/soil_cores_2020_lower.csv")
df4 = pd.read_csv("experimental_data/soil_cores_2021_lower.csv")

facility_ = df1["facility"].loc[:].values
plot_ = df1["plot"].loc[:].values
depth_ = df2["Depth"].loc[:].values
depth = np.unique(depth_)[::-1] 
facility = np.unique(facility_) 
plot = np.unique(plot_) 
marker = ['o','x']
cols = ['r','g','b']

data_d = []
for i in range(0, len(depth)): 
    data = []
    for j in range(0, len(plot)): 
        for k in range(0, len(facility)): 
            
            data3 = df3[(df3['Depth']==depth[i]) & (df3['plot']==plot[j])]
            if np.shape(data3)[0] != 0: 
                vrld_mean1 = np.mean(data3['RLD']) 
                vrld_SE1 = np.std(data3['RLD']) / (np.shape(data3)[0])**0.5
  
                data1 = df1[(df1['Depth']==depth[i]) & (df1['plot']==plot[j]) & (df1['facility']==facility[k])]
                prld_mean1 = np.mean(data1['pRLD']) 
                prld_SE1 = np.std(data1['pRLD']) / (np.shape(data1)[0])**0.5
            else: 
                vrld_mean1 = np.nan
                vrld_SE1 = np.nan
                prld_mean1 = np.nan
                prld_SE1 = np.nan
            
            data4 = df4[(df4['Depth']==depth[i]) & (df4['plot']==plot[j])]
            data2 = df2[(df2['Depth']==depth[i]) & (df2['plot']==plot[j]) & (df2['facility']==facility[k])]
            if np.shape(data4)[0] != 0: 
                vrld_mean2 = np.mean(data4['RLD']) 
                vrld_SE2 = np.std(data4['RLD']) / (np.shape(data4)[0])**0.5
                
                prld_mean2 = np.mean(data2['pRLD']) 
                if ~np.isnan(prld_mean2): 
                    prld_SE2 = np.std(data2['pRLD']) / (np.shape(data2)[0])**0.5
                else: 
                    prld_SE2 = np.nan
            else: 
                vrld_mean2 = np.nan
                vrld_SE2 = np.nan
                prld_mean2 = np.nan
                prld_SE2 = np.nan

            data.append([plot[j],k,prld_mean1, prld_SE1, vrld_mean1, vrld_SE1, prld_mean2, prld_SE2, vrld_mean2, vrld_SE2,depth[i]])
    data = np.array(data) 
    data_d.append(data) 
data_d = np.array(data_d)  

######PLOT########
for i in range(0, len(depth)-1): 
    data = data_d[i,:,:]
    mark = data[:,1]
    mark = mark.astype(int)
    ind_col = np.array(data[:,0])-1
    ind_col = ind_col.astype(int) 

    col = np.array(cols)[ind_col]; col_1 = col[mark==0]; col_2 = col[mark==1]
    X1 = data[:,2]; X1_1 = X1[mark==0]; X1_2 = X1[mark==1]
    X1_err = data[:,3]; X1_1err = X1_err[mark==0]; X1_2err = X1_err[mark==1]
    Y1 = data[:,4]; Y1_1 = Y1[mark==0]; Y1_2 = Y1[mark==1]
    Y1_err = data[:,5]; Y1_1err = Y1_err[mark==0]; Y1_2err = Y1_err[mark==1]
    X2 = data[:,6]; X2_1 = X2[mark==0]; X2_2 = X2[mark==1]
    X2_err = data[:,7]; X2_1err = X2_err[mark==0]; X2_2err = X2_err[mark==1]
    Y2 = data[:,8]; Y2_1 = Y2[mark==0]; Y2_2 = Y2[mark==1]
    Y2_err = data[:,9]; Y2_1err = Y2_err[mark==0]; Y2_2err = Y2_err[mark==1]

    #plot experimental data 
    for j in range(0, len(X2_2)): 
        axs[0,0].errorbar(X1_1[j], Y1_1[j], xerr = X1_1err[j], yerr=Y1_1err[j], ecolor = c[i], fmt=marker[0], mfc=c[i], mec=c[i], ms=10)
        axs[0,1].errorbar(X2_1[j], Y2_1[j], xerr = X2_1err[j], yerr=Y2_1err[j], ecolor = c[i], fmt=marker[0], mfc=c[i], mec=c[i], ms=10)
        axs[0,0].errorbar(X1_2[j], Y1_2[j], xerr = X1_2err[j], yerr=Y1_2err[j], ecolor = c[i], fmt=marker[0], mfc=c[i], mec=c[i], ms=10)
        axs[0,1].errorbar(X2_2[j], Y2_2[j], xerr = X2_2err[j], yerr=Y2_2err[j], ecolor = c[i], fmt=marker[0], mfc=c[i], mec=c[i], ms=10)

plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

