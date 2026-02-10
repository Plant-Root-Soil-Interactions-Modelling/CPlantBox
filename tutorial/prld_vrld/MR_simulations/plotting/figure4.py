import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import pandas as pd
from sklearn.preprocessing import StandardScaler
import sys
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt
import seaborn as sns

font = {'size'   : 20}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'

left  = 0.08  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.4   # the amount of height reserved for white space between subplots

df1 = pd.read_csv("../results/Zea_mays_3_Postma_2011_day120_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df2 = pd.read_csv("../results/wheat_Morandage_day215_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df_ = [df1, df2]
depth_ = df1["depth"].loc[:].values
depth = np.unique(depth_)[::-1] 
plant = ['Maize', 'Winter wheat']
num = ['(a)', '(b)']
cols = ['r', 'b', 'g', 'c', 'y','m', 'k']
annot = [['(a)', '(b)'], ['  (c)', '(d)'], ['(e)', '(f)']]

fac = ['rRLD', 'An', 'CV']
fac_new = ['CV of rRLD', 'CV of An']
fac_all = ['CV of rRLD', 'CV of An', 'CV']
fig, ax = plt.subplots(3,3) 
for i in range(0, len(df_)): 
    df = df_[i]
    
    target = df['vrld_sl']/df['prld_cont']
    df_corr = pd.DataFrame({"vRLD/pRLD":target})
    for k in range(0,2): 
        c = df[fac[k]].loc[:].values
        c_CV = np.std(c.reshape(-1, len(depth)), axis=1)/np.mean(c.reshape(-1, len(depth)), axis=1)
        c_CV_ = np.repeat(c_CV, len(depth))
        df_new = pd.DataFrame({fac_new[k]:c_CV_})
        df_corr = pd.concat([df_corr, df_new[fac_new[k]]], axis = 1)
    df_corr = pd.concat([df_corr, df["CV"]], axis = 1)
    df_corr_depth = pd.concat([df_corr, df["depth"]], axis = 1)
    
    for m in range(0, len(fac_all)): 
        dummymaxx = 0
        dummyminx = 0
        dummyy = 0
        for j in range(0, len(depth)): 
            df_corr_ = df_corr[(df_corr_depth['depth']==depth[j])]

            # Plotting 
            x = df_corr_[fac_all[m]].loc[:].values
            y = df_corr_['vRLD/pRLD'].loc[:].values
            z = np.polyfit(x, y, 1)
            p = np.poly1d(z)
            ax[m,i].scatter(x,y,marker = '.', color = cols[j], alpha = 0.1, edgecolor = 'None')
            ax[m,i].scatter(np.nan, np.nan, color = cols[j], label = str(int(depth[j]))+' cm')
            ax[m,i].plot(x, p(x), color = cols[j], linestyle = ':', linewidth = 2)
            ax[m,i].set_xlabel(fac_all[m]+' $(-)$')
            ax[m,i].set_ylabel('vRLD/pRLD $(cm^{-1})$')
            ax[m,i].set_ylim([0,3])
            if i == 0: 
                if m==0: 
                    ax[m,i].legend(bbox_to_anchor=(2.9,1.05))
            
            if m==0: 
                ax[m,i].set_title(plant[i])
            xmin, xmax = ax[m,i].get_xlim()
            ymin, ymax = ax[m,i].get_ylim()
            if xmax>dummymaxx: 
                dummymaxx= xmax
            if xmin>dummyminx: 
                dummyminx = xmin
            if ymax>dummyy: 
                dummyy= ymax
            if j == (len(depth)-1): 
                t = ax[m,i].text(0.8*(dummymaxx-dummyminx)+dummyminx, 0.8*dummyy, annot[m][i])
                t.set_bbox(dict(facecolor='w', edgecolor='w'))

for i in range(0,3): 
    ax[i,2].axis('off') 
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()
    
