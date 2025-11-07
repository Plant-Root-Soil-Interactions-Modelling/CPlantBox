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
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import scipy.stats as stats

font = {'size'   : 16}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'

left  = 0.08  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
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
xlim = [1,5]
ylim = [1,5]
c = ['r', 'b', 'g', 'c', 'y','m']
marker = ['x', '.']
alts = ['2x2','6x4', '2x20']
labels = ['2x2 cm', '6x4 cm', 'Entire tube surface']
annot = ['(a)', '(b)']

fig, axs = plt.subplots(2,3)
for i in range(0,len(df_)): 
    R2 = np.zeros((len(alts)))
    for k in range(0, len(alts)): 
        
        if i == 0: 
            df = pd.read_csv("../results/Zea_mays_3_Postma_2011_day120_reps100_imgsize"+alts[k]+"tubediam6.4_pdensity_base.csv")
        else: 
            df = pd.read_csv("../results/wheat_Morandage_day215_reps100_imgsize"+alts[k]+"tubediam6.4_pdensity_base.csv")
    
        axs[0,i].set_ylim(0, ylim[i]) 
        axs[0,i].set_xlim(0, xlim[i]) 
        data = df
        x = data['prld_stand'].loc[:].values
        y = data['vrld_sl'].loc[:].values
        popt, pcov = curve_fit(func, x, y)
        axs[0,i].scatter(x,y,marker = '.', color = c[k], alpha = 0.05, edgecolor = 'None')
        axs[0,i].plot(x, func(x, popt),color = c[k], linestyle = ':', linewidth = 2)
        
        #stats
        popt, pcov = curve_fit(func, x, y)
        y_pred = func(x, *popt)
        R2[k]= r2_score(y, y_pred)
    
    axs[0,i].scatter(1000, 1000, marker = '.', color = c[0], label = labels[0]+', $R^2$ = '+str(np.around(R2[0],2)))
    axs[0,i].scatter(1000, 1000, marker = '.', color = c[1],  label = labels[1]+', $R^2$ = '+str(np.around(R2[1],2)))
    axs[0,i].scatter(1000, 1000, marker = '.', color = c[2],  label = labels[2]+', \n$R^2$ = '+str(np.around(R2[2],2)))
    axs[0,i].plot([0,1000], [0,1000], 'k--')
    axs[0,i].set_ylabel('vRLD $(cm$ $cm^{-3})$')
    axs[0,i].set_title(plant[i])
    t = axs[0,i].text(0.9*xlim[i], 0.9*ylim[i], annot[i])
    t.set_bbox(dict(facecolor='w', edgecolor='w'))

    axs[0,i].legend()
    
for i in range(0,3): 
    axs[1,i].axis('off') 
axs[0,2].axis('off') 
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

