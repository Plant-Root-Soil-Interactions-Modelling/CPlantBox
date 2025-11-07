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

font = {'size'   : 20}
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

fig, axs = plt.subplots(2,3)
for i in range(0,len(df_)): 
    df = df_[i]
    x = df['vrld_sl'].loc[:].values
    y = df['prld_highr'].loc[:].values
    popt, pcov = curve_fit(func, x, y)
    y_pred = func(x, *popt)
    R2= r2_score(y, y_pred)
    
    
    # Statistics
    n = y.size                        
    m = popt.size                     
    dof = n - m                                   
    t = stats.t.ppf(0.975, n - m)                       
    
    # Estimates of Error in Data/Model
    resid = y - y_pred                                       
    s_err = np.sqrt(np.sum(resid**2) / dof)
    
    x2 = np.linspace(np.min(x), np.max(x), 100)
    y2 = func(x2, popt)
    pi = t * s_err * np.sqrt(1 + 1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))   
    axs[0,i].plot(func(x, popt),x, color = "0.3", linestyle = ':', linewidth = 1) 
    axs[0,i].plot(np.nan, np.nan,color = '0.3', linestyle = ':', label = "Regression line \nfor all depths")
    axs[0,i].plot([0,1000], [0,1000], 'k--', label = '1:1 line')
    slope = popt[0]
    print(plant[i], np.around(1/popt[0],2), np.around(pi[0], 2), np.around(R2,2)) 
    
    for j in range(0, len(depth)): 
        axs[0,i].set_ylim(0, ylim[i]) 
        axs[0,i].set_xlim(0, xlim[i]) 
        data = df[(df['depth']==depth[j])]
        x = data['vrld_sl'].loc[:].values
        y = data['prld_stand'].loc[:].values
        popt, pcov = curve_fit(func, x, y)
        axs[0,i].scatter(y,x,marker = '.', color = c[j], alpha = 0.1, edgecolor = 'None')
        axs[0,i].plot(func(x, popt),x, color = c[j], linestyle = ':', linewidth = 2) 
        
        y_pred = func(x, *popt)
        R2= r2_score(y, y_pred)
        print(plant[i], depth[j], np.around(1/popt[0],2), np.around(R2,2)) 
        axs[0,i].scatter(np.nan, np.nan, color = c[j], label = str(int(depth[j]))+' cm')


    axs[0,i].set_xlabel('pRLD $(cm$ $cm^{-2})$')
    axs[0,i].set_ylabel('vRLD $(cm$ $cm^{-3})$')
    axs[0,i].text(0.1*xlim[i], (0.85)*ylim[i], num[i]+' '+plant[i])
    if i == 0: 
        axs[0,i].legend(bbox_to_anchor=(3.15,1.05))
    
for i in range(0,3): 
    axs[1,i].axis('off') 
axs[0,2].axis('off') 
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()

