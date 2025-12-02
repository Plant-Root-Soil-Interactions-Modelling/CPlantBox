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

font = {'size'   : 16}
plt.rc('font', **font)
mpl.rcParams['mathtext.default'] = 'regular'

left  = 0.1  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.9      # the top of the subplots of the figure
wspace = 0.1   # the amount of width reserved for blank space between subplots
hspace = 0.15   # the amount of height reserved for white space between subplots

df1 = pd.read_csv("../results/Zea_mays_3_Postma_2011_day120_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df2 = pd.read_csv("../results/wheat_Morandage_day215_reps100_imgsize2x2tubediam6.4_pdensity_base.csv")
df_ = [df1, df2]
depth_ = df1["depth"].loc[:].values
depth = np.unique(depth_)[::-1] 
plant = ['Maize', 'Winter wheat']
num = ['(a)', '(b)']

fac = ['rRLD', 'An', 'CV']
fac_new = ['CV of rRLD', 'CV of An']
fig, ax = plt.subplots(2,6) 
cbar_ax = fig.add_axes([.93,.3,.01,.3])
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
    
    for j in range(0, len(depth)): 
        df_corr_ = df_corr[(df_corr_depth['depth']==depth[j])]

        # Plotting the correlation matrix
        mask = np.triu(np.ones_like(df_corr_.corr(), dtype=bool))
        np.fill_diagonal(mask, False)
        norm = plt.Normalize(-1,1)
        ytick = False 
        xtick = True
        if j ==0: 
            ytick = True
        if i == 0: 
            xtick = True 
        annot = False
        ticks = [-1, -0.5, 0, 0.5, 1]
        axcbar = sns.heatmap(df_corr_.corr(), mask = mask, square=True, annot=annot, ax = ax[i,j], cbar_kws={"shrink": 0.1, 'ticks': ticks}, cmap='coolwarm', norm = norm,xticklabels=xtick, yticklabels=ytick, cbar_ax = cbar_ax)
        cbar = axcbar.collections[0].colorbar
        cbar.set_label('$\hat{r}$', labelpad=-50, y=1.1, rotation=0)
        
        
        ax[i,j].tick_params('x', labelrotation=45) 
        if i == 0: 
            ax[i,j].tick_params('y', labelrotation=0) 
        if j==0: 
            ax[i,j].set_title(num[i] + '  '+ plant[i]+'\n'+str(int(depth[j]))+'cm \n', loc='left')
        elif j>0: 
            ax[i,j].set_title('\n'+'\n'+str(int(depth[j]))+' cm\n', loc='left')
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)    
plt.show()
    
