import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import pandas as pd

fontsize = 20
mpl.rcParams.update({
    'font.size': fontsize,
    'axes.titlesize': fontsize,
    'axes.labelsize': fontsize,
    'xtick.labelsize': fontsize,
    'ytick.labelsize': fontsize,
    'legend.fontsize': fontsize,
    'figure.titlesize': fontsize,
    'mathtext.default': 'regular'
})

left  = 0.1  # the left side of the subplots of the figure
right = 0.95    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.25   # the amount of width reserved for blank space between subplots
hspace = 0.5   # the amount of height reserved for white space between subplots

def get_axis_limits(ax, scalex=0.1, scaley=0.9):
    return ax.get_xlim()[1]*scalex, ax.get_ylim()[1]*scaley
    
def van_genuchten(h, theta_r, theta_s, alpha, n):

    m = 1 - 1 / n

    return theta_r + (
        (theta_s - theta_r) /
        (1 + (alpha * np.abs(h))**n)**m
    )

    
def residuals(params):
    theta_r, theta_s, alpha, n = params

    theta_pred = van_genuchten(
        h_fit,
        theta_r,
        theta_s,
        alpha,
        n
    )

    return theta_pred - theta_fit



soils = ['loam', 'sand']
cols = ['g', 'y', 'b', 'c']
# cols = ['blue', 'cornflowerblue']
marker = ['o', 'x', 'd', '*']
annots = ['(a)', '(b)']
path2exp = "data/"
depth_points = np.array([-10, -20, -40, -60])
exp_labels = [10,20,40,60]
theta_s_ = [0.411, 0.337]
tau_ = [ 0.5, 0.5]
Ks_ = [441, 1174]

fig, ax = plt.subplots(1,2,figsize=(16,8))
for j in range(0, len(soils)): 
    
    theta_s = theta_s_[j]
    tau = tau_[j]
    Ks = Ks_[j]
        
    #experiment 
    df_theta = pd.read_csv(path2exp+"WC_"+soils[j]+".csv")
    time_theta = df_theta['time'].values
    theta = []
    for k in range(0, len(depth_points)): 
        theta.append(df_theta['cm'+str(exp_labels[k])].values)
    theta = np.asarray(theta).T

    df_h = pd.read_csv(path2exp+"WP_"+soils[j]+".csv")
    time_psi = df_h['time'].values
    h = []
    for k in range(0, len(depth_points)): 
        h.append(df_h['cm'+str(exp_labels[k])].values)
    h = np.asarray(h).T
    
    #mask for values < 50 days
    mask = time_psi < 50
    theta_50 = theta[mask,:]
    h_50 = h[mask,:]
    
    theta_ = theta_50[:,:2].flatten()
    h_ = h_50[:,:2].flatten()

    #mask for outliers 
    mask = (
        np.isfinite(theta_) &
        np.isfinite(h_) &
        (theta_ > 0) &
        (theta_ < 0.5)
    )
    h_fit = h_[mask]
    theta_fit = theta_[mask]
    
    if j == 0:  # loam

        theta_r0 = 0.001
        alpha0 = 0.009
        n0 = 2
        theta_s0 = theta_s
        
        x0 = [theta_r0, theta_s0, alpha0, n0]

        bounds = (
            [0.001,0.411,  1e-4, 1.],   # lower
            [0.1, 0.45, 0.05,  5.0]   # upper
        )

    else: # sand
        theta_r0 = 0.0
        alpha0 = 0.01
        n0 = 5
        theta_s0 = theta_s
        
        x0 = [theta_r0, theta_s0, alpha0, n0]

        bounds = (
            [0.0, 0.337, 1e-3, 1.],   # lower
            [0.038, 0.338, 0.1,  5.5]    # upper
        )

    result = least_squares(
        residuals,
        x0,
        bounds=bounds
    )

    theta_r, theta_s, alpha, n = result.x

    print("theta_r =", theta_r)
    print("theta_s =", theta_s)
    print("alpha   =", alpha)
    print("n       =", n)
    print("RMSE    =", np.sqrt(np.mean(result.fun**2)))

    h_curve = np.logspace(-1, 5, 200)
    theta_curve = van_genuchten(h_curve, theta_r, theta_s, alpha, n)
          
    ax[j].semilogx(h_curve, theta_curve, color = 'k')
    
    #plot data points per depth
    # ax[j].scatter(abs(h_clean), theta_clean,color = 'k')
    for k in range(0, np.shape(theta)[1]): 
        # ax[j].scatter(abs(h_50[:,k]), theta_50[:,k], color = 'r', label=str(depth_points[k])+' cm')     
        ax[j].scatter(abs(h[:,k]), theta[:,k],  color = cols[k],  label=str(depth_points[k])+' cm')  
    ax[j].scatter(abs(h_fit), theta_fit,color = 'r', facecolor = 'None', label = 'points used for fitting')
    
    ax[j].set_xlim([10**0, 10**5])
    ax[j].set_ylim([0, 0.45])
    ax[j].set_xlabel("Soil water potential (cm)")
    ax[j].set_ylabel("Water content (cm³/cm³)")
    ax[j].grid()
    ax[j].annotate(annots[j], xy=get_axis_limits(ax[j]))
    ax[j].legend()
    
    
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()
