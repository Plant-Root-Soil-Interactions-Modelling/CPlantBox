import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

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

def vg_model(h, theta_r, alpha, n, theta_s):
    """
    h       : pressure head [cm] (positive values expected)
    theta_r : residual water content [cm3/cm3]
    theta_s : saturated water content [cm3/cm3]
    alpha   : [1/cm]
    n       : [-]
    """
    m = 1 - 1/n
    Se = (1 + (alpha * np.abs(h))**n)**(-m)
    return theta_r + (theta_s - theta_r) * Se
    
def vg_fit_full(h, theta_r, alpha, n, theta_s):
    return vg_model(h, theta_r, alpha, n, theta_s)

def vg_fit_fixed(h, theta_r, alpha, n):
    return vg_model(h, theta_r, alpha, n, theta_s)


# label_soil = [['0-20 cm', '20-75 cm'],['0-75cm']]
soils = ['loam', 'sand']
cols = ['r', 'm', 'b', 'c']
annots = ['(a)', '(b)']
data_Helena = ['data/retention_loam.csv', 'data/retention_sand.csv']
theta_s_ = [0.411, 0.337]
tau_ = [ 0.5, 0.5]
Ks_ = [441, 245]

fig, ax = plt.subplots(1,2,figsize=(16,8))
for j in range(0, len(soils)): 
    
    theta_s = theta_s_[j]
    tau = tau_[j]
    Ks = Ks_[j]
        
    data = np.loadtxt(data_Helena[j], delimiter=",")
    depth_levels = ['10 cm', '20 cm', '40 cm', '60 cm']
    depths = data[:, 0] 
    d_unique = np.unique(depths)
    h = 10**data[:, 1]     
    theta = data[:, 2] 

    if soils[j] == 'loam': 
        layers = 1
        layer_depths = [[1,2,3,4]]
    else: 
        layers = 1
        layer_depths = [[1,2,3,4]]

    for i in range(0, layers): 
        mask = np.isin(depths, layer_depths[i])
        h_mask = h[mask]  
        theta_mask = theta[mask]

    if j == 0:  # loam

        theta_r0 = 0.001
        alpha0 = 0.009
        n0 = 1.5
        theta_s0 = theta_s

        p0 = [theta_r0, alpha0, n0, theta_s0]

        bounds = (
            [0.0, 1e-4, 1.01, 0.411],   # lower
            [0.1, 0.05,  3.0,  0.60]   # upper
        )

        popt, pcov = curve_fit(
            vg_fit_full,
            h_mask,
            theta_mask,
            p0=p0,
            bounds=bounds,
            maxfev=10000
        )

        theta_r, alpha, n, theta_s_fit = popt

        print('loam')
        print('theta_r=', theta_r)
        print('alpha=', alpha)
        print('n=', n)
        print('theta_s=', theta_s_fit)

        h_fit = np.logspace(-1, 5, 200)
        theta_fit = vg_fit_full(h_fit, *popt)
       
    else:  # sand

        theta_r0 = 0.0
        alpha0 = 0.1
        n0 = 2.1

        p0 = [theta_r0, alpha0, n0]

        bounds = (
            [0.0, 1e-3, 1.0],   # lower
            [0.1, 0.5,  5.0]    # upper
        )

        popt, pcov = curve_fit(
            vg_fit_fixed,
            h_mask,
            theta_mask,
            p0=p0,
            bounds=bounds,
            maxfev=10000
        )

        theta_r, alpha, n = popt

        print('sand')
        print('theta_r=', theta_r)
        print('alpha=', alpha)
        print('n=', n)

        h_fit = np.logspace(-1, 5, 200)
        theta_fit = vg_fit_fixed(h_fit, *popt)
        
    ax[j].semilogx(h_fit, theta_fit, color = 'k')
    
    #plot Helenas data points 
    for m, depth in enumerate(d_unique):
        mask = depths == depth
        ax[j].scatter(h[mask], theta[mask], color = cols[m], label=depth_levels[m])     
    
    ax[j].set_xlim([10**-2, 10**5])
    ax[j].set_ylim([0, 0.45])
    ax[j].set_xlabel("Soil water potential (cm)")
    ax[j].set_ylabel("Water content (cm³/cm³)")
    ax[j].grid()
    ax[j].annotate(annots[j], xy=get_axis_limits(ax[j]))
    
    # #make 2 legends
    # lines = ax[j].get_lines()
    # if soils[j] == 'loam': 
        # legend1 = ax[j].legend([lines[i] for i in np.arange(0,1)], ['Optimized retention curve'],loc='lower left')
    # # else: 
        # # legend1 = ax[j].legend([lines[i] for i in np.arange(0,1)], ['0-75 cm'], loc='upper left')
        # ax[j].add_artist(legend1)
    
    # if j == 0: 
        # points = ax[j].collections
        # legend2 = ax[j].legend([points[i] for i in np.arange(0,4)], ['10 cm', '20 cm', '40 cm', '60 cm'], loc='center left')
        # ax[j].add_artist(legend2)
    
    
    
    
plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()
# plt.savefig('Jorda_retention_test', transparent=None, dpi='figure', format=None,
        # metadata=None, bbox_inches=None, pad_inches=0.1,
        # facecolor='auto', edgecolor='auto', backend=None)