""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import functional.van_genuchten as vg
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import timeit
from scipy import interpolate 
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors
from PIL import Image
import os
from pathlib import Path
from os import listdir
from os.path import isfile, join
import re

mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams.update({'font.size': 25})

left  = 0.1  # the left side of the subplots of the figure
right = 0.9    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.8    # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots

def params(plant):
    if plant == 'maize': 
        min_b = [-75/2, -15/2, -150] #[cm]
        max_b = [75/2, 15/2, 0.]  #[cm]
        param_name = "Zeamays_synMRI_modified.xml"
    elif plant == 'wheat': 
        min_b = [-15/2, -3/2, -150] #[cm]
        max_b = [15/2, 3/2, 0.]  #[cm]
        param_name = "wheat_Morandage.xml"
    else: 
        print('wrong plant')
    return min_b, max_b, param_name
    
    
    
scenarios = os.listdir('results/')
num_scenarios = len(scenarios)
day = 8
evap = False
regex = re.compile(r'\d+')
for i in range(0, len(scenarios)): 
    info = scenarios[i].split('_')
    simtime = int(regex.findall(scenarios[i])[0])
    plant = info[1]
    soil = info[5]
    if len(info)>6: 
        evap = True
    min_b, max_b, param_name = params(plant)

    fname = os.listdir('results/'+scenarios[i])[-1]
    npzfile = np.load('results/'+scenarios[i]+'/'+fname)
    wc_stitch =npzfile['arr_0'] #soil water content 
    C_frac_stitch =npzfile['arr_2'] #root volume fraction 
    C_rootvol_stitch =npzfile['arr_3'] #root volume
    cell_number_stitch = npzfile['arr_4'] #cell number 

    """visualisation of root systems and soil + root water content in x-y direction"""  
    # fig = plt.figure(figsize=plt.figaspect(0.5))
    fig, (ax2, ax1) = plt.subplots(ncols=2, sharex=True)
    cmap = plt.cm.jet_r
    cmap_name = 'jet_r'


    #stitched, y direction
    srwc_stitch = C_rootvol_stitch[:,int(cell_number_stitch[1]/2),:]*0.8 + wc_stitch[:,int(cell_number_stitch[1]/2),:] #assumption: mean root water content = 0.8
    srwc_stitch = np.swapaxes(srwc_stitch,0,1)
    srwc_stitch = srwc_stitch[::-1,:]
    colors = cmap(srwc_stitch)
    norm = mpl.colors.Normalize(vmin=np.min(srwc_stitch), vmax=np.max(srwc_stitch))
    ax1.imshow(srwc_stitch, cmap = cmap_name, norm = norm, aspect=1) #1.67
    m = cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax = ax1, label = 'Root and soil water content (-)',shrink=0.5)
    ax1.set_title("y direction")
    ax1.axis('off')
    
    #stitched, x direction
    srwc_stitch = C_rootvol_stitch[int(cell_number_stitch[0]/4),:,:]*0.8 + wc_stitch[int(cell_number_stitch[0]/4),:,:] #assumption: mean root water content = 0.8
    srwc_stitch = np.swapaxes(srwc_stitch,0,1)
    srwc_stitch = srwc_stitch[::-1,:]
    colors = cmap(srwc_stitch)
    # norm = mpl.colors.Normalize(vmin=np.min(srwc_stitch), vmax=np.max(srwc_stitch))
    ax2.imshow(srwc_stitch, cmap = cmap_name, norm = norm, aspect=1) #1.67
    m = cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array([])
    plt.colorbar(m, ax = ax2, label = 'Root and soil water content (-)',shrink=0.5)
    ax2.set_title("x direction")
    ax2.axis('off')

    if evap: 
        title = plant+', '+soil+', with evaporation'
    else: 
        title = plant+', '+soil+', '
    plt.suptitle(title)
    
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    plt.show()