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


mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams.update({'font.size': 12})

left  = 0.02  # the left side of the subplots of the figure
right = 0.93    # the right side of the subplots of the figure
bottom = 0.1   # the bottom of the subplots of the figure
top = 0.95      # the top of the subplots of the figure
wspace = 0.3   # the amount of width reserved for blank space between subplots
hspace = 0.3   # the amount of height reserved for white space between subplots

day = 7


npzfile = np.load('results/maize_high_resolution_hydrus_loam_day'+str(day)+'.npz')
wc =npzfile['arr_0'] #soil water content 
C_frac =npzfile['arr_1'] #root volume fraction 
C_rootvol =npzfile['arr_2'] #root volume
cell_number = npzfile['arr_3'] #cell number 

npzfile = np.load('results/maize_high_resolution_hydrus_loam_stitched_day'+str(day)+'.npz')
wc_stitch =npzfile['arr_0'] #soil water content 
C_frac_stitch =npzfile['arr_1'] #root volume fraction 
C_rootvol_stitch =npzfile['arr_2'] #root volume
cell_number_stitch = npzfile['arr_3'] #cell number 


"""visualisation of last soil + root water content to check if correct"""  
fig = plt.figure(figsize=plt.figaspect(0.5))
cmap = plt.cm.jet_r
cmap_name = 'jet_r'

#single plant, x direction
ax = fig.add_subplot(1, 4, 1)
# srwc = wc[:,int(cell_number[1]/2),:] 
srwc = C_rootvol[int(cell_number[0]/2),:,:]*0.8 +wc[int(cell_number[0]/2),:,:]  #assumption: mean root water content = 0.8
srwc = np.swapaxes(srwc,0,1)
srwc = srwc[::-1,:]
colors = cmap(srwc)
norm = mpl.colors.Normalize(vmin=np.min(srwc), vmax=np.max(srwc)-0.04)
ax.imshow(srwc, cmap = cmap_name, norm = norm) #1.67
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax, label = 'Root and soil water content (-)',shrink=0.5)
ax.set_title("RWC (-) and $\\theta$ (-), single plant,\n x direction")
ax.axis('off')

#single plant, y direction
ax = fig.add_subplot(1, 4, 2)
srwc = C_rootvol[:,int(cell_number[1]/2),:]*0.8 +wc[:,int(cell_number[1]/2),:]  #assumption: mean root water content = 0.8
srwc = np.swapaxes(srwc,0,1)
srwc = srwc[::-1,:]
colors = cmap(srwc)
norm = mpl.colors.Normalize(vmin=np.min(srwc), vmax=np.max(srwc)-0.04)
ax.imshow(srwc, cmap = cmap_name, norm = norm) #1.67
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax, label = 'Root and soil water content (-)',shrink=0.5)
ax.set_title("RWC (-) and $\\theta$ (-), single plant,\n y direction")
ax.axis('off')

#stitched, x direction
ax = fig.add_subplot(1, 4, 3)
srwc_stitch = C_rootvol_stitch[int(cell_number_stitch[0]/4),:,:]*0.8 + wc_stitch[int(cell_number_stitch[0]/4),:,:] #assumption: mean root water content = 0.8
srwc_stitch = np.swapaxes(srwc_stitch,0,1)
srwc_stitch = srwc_stitch[::-1,:]
colors = cmap(srwc_stitch)
# norm = mpl.colors.Normalize(vmin=np.min(srwc_stitch), vmax=np.max(srwc_stitch))
ax.imshow(srwc_stitch, cmap = cmap_name, norm = norm, aspect=1) #1.67
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax, label = 'Root and soil water content (-)',shrink=0.5)
ax.set_title("RWC (-) and $\\theta$ (-), stitched,\n x direction")
ax.axis('off')

#stitched, y direction
ax = fig.add_subplot(1, 4, 4)
srwc_stitch = C_rootvol_stitch[:,int(cell_number_stitch[1]/2),:]*0.8 + wc_stitch[:,int(cell_number_stitch[1]/2),:] #assumption: mean root water content = 0.8
srwc_stitch = np.swapaxes(srwc_stitch,0,1)
srwc_stitch = srwc_stitch[::-1,:]
colors = cmap(srwc_stitch)
# norm = mpl.colors.Normalize(vmin=np.min(srwc_stitch), vmax=np.max(srwc_stitch))
ax.imshow(srwc_stitch, cmap = cmap_name, norm = norm, aspect=1) #1.67
m = cm.ScalarMappable(cmap=cmap, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax, label = 'Root and soil water content (-)',shrink=0.5)
ax.set_title("RWC (-) and $\\theta$ (-), stitched,\n y direction")
ax.axis('off')


plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
plt.show()