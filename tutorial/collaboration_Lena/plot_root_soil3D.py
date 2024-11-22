import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.root_conductivities import *  # hard coded conductivities
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import math
import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors
import timeit

npzfile = np.load('sims_root_soil2.npz')
swc =npzfile['arr_0']
rwc =npzfile['arr_1']
cn = npzfile['arr_2']

#PLOT
fig = plt.figure(figsize=plt.figaspect(0.5))

#visualise soil water content
ax1 = fig.add_subplot(1, 3, 1, projection='3d')
swc = np.swapaxes(swc, 0, 2)
colors = plt.cm.seismic_r(swc)
colors[:,:,:,3] = 0.5
print(np.min(swc), np.max(swc))
norm = matplotlib.colors.Normalize(vmin=np.min(swc), vmax=np.max(swc))
ax1.voxels(swc, facecolors=colors, edgecolor='none',norm = norm)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax1, label = 'Water content (-)')
ax1.set_title("Soil water content")
ax1.axis('off')

#visualise root water content
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
rwc = np.swapaxes(rwc, 0, 2)
colors = plt.cm.seismic_r(rwc)
colors[:,:,:,3] = 0.5
norm = matplotlib.colors.Normalize(vmin=np.min(rwc), vmax=np.max(rwc))
ax2.voxels(rwc, facecolors=colors, edgecolor='none',norm = norm)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax2, label = 'Water content (-)')
ax2.set_title("Root water content")
ax2.axis('off')

#visualise soil + root water content
srwc = swc + rwc
ax3 = fig.add_subplot(1, 3, 3, projection='3d')
colors = plt.cm.seismic_r(srwc)
colors[:,:,:,3] = 0.5
norm = matplotlib.colors.Normalize(vmin=np.min(srwc), vmax=np.max(srwc))
ax3.voxels(srwc, facecolors=colors, edgecolor='none',norm = norm)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax3, label = 'Water content (-)')
ax3.set_title("Combined soil and root water content")
ax3.axis('off')

plt.show()
