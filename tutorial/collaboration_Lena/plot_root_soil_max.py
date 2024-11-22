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

matplotlib.rcParams.update({'font.size': 18})

#root system
#""" Parameters """
path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI"
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
sim_time = 1  # [day] for task b
rs_age = 10  # root system initial age

""" Initialize xylem model """
rs = pb.MappedRootSystem()
rs.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
rs.setGeometry(sdf)
rs.setSeed(1)
rs.initialize()
rs.simulate(rs_age, False)
for i in range(0, sim_time):
    rs.simulate(i)
    
ana = pb.SegmentAnalyser(rs.mappedSegments())
segs = ana.segments
nodes = ana.nodes
radius = ana.getParameter("radius")
seglen = ana.getParameter("length")

npzfile = np.load('sims_root_soil_fine.npz')
swc =npzfile['arr_0']
rwc =npzfile['arr_1']
cn = npzfile['arr_2']

#PLOT
fig = plt.figure(figsize=plt.figaspect(0.5))

#plot root system 
ax1 = fig.add_subplot(1, 4, 1)
cmap = cm.seismic
norm = matplotlib.colors.Normalize(vmin=np.min(radius), vmax=np.max(radius))
fc = cmap(norm(radius))
for k, s in enumerate(segs):
    s1 = segs[k]
    n1, n2 = nodes[s1.x], nodes[s1.y]
    ax1.plot([n1.x, n2.x], [n1.z, n2.z], color = fc[k,:])
m = cm.ScalarMappable(cmap=plt.cm.seismic)
m.set_array([])
plt.colorbar(m, ax = ax1, label = 'Radius (cm)')
ax1.set_aspect('equal')
ax1.axis('off')


#visualise soil water content
ax2 = fig.add_subplot(1, 4, 2)
swc_med_ = np.median(swc)
swc_med = np.ones(np.shape(swc))*swc_med_
swc_ = np.amax(np.absolute(swc-swc_med_),1)
colors = plt.cm.seismic_r(swc_)
#norm = matplotlib.colors.Normalize(vmin=np.min(swc_), vmax=np.max(swc_))
norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(swc_))
ax2.imshow(swc_, cmap = 'seismic_r', norm = norm, aspect = 1.67)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
#m.set_array([])
#plt.colorbar(m, ax = ax2, label = 'Water content (-)')
ax2.set_title("Soil")
ax2.axis('off')

#visualise root water content
ax3 = fig.add_subplot(1, 4, 3)
rwc_med_ = np.median(rwc)
rwc_med = np.ones(np.shape(rwc))*rwc_med_
rwc_ = np.amax(np.absolute(rwc-rwc_med_),1)
colors = plt.cm.seismic_r(rwc_)
#norm = matplotlib.colors.Normalize(vmin=np.min(rwc_), vmax=np.max(rwc_))
norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(swc_))
ax3.imshow(rwc_, cmap = 'seismic_r', norm = norm, aspect = 1.67)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
#m.set_array([])
#plt.colorbar(m, ax = ax3, label = 'Water content (-)')
ax3.set_title("Root")
ax3.axis('off')

#visualise soil + root water content
srwc = swc + rwc
ax4 = fig.add_subplot(1, 4, 4)
srwc_med_ = np.median(srwc)
srwc_med = np.ones(np.shape(srwc))*srwc_med_
srwc_ = np.amax(np.absolute(srwc-srwc_med_),1)
colors = plt.cm.seismic_r(srwc_)
#norm = matplotlib.colors.Normalize(vmin=np.min(srwc_), vmax=np.max(srwc_))
norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(swc_))
ax4.imshow(srwc_, cmap = 'seismic_r',norm = norm, aspect = 1.67)
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax4, label = 'max deviation in water content (-)')
ax4.set_title("Soil and root")
ax4.axis('off')

plt.show()
