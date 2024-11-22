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
import csv

mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams.update({'font.size': 18})

npzfile = np.load('sims_root_soil.npz')
swc =npzfile['arr_0']
rwc =npzfile['arr_1']
cn = npzfile['arr_2']
name = npzfile['arr_3']
min_b = npzfile['arr_4']
max_b = npzfile['arr_5']
sim_time = npzfile['arr_6']
rs_age = npzfile['arr_7']

""" resimulate root system """
path = "../../modelparameter/structural/rootsystem/"
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

#PLOT
fig = plt.figure(figsize=plt.figaspect(0.5))

#plot root system 
ax1 = fig.add_subplot(1, 3, 1)
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
ax2 = fig.add_subplot(1, 3, 2)
swc_ = swc[:,int(cn[1]/2),:]
colors = plt.cm.seismic_r(swc_)
norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(swc_))
ax2.imshow(swc_, cmap = 'seismic_r', norm = norm, aspect=1)#1.67
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax2, label = 'Soil water content (-)')
ax2.set_title("$\\theta$ (-)")
ax2.axis('off')

#visualise root volume fraction
ax3 = fig.add_subplot(1, 3, 3)
rwc_ = rwc[:,int(cn[1]/2),:]
colors = plt.cm.seismic_r(rwc_)
norm = matplotlib.colors.Normalize(vmin=0, vmax=np.max(rwc_))
ax3.imshow(rwc_, cmap = 'seismic_r', norm = norm, aspect=1) #1.67
m = cm.ScalarMappable(cmap=plt.cm.seismic_r, norm=norm)
m.set_array([])
plt.colorbar(m, ax = ax3, label = 'Root volume fraction (-)')
ax3.set_title("RVF")
ax3.axis('off')

plt.show()
