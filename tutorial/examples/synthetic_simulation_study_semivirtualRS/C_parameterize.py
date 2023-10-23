import sys
sys.path.append("../../../")
sys.path.append("../../../src/");
sys.path.append("../../../gui/viewer/");
import plantbox as pb
from functional.xylem_flux import XylemFluxPython 
import visualisation.vtk_plot as vp
from viewer_data import ViewerDataModel
from visualisation.vtk_tools import *
import rsml.rsml_writer
import rsml.rsml_reader as rsml
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import os
import math
import random
#########################################################
params = {'legend.fontsize': 18,
          'axes.labelsize': 18,
          'xtick.labelsize': 18,
          'ytick.labelsize': 18,
          'legend.handlelength': 2}
plt.rcParams.update(params)


path  = 'results/'
fname = 'Faba_day10'

#read in the rsml scaffold 
rs1 = XylemFluxPython.read_rsml(path + 'GT_'+ fname +".rsml")
polylines1, props1, funcs1, _ = rsml.read_rsml(path +'GT_'+  fname +".rsml")
rs2 = XylemFluxPython.read_rsml(path +'REC_'+ fname + ".rsml")
polylines2, props2, funcs2, _ = rsml.read_rsml(path + 'REC_'+ fname +".rsml")

ana1 = pb.SegmentAnalyser(rs1)
ana2 = pb.SegmentAnalyser(rs2)

#determine the missing root length
radii1 = ana1.data["radius"]
weights1 = [ana1.getSegmentLength(i) for i in range(0, len(ana1.segments))]
weights1a = [radii1[i]*ana1.getSegmentLength(i) for i in range(0, len(ana1.segments))]
cts1 = np.array(ana1.data["creationTime"])
l1_, t1_ = np.histogram(cts1, 100, weights = weights1)
r1_, t1_ = np.histogram(cts1, 100, weights = weights1a)

radii2 = ana2.data["radius"] 
weights2 = [ana2.getSegmentLength(i) for i in range(0, len(ana2.segments))]
weights2a = [radii2[i]*ana2.getSegmentLength(i) for i in range(0, len(ana2.segments))]
cts2 = np.array(ana2.data["creationTime"])
l2_, t2_ = np.histogram(cts2, 100, weights = weights2)
r2_, t2_ = np.histogram(cts2, 100, weights = weights2a)

mean_rad1 = np.cumsum(r1_)/np.cumsum(l1_)
mean_rad2 = np.cumsum(r2_)/np.cumsum(l2_)
l1 = np.cumsum(l1_)
l2 = np.cumsum(l2_)
misslen = l1-l2
recovery_rate = l2/l1*100
rad_virt = ((mean_rad1*l1)-(mean_rad2*l2))/(l1-l2)
rad_virt[rad_virt<0] = np.nan
rad_virt[rad_virt>1] = np.nan

print('The recovery rate on day 10 is ', recovery_rate[-1], '%')

print('The missing root length on day 10 is ', misslen[-1], 'cm and should have a mean radius of ', rad_virt[-1], 'cm')

print('is the radius correct?', (misslen[-1]*rad_virt[-1]+l2[-1]*mean_rad2[-1])/(misslen[-1]+l2[-1]), mean_rad1[-1])


fig, ax1 = plt.subplots()
ax1.plot(np.nan, np.nan, 'k-', label = 'Length')
ax1.plot(np.nan, np.nan, 'k--', label = 'Radius')
ax1.plot(0.5 * (t1_[1:] + t1_[:-1]), l1, label = 'Original')
ax1.plot(0.5 * (t1_[1:] + t1_[:-1]), l2, label = 'Reconstructed')
ax1.set_ylabel('Total root length (cm)')
ax2 = ax1.twinx()
ax2.plot(0.5 * (t1_[1:] + t1_[:-1]), mean_rad1, '--')
ax2.plot(0.5 * (t1_[1:] + t1_[:-1]), mean_rad2, '--')
ax2.set_ylim([0, 0.2])
ax2.set_ylabel('Mean root radius (cm)')
ax1.set_xlabel('Time (d)')
ax1.legend()
plt.show()
 
