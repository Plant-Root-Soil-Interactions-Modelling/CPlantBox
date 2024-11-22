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

n = [[8, 8, 25],[16, 16, 30]]
name = ['sims_root_soil_coarse', 'sims_root_soil_fine']
RSage = 11

for i in range(0,len(name)): 
    nx = n[i][0]
    ny =n[i][1]
    nz = n[i][2]

    name_ = name[i]
    npzfile = np.load(name_ + '.npz')
    swc =npzfile['arr_0']
    rwc =npzfile['arr_1']
    cn = npzfile['arr_2']

    namep = name_.split('_')

    #visualise soil water content
    swc = swc*int(255) #int(8)
    #swc = np.swapaxes(swc, 0, 2)
    swc.astype('int8').tofile('results/swc_' + namep[3] + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')



    #visualise root water content
    rwc = rwc*int(255) #int(8)
    #rwc = np.swapaxes(rwc, 0, 2)
    rwc.astype('int8').tofile('results/rwc_' + namep[3] + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')


    #visualise soil + root water content
    srwc = swc + rwc
    #srwc = srwc*int(255) #int(8)
    #srwc = np.swapaxes(srwc, 0, 2)
    srwc.astype('int8').tofile('results/srwc_' + namep[3] + str(nx) + 'x' + str(ny) + 'x' + str(nz) + '.raw')


