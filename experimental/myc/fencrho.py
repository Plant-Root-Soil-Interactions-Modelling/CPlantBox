import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import time
import matplotlib as mpl
import AMFAnalysis as amf

# I changed the c++ code, so set might not be valid anymore
good_seeds = [10, 15, 25, 30, 45, 50, 70, 105, 115, 125]

path = "tomatoparameters/"
name = "TwoHyphaePlusBAS"

## Setting up petri dish
diameter = 9.4
radius = diameter /2
height = 1.6
    # introduce parameters for barrier and opening
barrier_thickness = 0.16
barrier_height = height 
opening_length = 5.0
opening_height = 0.2

nRings = 15

petri_dish, small_hyphae_dish, half_dish, rings = amf.makedishes(diameter, height, barrier_thickness, barrier_height, opening_height, opening_length, 0.1, nRings)

dt, nSteps = amf.setUpSimulationTime(5,24)
seed = good_seeds[0]
mycp, crossed_time, times = amf.makeSimulation(seed, path, name, height,petri_dish, small_hyphae_dish, half_dish,dt,nSteps, 20, "f_enc_vs_rho_h",animation=False)
ana = pb.SegmentAnalyser()
ana = amf.getMycSegmentAnalyser(mycp)
tipMat = amf.getParaSumperRing("nodeTips",times,ana,rings)
print(tipMat)