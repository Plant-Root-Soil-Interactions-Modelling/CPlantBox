import sys
from cmath import pi
sys.path.append("../../..")
sys.path.append("../../../src/");
sys.path.append("../../../gui/viewer")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython
from viewer_data import ViewerDataModel
from viewer_plots import plot_depth_profile
from visualisation.vtk_tools import *
import rsml.rsml_writer
import rsml.rsml_reader as rsml
import numpy as np
import matplotlib.pyplot as plt

depth = 19
width = 2.75
layers = 8

plt.rcParams.update({'font.size': 18})

#read in the rsmls
path = 'results/'
fname1 = path + 'GT_Faba_day10'
fname2 = path + 'REC_Faba_day10'
fname3 = path + 'semivirt_Faba_day10'
rs1 = XylemFluxPython.read_rsml(fname1 +".rsml")
rs2 = XylemFluxPython.read_rsml(fname2 +".rsml")
rs3 = XylemFluxPython.read_rsml(fname3 +".rsml")

# Make an RLD distribution 
z_ = np.linspace(0, -depth, layers)  # z - axis
soilvolume = np.array((depth / layers) * (width)**2*pi)

ana1 = pb.SegmentAnalyser(rs1)
ana2 = pb.SegmentAnalyser(rs2)
ana3 = pb.SegmentAnalyser(rs3)
rld1 = ana1.distribution("length", 0., -depth, layers, True)/soilvolume
rld2 = ana2.distribution("length", 0., -depth, layers, True)/soilvolume
rld3 = ana3.distribution("length", 0., -depth, layers, True)/soilvolume
plt.plot(rld1, z_, "k-", label = 'Original RS')
plt.plot(rld2, z_, "k--",label = 'Reconstructed RS')
plt.plot(rld3, z_, "k:",label = 'Semivirtual RS')


plt.xlabel("RLD (cm / $cm^3$)")
plt.ylabel("Depth (cm)")
plt.legend()
plt.tight_layout()
plt.show()



