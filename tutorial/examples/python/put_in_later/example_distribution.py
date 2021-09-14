"""small example"""
import sys
base = "../../../.."
sys.path.append(base); sys.path.append(base+"/src/python_modules")
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = base+"/modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initialize()

# Simulate
rs.simulate(30, True)

Nz = 10
z_ = np.linspace(0, -50, Nz)+(-50)/Nz/2 # cell mid points
ana = pb.SegmentAnalyser(rs)
l_ = ana.distribution("length", 0, -50, Nz, True) #  distribution(std::string name, double top, double bot, int n, bool exact=false)
plt.plot(l_,z_)
plt.show()

Nz = 20
Nx = 10
ana = pb.SegmentAnalyser(rs)
l_ = ana.distribution2("surface", 0,-50, -10, 10, Nz, Nx, True) # distribution2(std::string name, double top, double bot, double left, double right, int n, int m, bool exact=false) const; ///< 2d distribution (x,z) of a parameter
l_ = np.array([np.array(l) for l in l_]) # convert list of list to a numpy array
print(l_.shape)
plt.imshow(l_)
plt.colorbar()
plt.show()

# Export final result (as vtp)
#rs.write("results/example_distribution.vtp")

# # Plot, using vtk
# vp.plot_roots(rs, "creationTime")