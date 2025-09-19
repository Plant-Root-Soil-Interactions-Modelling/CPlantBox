import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 11:05:04 2024

@author: srao
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
##############################################################################
# from agrocpb_core import AgroCPB
# agrocpb = AgroCPB(__file__)
##############################################################################
import plantbox as pb

# Path to your maize XML parameter file
path = "../../../../CPlantBox/modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"  # Assuming the XML file is named maize.xml

# Initialize the root system
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.initialize(False)

# Simulation parameters
depth = 150  # Max rooting depth in cm (-1480 mm)
layers = 15  # Number of soil layers
simulation_time = 30  # Simulate for 15 days

rrd = []  # Initialize as an empty list for relative root density
ald = []  # Initialize as an empty list for absolute root length density

# Simulation loop
for t in range(simulation_time):
    rs.simulate(1, False)

    # Calculate root density
    ana = pb.SegmentAnalyser(rs)
    root_length = ana.distribution("length", 0., -depth, layers, True)
    root_length = np.array(root_length)

    soil_volume = (depth / layers) * 10 * 10  # Volume of each soil layer (cm^3)
    root_length_density = root_length / soil_volume  # Convert to density (cm/cm^3)

    # Calculate relative root density
    total_root_length = np.sum(root_length)
    relative_root_density = root_length / total_root_length

    # Ensure relative root densities sum to 1 (allowing for small floating-point errors)
    assert np.isclose(np.sum(relative_root_density), 1.0), f"Relative root densities at day {t} do not sum to 1"

    rrd.append(relative_root_density)
    ald.append(root_length_density)

# Convert lists to numpy arrays for easier plotting
rrd = np.array(rrd)
ald = np.array(ald)

# Plotting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (20, 8))

# Create a mesh for the heatmap
x = np.arange(simulation_time)
y = np.linspace(0, depth, layers)
X, Y = np.meshgrid(x, y)

# Plot absolute root length density
c1 = ax1.pcolormesh(X, Y, ald.T, cmap = 'YlOrBr', shading = 'auto')
ax1.set_title('Absolute Root Length Density Over Time', fontsize = 16)
ax1.set_xlabel('Time (days)', fontsize = 12)
ax1.set_ylabel('Depth (cm)', fontsize = 12)
ax1.invert_yaxis()  # Invert y-axis to have 0 at the top
cbar1 = fig.colorbar(c1, ax = ax1)
cbar1.set_label('Root Length Density (cm/cm³)', fontsize = 12)

# Plot relative root density
c2 = ax2.pcolormesh(X, Y, rrd.T, cmap = 'YlOrBr', shading = 'auto')
ax2.set_title('Relative Root Density Over Time', fontsize = 16)
ax2.set_xlabel('Time (days)', fontsize = 12)
ax2.set_ylabel('Depth (cm)', fontsize = 12)
ax2.invert_yaxis()  # Invert y-axis to have 0 at the top
cbar2 = fig.colorbar(c2, ax = ax2)
cbar2.set_label('Relative Root Density', fontsize = 12)

plt.tight_layout()
plt.show()

# Optional: Print some statistics
print(f"Maximum absolute root length density: {np.max(ald):.4f} cm/cm³")
print(f"Maximum relative root density: {np.max(rrd):.4f}")
print(f"Total root length at end of simulation: {np.sum(root_length):.2f} cm")
