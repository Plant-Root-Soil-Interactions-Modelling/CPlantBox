import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import matplotlib.pyplot as plt

import plantbox as pb

# Path to your maize XML parameter file
path = "../../modelparameter/structural/rootsystem/"
name = "Glycine_max"

# Initialize the root system
mycp = pb.MycorrhizalPlant()
mycp.readParameters(path + name + ".xml")
mycp.initialize(True)

# Simulation parameters
simtime = 10
fpd = 30
anim_time = simtime
N = fpd * anim_time
dt = 1 / fpd

# Set up soil geometry
depth, layers = 150, 30
z_ = np.linspace(0, -1 * depth, layers)

# Create and set hyphal parameters
hyphae_parameter = pb.HyphaeRandomParameter(mycp)
hyphae_parameter.subType = 1
hyphae_parameter.a = 0.01
hyphae_parameter.dx = 0.01
mycp.setOrganRandomParameter(hyphae_parameter)

# Make a hyphal length distribution

# rhd = []  # Initialize as an empty list for relative hyphae density
# ald = []  # Initialize as an empty list for absolute hyphae length density
# Simulation loop
for t in range(1,N):
    mycp.simulate(dt, True)

 # Calculate hyphal density
ana = pb.SegmentAnalyser(mycp)
ana.filter("organType",5)
# hyphae_length = ana.distribution("length", 0., -depth, layers, True)
# hyphae_length = np.array(hyphae_length)

soil_volume = (depth / layers) * 10 * 10  # Volume of each soil layer (cm^3)
# hyphae_length_density = hyphae_length / soil_volume  # Convert to density (cm/cm^3)

#     # Calculate relative hyphal density
#     # total_hyphae_length = np.sum(hyphae_length)
#     # relative_hyphae_density = hyphae_length / total_hyphae_length

#     # Ensure relative hyphal densities sum to 1 (allowing for small floating-point errors)
#     # assert np.isclose(np.sum(relative_hyphae_density), 1.0), f"Relative root densities at day {t} do not sum to 1"

#     # rhd.append(relative_hyphae_density)
# ald.append(hyphae_length_density)

# # Convert lists to numpy arrays for easier plotting
# # rhd = np.array(rhd)
# ald = np.array(ald)
# { 
# # Plotting
# # fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (20, 8))

# # # Create a mesh for the heatmap
# # # x = np.arange(simtime)
# # x = np.linspace(0, simtime, N)
# # y = np.linspace(0, depth, layers)
# # X, Y = np.meshgrid(x, y)

# # # Plot absolute hyphae length density
# # c1 = ax1.pcolormesh(X, Y, ald.T, cmap = 'YlOrBr', shading = 'auto')
# # ax1.set_title('Absolute Hyphal Length Density Over Time', fontsize = 16)
# # ax1.set_xlabel('Time (days)', fontsize = 12)
# # ax1.set_ylabel('Depth (cm)', fontsize = 12)
# # ax1.invert_yaxis()  # Invert y-axis to have 0 at the top
# # cbar1 = fig.colorbar(c1, ax = ax1)
# # cbar1.set_label('Hyphal Length Density (cm/cm³)', fontsize = 12)

# # # Plot relative hyphae density
# # # c2 = ax2.pcolormesh(X, Y, rhd.T, cmap = 'YlOrBr', shading = 'auto')
# # # ax2.set_title('Relative Hyphal Density Over Time', fontsize = 16)
# # # ax2.set_xlabel('Time (days)', fontsize = 12)
# # # ax2.set_ylabel('Depth (cm)', fontsize = 12)
# # # ax2.invert_yaxis()  # Invert y-axis to have 0 at the top
# # # cbar2 = fig.colorbar(c2, ax = ax2)
# # # cbar2.set_label('Relative Hyphal Density', fontsize = 12)

# # plt.tight_layout()
# # plt.show()
# }


# Optional: Print some statistics
# print(f"Maximum absolute hyphal length density: {np.max(ald):.4f} cm/cm³")
# print(f"Total hyphal length at end of simulation: {np.sum(hyphae_length):.2f} cm")

# layerVolume = depth / layers * 37 * 6  # actually the only thing that changes
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 8))
# ana.write("results/" + name + "/" + name + "_periodic_154days.vtp")
rl0_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, simtime)
rl1_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, simtime*0.75)
rl2_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, simtime*0.5)
rl3_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, simtime*0.25 )
rl4_ = ana.distribution("length", 0., -depth, layers, True)
axes.set_title('All hyphae')
axes.plot(np.array(rl4_) / soil_volume * 10, z_)
axes.plot(np.array(rl3_) / soil_volume * 10, z_)
axes.plot(np.array(rl2_) / soil_volume * 10, z_)
axes.plot(np.array(rl1_) / soil_volume * 10, z_)
axes.plot(np.array(rl0_) / soil_volume * 10, z_, color = 'goldenrod')
axes.set_xlabel('$\dfrac{1}{10}$ HLD (cm/cm^3)')
axes.set_ylabel('Depth (cm)')
axes.legend([f"{simtime*t:.0f} days" for t in [0.25, 0.5, 0.75, 1.0]], loc = 'lower right')
axes.set_xlim(0,2)
axes.set_ylim(-depth,0)
fig.subplots_adjust()
# plt.savefig("results/" + name + "/" + name + "_RLDperiodicSoil.pdf", dpi = 300)

fig, ax = plt.subplots()
ax.plot(np.array(rl4_) / soil_volume, z_)
ax.plot(np.array(rl3_) / soil_volume, z_)
ax.plot(np.array(rl2_) / soil_volume, z_)
ax.plot(np.array(rl1_) / soil_volume, z_)
ax.plot(np.array(rl0_) / soil_volume, z_)
ax.set_xlabel('HLD (cm/cm^3)')
ax.set_ylabel('Depth (cm)')
ax.legend([f"{simtime*t:.0f} days" for t in [0.25, 0.5, 0.75, 1.0]], loc = 'lower right')
ax.set_ylim(-depth,0)
fig.subplots_adjust()
# plt.savefig("results/" + name + "/" + name + "_RLDperiodic.pdf", dpi = 300)

plt.show()   
