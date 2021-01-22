"""Spring barley Multiple Rootsystems Simulation PW"""
import sys; sys.path.append("../.."); sys.path.append("../../src/python_modules")
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import vtk_plot as vp
import plantbox as pb
import math


path = "../../modelparameter/rootsystem/"
# name = "spring_barley_Fernand"
name =  "wheat"

simtime = 107
N = 8                      # number of rows
M = 2                      # number of plants in a row 
distp = 3.                 # distance between the root systems along row[cm]
distr = 12.5               # distance between the rows[cm]
distTr = N*distr           # total row spacing
distTp = M*distp           # total distance between plants 
top = 0.                   # vertical top position (cm) (normally = 0)
bot = -100.                # vertical bot position (cm) (e.g. = -100 cm)
left = -47.5               # left along x-axis (cm)
right = 47.5               # right along x-axis (cm)
n = 20                     # number of layers, each with a height of (top-bot)/n
m = 20                     # number of horizontal grid elements (each with length of (right-left)/m)
exact = True               # calculates the intersection with the layer boundaries (true), only based on segment midpoints (false)


soilVolume =   5*5*5

# Time points during vegetation period in days
times = [73, 72, 71, 46, 45]    # (shooting=45,46) and (flowering= 71,72,73) 


# Make a root length distribution along the soil profile wall 

# Initializes N*M root systems
allRS = []
for i in range(0, N):
    for j in range(0, M):
        rs = pb.RootSystem()
        rs.readParameters(path + name + ".xml")
        rs.getRootSystemParameter().seedPos = pb.Vector3d(distr * i, distp * j, -3.)  # cm
        rs.initialize(False)
        allRS.append(rs)
        
# Simulate
time = 0
dt = 1 # day
while time < simtime: # for future coupling with dynamic water movement 
    print("day", time)
    for rs in allRS:
        rs.simulate(dt)
    time += dt
        
# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()                           # see example 3b
for z, rs in enumerate(allRS):
    # vtpname = "results/plantsb" + str(z) + ".vtp"
    # rs.write(vtpname)
    ana.addSegments(rs)  # collect all

# Write all into single file (segments)
ana.write("results/plantsb_all.vtp")

# Set periodic domain
ana.mapPeriodic(distTr, distTp)                         
ana.write("results/plantsb_periodic.vtp")
# vp.plot_roots(ana, "length", "length (cm)")
ana.pack()           
rl_ = []

for j in range(len(times)):
    ana.filter("creationTime", 0, times[j])
    rl_.append(ana.distribution2("length", top, bot, -50, 50, n, m, True))
    rld_ =  np.array(rl_)/soilVolume
print(rld_.shape)

# df = rld_.reshape(-1, rld_.shape[1])

# np.savetxt("sim_data20x20.txt", np.around(df, decimals = 4),fmt='%.9f',delimiter=',')

# Flowering

rld = rld_[0]     						# day 73 during vegetation period                                      
data_flowering3 = np.array([np.array(l) for l in rld])     

rld = rld_[1]							# day 72 during vegetation period
data_flowering2 = np.array([np.array(l) for l in rld])     

rld = rld_[2]							# day 71 during vegetation period
data_flowering1 = np.array([np.array(l) for l in rld])     
# print(data_flowering1.shape)



# Shooting

rld = rld_[3]                                            	# day 46 during vegetation period
data_shooting2 = np.array([np.array(l) for l in rld])     
# print(data_shooting2.shape)

rld = rld_[4]                                             	# day 45 during vegetation period
data_shooting1 = np.array([np.array(l) for l in rld])     
# print(data_shooting1.shape)



x = np.linspace(-47.5, 47.5, 20)
y = np.linspace(-2.5, -97.5, 20)
levels = np.linspace(0, 4, 21)
X, Y = np.meshgrid(x, y)



# Plot during shooting
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
fig.suptitle("RLD during plant shooting[45-46]", fontsize=16)

fig12 = ax1.contourf(X, Y, data_shooting1, cmap='Greens', levels=levels)
ax1.set_xlabel("Distance of the row [cm]")
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position("top")
ax1.set_ylabel("Depth [cm]")
fig.colorbar(fig12, ax=ax1)

fig13 = ax2.contourf(X, Y, data_shooting2, cmap='Reds', levels=levels)
ax2.set_xlabel("Distance of the row [cm]")
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position("top")
fig.colorbar(fig13, ax=ax2)

# plt.savefig("RLD during shooting[45, 46].png")
plt.show()

# Plot during flowering
fig2, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))
fig2.suptitle("RLD during plant flowering[71-73]", fontsize=16)

fig21 = ax1.contourf(X, Y, data_flowering1, cmap='Blues', levels=levels)
ax1.set_xlabel("Distance of the row [cm]")
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position("top")
ax1.set_ylabel("Depth [cm]")
fig2.colorbar(fig21, ax=ax1)

fig22 = ax2.contourf(X, Y, data_flowering2, cmap='Greens', levels=levels)
ax2.set_xlabel("Distance of the row [cm]")
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position("top")
fig2.colorbar(fig22, ax=ax2)

fig23 = ax3.contourf(X, Y, data_flowering3, cmap='Reds', levels=levels)
ax3.set_xlabel("Distance of the row [cm]")
ax3.xaxis.tick_top()
ax3.xaxis.set_label_position("top")
fig2.colorbar(fig23, ax=ax3)

# plt.savefig("RLD during flowering[71, 72, 73].png")
plt.show()

