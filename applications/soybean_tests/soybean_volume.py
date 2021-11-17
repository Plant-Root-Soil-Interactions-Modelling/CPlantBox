"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../..")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt
import csv
from itertools import zip_longest
import math

path = "../../modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"  # ""

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Create and set geometry
rs.setMinDx(1.e-3)
x0 = pb.Vector3d(0., 0., -1.)
nx = pb.Vector3d(1., 0., -1.)
ny = pb.Vector3d(0., 1., -1.)
soil_layer = pb.SDF_HalfPlane(x0, nx, ny)  # there was bug, with updated CPlantBox
rs.setGeometry(soil_layer)

rs.setSeed(2)
rs.initialize()

simtime = 154.  # days
dt = 1.
N = round(simtime / dt)  # steps
v_ = np.zeros(N)
t_ = np.linspace(dt, N * dt, N)

for i in range(0, N):
    rs.simulate(dt)
    l = rs.getSummed("length")
    v = rs.getSummed("volume")
    v_[i] = np.sum(v)
    
print(v_)

volume_ = [t_, v_]
export_data = zip_longest(*volume_, fillvalue = '')
with open("results/" + name + "/" + name + "_volume.csv", "w") as f:
	writer = csv.writer(f)
	writer.writerow(("x", "y"))
	writer.writerows(export_data)
f.close()

fig, ax1 = plt.subplots()
ax1.plot(t_, v_)
ax1.set_xlabel("Time [days]")
ax1.set_ylabel(" Rootsystem volume [cm^3]")
ax1.legend(["total"])
plt.savefig("results/" + name + "/" + name + "_volume.pdf")

plt.show()
