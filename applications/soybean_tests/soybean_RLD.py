"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../..")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"

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
rs.simulate(154, True)

depth, layers = 150, 30
z_ = np.linspace(0, -1 * depth, layers)
#fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 8))
img = plt.imread("soybean_2020.png")
fig, axes = plt.subplots(figsize = (6, 8))
axes.imshow(img, extent=[0, 37, -100, 0])
# Make a root length distribution
ana = pb.SegmentAnalyser(rs)
ana.mapPeriodic(37, 6)
layerVolume = depth / layers * 37 * 6  # actually the only thing that changes
ana.write("results/" + name + "/" + name + "_periodic_154days.vtp")
rl0_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 120)
rl1_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 60)
rl2_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 40)
rl3_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 20)
rl4_ = ana.distribution("length", 0., -depth, layers, True)
axes.set_title('All roots (periodic)')
axes.plot(np.array(rl4_) / layerVolume * 10, z_)
axes.plot(np.array(rl3_) / layerVolume * 10, z_)
axes.plot(np.array(rl2_) / layerVolume * 10, z_)
axes.plot(np.array(rl1_) / layerVolume * 10, z_)
axes.plot(np.array(rl0_) / layerVolume * 10, z_, color = 'goldenrod')
axes.set_xlabel('$\dfrac{1}{10}$ RLD (cm/cm^3)')
axes.set_ylabel('Depth (cm)')
axes.legend(["20 days", "40 days", "60 days", "120 days", "154 days"], loc = 'lower right')
#axes.set_xlim(0,2)
axes.set_ylim(-150,0)
fig.subplots_adjust()
plt.savefig("results/" + name + "/" + name + "_RLDperiodicSoil.pdf", dpi = 300)

fig, ax = plt.subplots()
ax.plot(np.array(rl4_) / layerVolume, z_)
ax.plot(np.array(rl3_) / layerVolume, z_)
ax.plot(np.array(rl2_) / layerVolume, z_)
ax.plot(np.array(rl1_) / layerVolume, z_)
ax.plot(np.array(rl0_) / layerVolume, z_)
ax.set_xlabel('RLD (cm/cm^3)')
ax.set_ylabel('Depth (cm)')
ax.legend(["20 days", "40 days", "60 days", "120 days", "154 days"], loc = 'lower right')
ax.set_ylim(-150,0)
fig.subplots_adjust()
plt.savefig("results/" + name + "/" + name + "_RLDperiodic.pdf", dpi = 300)

plt.show()   
