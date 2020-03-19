"""analysis of results using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../../..")
import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../../modelparameter/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(120)
rs.write("results/example_3d.vtp")

# Soil core analysis
r, depth, layers = 10, 100., 100

soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
soilcolumn2 = pb.SDF_RotateTranslate(soilcolumn, 0, 0, pb.Vector3d(12, 0, 0))  # shift 10 cm

# pick one geometry for further analysis
geom = soilcolumn

z_ = np.linspace(0, -1 * depth, layers)
fig, axes = plt.subplots(nrows = 1, ncols = 4, figsize = (16, 8))
for a in axes:
    a.set_xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
    a.set_ylabel('Depth (cm)')

# Make a root length distribution
ana = pb.SegmentAnalyser(rs)
ana.cropDomain(20, 20, depth)
layerVolume = depth / layers * 20 * 20
rl0_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 60)
rl1_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 30)
rl2_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 10)
rl3_ = ana.distribution("length", 0., -depth, layers, True)
axes[0].set_title('All roots in 20*20*100')
axes[0].plot(np.array(rl3_) / layerVolume, z_)
axes[0].plot(np.array(rl2_) / layerVolume, z_)
axes[0].plot(np.array(rl1_) / layerVolume, z_)
axes[0].plot(np.array(rl0_) / layerVolume, z_)
axes[0].legend(["10 days", "30 days", "60 days", "120 days"])

# Make a root length distribution along the soil core
layerVolume = depth / layers * r * r * pi
ana = pb.SegmentAnalyser(rs)
ana.crop(geom)
ana.pack()
rl0_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 60)
rl1_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 30)
rl2_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 10)
rl3_ = ana.distribution("length", 0., -depth, layers, True)
axes[1].set_title('Soil core')
axes[1].plot(np.array(rl3_) / layerVolume, z_)
axes[1].plot(np.array(rl2_) / layerVolume, z_)
axes[1].plot(np.array(rl1_) / layerVolume, z_)
axes[1].plot(np.array(rl0_) / layerVolume, z_)
axes[1].legend(["10 days", "30 days", "60 days", "120 days"])

# Only laterals
ana = pb.SegmentAnalyser(rs)
ana.cropDomain(20, 20, depth)
layerVolume = depth / layers * 20 * 20
ana.pack()
a = pb.SegmentAnalyser(ana)  # copy
a.filter("subType", 4)  # basal
rl4_ = a.distribution("length", 0., -depth, layers, True)
a = pb.SegmentAnalyser(ana)  # copy
a.filter("subType", 3)  # 2nd order
rl3_ = a.distribution("length", 0., -depth, layers, True)
a = pb.SegmentAnalyser(ana)  # copy
a.filter("subType", 2)  # 1st order
rl2_ = a.distribution("length", 0., -depth, layers, True)
a = pb.SegmentAnalyser(ana)  # copy
a.filter("subType", 1)  # tap
rl1_ = a.distribution("length", 0., -depth, layers, True)
axes[2].set_title('Soil core')
axes[2].plot((np.array(rl1_) + np.array(rl4_)) / layerVolume, z_)
axes[2].plot(np.array(rl2_) / layerVolume, z_)
axes[2].plot(np.array(rl3_) / layerVolume, z_)
axes[2].legend(["basal roots", "first order roots", "second order roots"])

# Make a root length distribution
ana = pb.SegmentAnalyser(rs)
ana.mapPeriodic(20, 20)
layerVolume = depth / layers * 20 * 20  # actually the only thing that changes
ana.write("results/example_3d_periodic.vtp")
rl0_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 60)
rl1_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 30)
rl2_ = ana.distribution("length", 0., -depth, layers, True)
ana.filter("creationTime", 0, 10)
rl3_ = ana.distribution("length", 0., -depth, layers, True)
axes[3].set_title('All roots (periodic)')
axes[3].plot(np.array(rl3_) / layerVolume, z_)
axes[3].plot(np.array(rl2_) / layerVolume, z_)
axes[3].plot(np.array(rl1_) / layerVolume, z_)
axes[3].plot(np.array(rl0_) / layerVolume, z_)
axes[3].legend(["10 days", "30 days", "60 days", "120 days"])

fig.subplots_adjust()
plt.savefig("results/example_3d.png")
plt.show()
