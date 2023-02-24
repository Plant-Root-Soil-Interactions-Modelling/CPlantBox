"""analysis of results using signed distance functions"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_1_Leitner_2010"  # Zea_mays_1_Leitner_2010, Brassica_napus_a_Leitner_2010

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(120)
rs.write("results/example_3b.vtp")

r, depth, layers = 5, 100., 100  # Soil core analysis
soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
soilcolumn2 = pb.SDF_RotateTranslate(soilcolumn, 0, 0, pb.Vector3d(10, 0, 0))  # shift 10 cm

# pick one geometry for further analysis
geom = soilcolumn

z_ = np.linspace(0, -1 * depth, layers)
fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (16, 8))
for a in axes:
    a.set_xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
    a.set_ylabel('Depth (cm)')

# Make a root length distribution
layerVolume = depth / layers * 20 * 20
times = [120, 60, 30, 10]
ana = pb.SegmentAnalyser(rs)
ana.cropDomain(20, 20, depth)  # ana.mapPeriodic(20, 20)
rl_ = []
axes[0].set_title('All roots in 20*20*100')
for t in times:
    ana.filter("creationTime", 0, t)
    rl_.append(ana.distribution("length", 0., -depth, layers, True))
    axes[0].plot(np.array(rl_[-1]) / layerVolume, z_)
axes[0].legend(["10 days", "30 days", "60 days", "120 days"])

# Make a root length distribution along the soil core
layerVolume = depth / layers * r * r * np.pi
ana = pb.SegmentAnalyser(rs)
ana.crop(geom)
ana.pack()
rl_ = []
axes[1].set_title('Soil core')
for t in times:
    ana.filter("creationTime", 0, t)
    rl_.append(ana.distribution("length", 0., -depth, layers, True))
    axes[1].plot(np.array(rl_[-1]) / layerVolume, z_)
axes[1].legend(["10 days", "30 days", "60 days", "120 days"])

# distributions per root type
ana = pb.SegmentAnalyser(rs)
ana.crop(geom)
ana.pack()
rl_ = []
for i in range(1, 5):
    a = pb.SegmentAnalyser(ana)  # copy
    a.filter("subType", i)
    rl_.append(a.distribution("length", 0., -depth, layers, True))
axes[2].set_title('Soil core')
axes[2].plot((np.array(rl_[0]) + np.array(rl_[3])) / layerVolume, z_)
axes[2].plot(np.array(rl_[1]) / layerVolume, z_)
axes[2].plot(np.array(rl_[2]) / layerVolume, z_)
axes[2].legend(["basal roots", "first order roots", "second order roots"])

fig.subplots_adjust()
plt.savefig("results/example_3b.png")
plt.show()
