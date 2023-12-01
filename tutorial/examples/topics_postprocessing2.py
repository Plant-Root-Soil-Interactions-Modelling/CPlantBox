"""analysis of results using signed distance functions to crop the segments to a certain geometry"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_1_Leitner_2010"

plant = pb.Plant()
plant.readParameters(path + name + ".xml")
# plant.setGeometry(pb.SDF_PlantContainer(10, 10, 200, True))
plant.initialize()
plant.simulate(120)
plant.write("results/topics_postprocessing2.vtp")

ana = pb.SegmentAnalyser(plant)
print(ana.getMinBounds(), ana.getMaxBounds())

r, depth, layers = 5, 100., 100  # Soil core analysis
soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
soilcolumn2 = pb.SDF_RotateTranslate(soilcolumn, 0, 0, pb.Vector3d(10, 0, 0))  # shift for 10 cm

geom = soilcolumn  # pick one geometry for further analysis

z_ = np.linspace(0, -1 * depth, layers)
fig, axes = plt.subplots(nrows = 1, ncols = 3, figsize = (16, 8))
for a in axes:
    a.set_xlabel('RLD (cm/cm^3)')  # layer size is 1 cm
    a.set_ylabel('Depth (cm)')

# Make a root length distribution
layerVolume = depth / layers * 20 * 20
times = [120, 60, 30, 10]
ana = pb.SegmentAnalyser(plant)
ana.cropDomain(20, 20, depth)  # ana.mapPeriodic(20, 20)
axes[0].set_title('All roots in 20*20*100')
rl_ = []
for t in times:
    ana.filter("creationTime", 0, t)
    rl_.append(ana.distribution("length", 0., -depth, layers, True))
    axes[0].plot(np.array(rl_[-1]) / layerVolume, z_, label = "{:g} days".format(t))
axes[0].legend()  # ["10 days", "30 days", "60 days", "120 days"]

# Make a root length distribution along the soil core
layerVolume = depth / layers * r * r * np.pi
ana = pb.SegmentAnalyser(plant)
ana.crop(geom)
ana.pack()
ana_vis = pb.SegmentAnalyser(ana)
axes[1].set_title('Soil core (over time)')
rl_ = []
for t in times:
    ana.filter("creationTime", 0, t)
    rl_.append(ana.distribution("length", 0., -depth, layers, True))
    axes[1].plot(np.array(rl_[-1]) / layerVolume, z_, label = "{:g} days".format(t))
axes[1].legend()

# distributions per root type
ana = pb.SegmentAnalyser(plant)
ana.crop(geom)
ana.pack()
axes[2].set_title('Soil core (per subType)')
rl_ = []
for i in range(1, 5):
    a = pb.SegmentAnalyser(ana)  # copy
    a.filter("subType", i)
    rl_.append(a.distribution("length", 0., -depth, layers, True))

axes[2].plot((np.array(rl_[0]) + np.array(rl_[3])) / layerVolume, z_)
axes[2].plot(np.array(rl_[1]) / layerVolume, z_)
axes[2].plot(np.array(rl_[2]) / layerVolume, z_)
axes[2].legend(["basal roots", "first order roots", "second order roots"])

fig.subplots_adjust()
plt.savefig("results/topics_postprocessing2.png")
plt.show()

vp.plot_roots(ana_vis, "subType")
