"""root system surface density"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np
import matplotlib.pyplot as plt

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"
plant.readParameters(path + name + ".xml")

depth = 220
layers = 50
runs = 10

box = pb.SDF_PlantContainer(5, 5, depth, True)  # [cm3]

rs_ = []
for i in range(0, runs):
    plant.setGeometry(box)
    plant.initialize(False)
    plant.simulate(90, False)
    ana = pb.SegmentAnalyser(plant)
    rs_.append(ana.distribution("surface", 0., -depth, layers, True))

soilvolume = (depth / layers) * 10 * 10  # [cm3]
rs_ = np.array(rs_) / soilvolume  # convert to density [cm2/cm3]
rs_mean = np.mean(rs_, axis = 0)
rs_err = np.std(rs_, axis = 0)

dx2 = 0.5 * (depth / layers)  # half layer width
z_ = np.linspace(-dx2, -depth + dx2, layers)  # layer mid points
plt.plot(rs_mean, z_, "b*")
plt.plot(rs_mean + rs_err, z_, "b:")
plt.plot(rs_mean - rs_err, z_, "b:")
plt.xlabel("root surface (cm^2 / cm^3)")
plt.ylabel("z-coordinate (cm)")
plt.legend(["mean value (" + str(runs) + " runs)", "std"])
plt.savefig("results/topics_postprocessing.png")
plt.show()

print(ana.getMinBounds(), ana.getMaxBounds())
vp.plot_roots(ana, "subType")
