"""root system surface density"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import numpy as np

import matplotlib.pyplot as plt

path = "../../modelparameter/structural/rootsystem/"
name = "Brassica_napus_a_Leitner_2010"  # "Crypsis_aculeata_Clausnitzer_1994"

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

depth = 220
layers = 50
runs = 10

rl_ = []
for i in range(0, runs):
    rs.initialize(False)
    rs.simulate(120, False)
    ana = pb.SegmentAnalyser(rs)
    rl_.append(ana.distribution("length", 0., -depth, layers, True))

soilvolume = (depth / layers) * 10 * 10
rl_ = np.array(rl_) / soilvolume  # convert to density
rl_mean = np.mean(rl_, axis = 0)
rl_err = np.std(rl_, axis = 0) / np.sqrt(runs)

z_ = np.linspace(0, -depth, layers)  # z - axis
plt.plot(rl_mean, z_, "b")
plt.plot(rl_mean + rl_err, z_, "b:")
plt.plot(rl_mean - rl_err, z_, "b:")

plt.xlabel("root surface (cm^2 / cm^3)")
plt.ylabel("z-coordinate (cm)")
plt.legend(["mean value (" + str(runs) + " runs)", "error"])
plt.savefig("results/example_3a.png")
plt.show()
