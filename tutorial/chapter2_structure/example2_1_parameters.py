"""simple root system from scratch (without parameter files)"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import matplotlib.pyplot as plt
import numpy as np

plant = pb.Plant()
p0 = pb.RootRandomParameter(plant)  # with default values,
p1 = pb.RootRandomParameter(plant)  # all standard deviations are 0

p0.name = "taproot"
p0.a = 0.2  # [cm] radius
p0.subType = 1  # [-] index starts at 1
p0.lb = 5  # [cm] basal zone
p0.la = 10  # [cm] apical zone
p0.lmax = 30  # [cm] maximal root length, number of lateral branching nodes = round((lmax-lb-la)/ln) + 1
p0.ln = 1.  # [cm] inter-lateral distance (16 branching nodes)
p0.theta = 0.  # [rad]
p0.r = 1  # [cm/day] initial growth rate
p0.dx = 10  # [cm] axial resolution
p0.successor = [[2]]  # add successors
p0.successorP = [[1]]  # probability that successor emerges
p0.tropismT = pb.TropismType.gravi  #
p0.tropismN = 1.8  # [-] strength of tropism
p0.tropismS = 0.2  # [rad/cm] maximal bending

p1.name = "lateral"
p1.a = 0.1  # [cm] radius
p1.subType = 2  # [1] index starts at 1
p1.lmax = 15  # # [cm] apical zone
p1.lmaxs = 0.15  # [cm] standard deviation of the apical zone
p1.theta = 90. / 180.*np.pi  # [rad]
p1.r = 2  # initial growth rate
p1.dx = 1  # [cm] axial resolution
p1.tropismT = pb.TropismType.gravi  # exo
p1.tropismN = 2  # [-] strength of tropism
p1.tropismS = 0.1  # [rad/cm] maximal bending

plant.setOrganRandomParameter(p0)
plant.setOrganRandomParameter(p1)

srp = pb.SeedRandomParameter(plant)  # with default values
srp.seedPos = pb.Vector3d(0., 0., -3.)  # [cm] seed position
srp.maxB = 0  # [-] number of basal roots (neglecting basal roots and shoot borne)
srp.firstB = 10.  # [day] first emergence of a basal root
srp.delayB = 3.  # [day] delay between the emergence of basal roots
plant.setOrganRandomParameter(srp)

plant.initialize()

fig, axes = plt.subplots(1, 3, figsize = (15, 7))
simtimes = [0, 30, 60, 125]  # the last lateral will emerge at
for i in range(0, 3):

    plant.simulate(np.diff(simtimes)[i])  #  [day]
    a = axes[i]
    a.set_xlim([-15, 15.])
    a.set_ylim([-35., 0.])  # starts at -3 cm, max lengthsimple root system from scratch (without parameter files) 30 cm
    a.set_title("after {} days".format(plant.getSimTime()))
    roots = plant.getPolylines()
    for root in roots:
        for j, n in enumerate(root[:-1]):
            n2 = root[j + 1]
            a.plot([n.x, n2.x], [n.z, n2.z], "g")

fig.tight_layout()
plt.show()

plant.write("/results/topics_parameters2.vtp")

# Some outputs....
print(" length", plant.getParameter("length"))
print("    age", plant.getParameter("age"))
print("subType", plant.getParameter("subType"))
print("     la", plant.getParameter("la"))
print("la_mean", plant.getParameter("la_mean"))
print(" radius", plant.getParameter("radius"))
