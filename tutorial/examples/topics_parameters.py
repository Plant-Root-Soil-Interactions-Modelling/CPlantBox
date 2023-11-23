""" something basic"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import numpy as np

plant = pb.Plant()

stem = pb.StemRandomParameter(plant)
stem.subType = 1

leaf = pb.LeafRandomParameter(plant)

seed = pb.SeedRandomParameter(plant)
# seed.seedPos = pb.Vector3d(0., 0., -3.)  # [cm] seed position
# seed.maxB = 0  # [-] number of basal roots (neglecting basal roots and shoot borne)
# seed.firstB = 10.  # [day] first emergence of a basal root
# seed.delayB = 3.  # [day] delay between the emergence of basal roots

plant.setOrganRandomParameter(stem)
# plant.setOrganRandomParameter(leaf)
plant.setOrganRandomParameter(seed)

p0 = pb.RootRandomParameter(plant)
p1 = pb.RootRandomParameter(plant)
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
# p0.successor = [[]]  # add successors
# p0.successorP = [[1]]  # probability that successor emerges
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
print("successor", p0.successor)
print("successorP", p0.successorP)
plant.setOrganRandomParameter(p1)

plant.initialize()

plant.simulate(10)

vp.plot_plant(plant, "subType")

