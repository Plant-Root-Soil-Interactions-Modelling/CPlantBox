"""simple root system from scratch (without parameter files)"""

import matplotlib.pyplot as plt  # |\label{l2_1:matplotlib}|
import numpy as np  # |\label{l2_1:numpy}|

import plantbox as pb
from plantbox.visualisation import figure_style
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
p0 = pb.RootRandomParameter(plant)  # |\label{l2_1:p0}|
p1 = pb.RootRandomParameter(plant)  # |\label{l2_1:p1}|
s1 = pb.StemRandomParameter(plant)  # |\label{l2_1:s1}|
l1 = pb.LeafRandomParameter(plant)  # |\label{l2_1:l1}|

p0.name = "taproot"  # |\label{l2_1:tap_start}|
p0.a = 0.2  # radius (cm)
p0.subType = 1  # index starts at 1
p0.lb = 5  # basal zone (cm)
p0.la = 10  # apical zone (cm)
p0.lmax = 30  # maximal root length (cm)
p0.ln = 1.0  # inter-lateral distance (cm)
p0.theta = 0.0  # (rad)
p0.r = 1  # initial growth rate (cm day-1)
p0.dx = 10  # axial resolution (cm)
p0.successorST = [[2]]  # add successors
p0.successorP = [[1]]  # probability that successor emerges
p0.tropismT = pb.TropismType.gravi  #
p0.tropismN = 1.8  # strength of tropism (1)
p0.tropismS = 0.2  # maximal bending (rad/cm) |\label{l2_1:tap_end}|

p1.name = "lateral"  # |\label{l2_1:lat_start}|
p1.a = 0.1  # radius (cm)
p1.subType = 2  # index starts at 1
p1.lmax = 15  # maximal root length (cm)
p1.lmaxs = 0.15  # standard deviation of the maximal root length (cm) |\label{l2_1:lat_lmaxs}|
p1.theta = 90.0 / 180.0 * np.pi  # (rad)
p1.r = 2  # initial growth rate (cm day-1)
p1.dx = 1  # axial resolution (cm)
p1.tropismT = pb.TropismType.gravi  #
p1.tropismN = 2  # strength of tropism (1)
p1.tropismS = 0.1  # maximal bending (rad cm-1) |\label{l2_1:lat_end}|

s1.name = "stem"  # |\label{l2_1:stem_start}|
s1.subType = 1  # radius (cm)
s1.a = 0.2  # radius (cm)
s1.ldelay = -1  # delay between lateral creation and start of growth
s1.nodalGrowth = 0  # inter-lateral distance (cm)
s1.lb = 5  # basal zone (cm)
s1.la = 10  # apical zone (cm)
s1.lmax = 30  # maximal root length (cm)
s1.ln = 1.0  # inter-lateral distance (cm)
s1.r = 2  # growth rate (cm)
s1.successorOT = [[4]]
s1.successorST = [[1]]
s1.successorP = [[1]]  # |\label{l2_1:stem_end}|

l1.name = "leaf"  # |\label{l2_1:leaf_start}|
l1.subType = 1  # radius (cm)
l1.ldelay = -1  # delay between lateral creation and start of growth
l1.r = 5  # growth rate (cm)
l1.a = 0.05
l1.la, l1.lb, l1.lmax, l1.ln = 3.5, 1.0, 7.5, 3
l1.shapeType = 2
l1.areaMax = 10
l1.leafGeometryPhi = np.array([-90.0, -45, 0.0, 45.0, 90.0]) / 180.0 * np.pi
l1.leafGeometryX = np.array([3.0, 2.2, 1.7, 2.0, 3.5])
N = 101
l1.createLeafRadialGeometry(l1.leafGeometryPhi, l1.leafGeometryX, N)  # |\label{l2_1:leaf_end}|

plant.setOrganRandomParameter(p0)  # |\label{l2_1:set_p0}|
plant.setOrganRandomParameter(p1)  # |\label{l2_1:set_p1}|
plant.setOrganRandomParameter(l1)  # |\label{l2_1:set_L1}|
plant.setOrganRandomParameter(s1)  # |\label{l2_1:set_S1}|

srp = pb.SeedRandomParameter(plant)  # with default values |\label{l2_1:srp_start}|
srp.seedPos = pb.Vector3d(0.0, 0.0, -3.0)  # seed position (cm)
plant.setOrganRandomParameter(srp)  # |\label{l2_1:srp_end}|

plant.initialize(False)

plant.simulate(50, False)  # |\label{l2_1:simulation_start}|
vp.plot_plant(plant, "creationTime")
plant.write("results/example2_1_parameters.vtp")  # |\label{l2_1:simulation_end}|

# Some outputs....
print(srp)  # |\label{l2_1:print_srp}|
print(s1)  # |\label{l2_1:print_S1}|
