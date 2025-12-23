"""simple root system from scratch (without parameter files)"""

import matplotlib.pyplot as plt  # |\label{l2_1:matplotlib}|
import numpy as np  # |\label{l2_1:numpy}|
import plantbox.visualisation.vtk_plot as vp
import plantbox as pb
from plantbox.visualisation import figure_style

plant = pb.Plant()
p0 = pb.RootRandomParameter(plant)  # with default values  |\label{l2_1:p0}|
p1 = pb.RootRandomParameter(plant)  # all standard deviations are 0  |\label{l2_1:p1}|
S1 = pb.StemRandomParameter(plant)  
L1 = pb.LeafRandomParameter(plant)  #|\label{l2_1:p1}|

p0.name = "taproot" # |\label{l2_1:tap_start}|
p0.a = 0.2  # radius (cm) |\label{l2_1:tap_a}|
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

p1.name = "lateral" # |\label{l2_1:lat_start}|
p1.a = 0.1  # radius (cm)
p1.subType = 2  # index starts at 1
p1.lmax = 15  # maximal root length (cm) |\label{l2_1:lat_lmax}|
p1.lmaxs = 0.15  # standard deviation of the maximal root length (cm) # |\label{l2_1:lat_lmaxs}|
p1.theta = 90.0 / 180.0 * np.pi  # (rad)
p1.r = 2  # initial growth rate (cm day-1)
p1.dx = 1  # axial resolution (cm)
p1.tropismT = pb.TropismType.gravi  #
p1.tropismN = 2  # strength of tropism (1)
p1.tropismS = 0.1  # maximal bending (rad cm-1) |\label{l2_1:lat_end}|

S1.name = "stem" # |\label{l2_1:stem2_start}|
S1.subType = 1  # radius (cm)
S1.a = 0.2  # radius (cm)
S1.ldelay = -1 #delay between lateral creation and start of growth
S1.nodalGrowth = 0  # inter-lateral distance (cm)
S1.lb = 5  # basal zone (cm)
S1.la = 10  # apical zone (cm)
S1.lmax = 30  # maximal root length (cm)
S1.ln = 1.0  # inter-lateral distance (cm)
S1.r = 2  # growth rate (cm) # |\label{l2_1:stem2_end}|
S1.successorOT = [[4]]
S1.successorST = [[1]]
S1.successorP = [[1]]

L1.name = "leaf" # |\label{l2_1:leaf_start}|
L1.subType = 1  # radius (cm)
L1.ldelay = -1 #delay between lateral creation and start of growth
L1.r = 5  # growth rate (cm)# |\label{l2_1:leaf_end}|
L1.a = 0.05
L1.la, L1.lb, L1.lmax, L1.ln = 3.5, 1., 7.5, 3
L1.shapeType = 2
L1.areaMax = 10
L1.leafGeometryPhi = np.array([-90., -45, 0., 45., 90.]) / 180. * np.pi
L1.leafGeometryX = np.array([3., 2.2, 1.7, 2., 3.5])
N = 101  # N is rather high for testing
L1.createLeafRadialGeometry(L1.leafGeometryPhi, L1.leafGeometryX, N)

plant.setOrganRandomParameter(p0) # |\label{l2_1:setrp_p0}|
plant.setOrganRandomParameter(p1) # |\label{l2_1:sterp_p1}|
plant.setOrganRandomParameter(L1) # |\label{l2_1:sterp_L1}|
plant.setOrganRandomParameter(S1) # |\label{l2_1:sterp_S1}|


srp = pb.SeedRandomParameter(plant)  # with default values |\label{l2_1:seed_rp}|
srp.seedPos = pb.Vector3d(0.0, 0.0, -3.0)  # seed position (cm)
plant.setOrganRandomParameter(srp) # |\label{l2_1:srp}|

plant.initialize(False)

plant.simulate(50, False)  # |\label{l2_1:simulation_start}|
vp.plot_plant(plant, "creationTime")
plant.write("results/example2_1_parameters.vtp")# |\label{l2_1:simulation_end}|

# Some outputs....
print(srp) # |\label{l2_1:print_srp}|
print(p1) # |\label{l2_1:print_p1}|
