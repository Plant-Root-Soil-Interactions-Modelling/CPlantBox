"""simple root system from scratch (without parameter files)"""

import matplotlib.pyplot as plt  # |\label{l2_1:matplotlib}|
import numpy as np  # |\label{l2_1:numpy}|

import plantbox as pb

plant = pb.Plant()
p0 = pb.RootRandomParameter(plant)  # with default values  |\label{l2_1:p0}|
p1 = pb.RootRandomParameter(plant)  # all standard deviations are 0  |\label{l2_1:p1}|

p0.name = "taproot" # |\label{l2_1:tap_start}|
p0.a = 0.2  # radius (cm) |\label{l2_1:tap_a}|
p0.subType = 1  # index starts at 1 
p0.lb = 5  # basal zone (cm)
p0.la = 10  # apical zone (cm)
p0.lmax = 30  # maximal root length (cm)
p0.ln = 1.0  # inter-lateral distance (cm)
p0.theta = 0.0  # (rad)
p0.r = 1  # initial growth rate (cm/day)
p0.dx = 10  # axial resolution (cm)
p0.successor = [[2]]  # add successors
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
p1.r = 2  # initial growth rate (cm/day)
p1.dx = 1  # axial resolution (cm)
p1.tropismT = pb.TropismType.gravi  # exo
p1.tropismN = 2  # strength of tropism (1)
p1.tropismS = 0.1  # maximal bending (rad/cm) |\label{l2_1:lat_end}|

plant.setOrganRandomParameter(p0) # |\label{l2_1:setrp_p0}|
plant.setOrganRandomParameter(p1) # |\label{l2_1:sterp_p1}|

srp = pb.SeedRandomParameter(plant)  # with default values |\label{l2_1:seed_rp}|
srp.seedPos = pb.Vector3d(0.0, 0.0, -3.0)  # seed position (cm)
srp.maxB = 0  # number of basal roots (neglecting basal roots and shoot borne) (1)
srp.firstB = 10.0  # first emergence of a basal root (day)
srp.delayB = 3.0  # delay between the emergence of basal roots (day)
plant.setOrganRandomParameter(srp) # |\label{l2_1:srp}|

plant.initialize()

fig, axes = plt.subplots(1, 3, figsize=(15, 7))
simtimes = [0, 30, 60, 125]  # the last lateral will emerge at last time step |\label{l2_1:simtime}|
for i in range(0, 3): # |\label{l2_1:simloop_start}|
    plant.simulate(np.diff(simtimes)[i])  # (day)
    a = axes[i]
    a.set_xlim([-15, 15.0])  # x axis limit |\label{l2_1:x_lim}|
    a.set_ylim([-35.0, 0.0]) # y axis limit |\label{l2_1:y_lim}|
    a.set_title(f"after {plant.getSimTime()} days") # plot title |\label{l2_1:title}|
    roots = plant.getPolylines() # |\label{l2_1:polylines}|
    for root in roots: # |\label{l2_1:coord_start}|
        for j, n in enumerate(root[:-1]): 
            n2 = root[j + 1]
            a.plot([n.x, n2.x], [n.z, n2.z], "g") # |\label{l2_1:simloop_end}|

fig.tight_layout()
plt.show()

plant.write("/results/topics_parameters2.vtp")

# Some outputs....
print(" length", plant.getParameter("length")) # get length of roots |\label{l2_1:root_len}|
print("    age", plant.getParameter("age"))  # get age of roots |\label{l2_1:root_age}|
print("subType", plant.getParameter("subType")) # get subTypes of roots |\label{l2_1:root_subType}|
print("     la", plant.getParameter("la")) # get apical length of roots |\label{l2_1:root_la}|
print("la_mean", plant.getParameter("la_mean")) # get apical length mean of roots |\label{l2_1:root_la_mean}|
print(" radius", plant.getParameter("radius")) # get radius of roots |\label{l2_1:root_a}|
