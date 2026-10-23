"""Define the growth delays of organs"""

import numpy as np

import plantbox as pb
import plantbox.visualisation.vtk_animate as va
import plantbox.visualisation.vtk_plot as vp


def show_me(plant, simTime=10):
    """makes an animation and a plot visualising plant growth delays"""
    dt = 0.1
    N = int(simTime / dt)
    min_ = np.array([-20, -20, -50])
    max_ = np.array([20, 20, 30.0])

    plant.initialize(False)

    anim = va.AnimateRoots(plant)
    anim.min = min_
    anim.max = max_
    anim.res = [1, 1, 1]
    anim.plant = True
    anim.param_name = "creationTime"
    anim.start()

    for i in range(0, N):
        plant.simulate(dt, False)
        anim.update()
    vp.plot_plant(plant, "creationTime")  # "creationTime")


path = pb.data_path() + "/structural/plant/"

# ---------------------------------------------------------------------------
# Part A - Mixed distance-based delays for roots, time-based for shoot (self-defined fixed delay)
# ---------------------------------------------------------------------------
plant = pb.MappedPlant(2)  # |\label{l2_1_3:plantStart}|
plant.readParameters(path + "example2_1_3.xml")
srp = plant.getOrganRandomParameter(pb.seed)[0]

srp.delayDefinitionRoot = pb.dd_distance  # roots: laterals wait for the apical zone to emerge  |\label{l2_1_3:ddRoot}|
srp.delayDefinitionShoot = pb.dd_time_self  # shoot: time-based value for all laterals |\label{l2_1_3:ddShoot}|

# Leaf subType 1
lrp = plant.getOrganRandomParameter(pb.leaf)[1]
lrp.ldelay = 3  # |\label{l2_1_3:leaf_ldelayA}|
lrp.ldelays = 0.0  # |\label{l2_1_3:leaf_ldelaysA}|

show_me(plant, simTime=15)

tap_root = plant.getOrgans(pb.root)[0]  # |\label{l2_1_3:analyticald}|
meanLn = tap_root.getParameter("lnMean")  # mean inter-lateral distance
effectiveLa = max(tap_root.getParameter("la") - meanLn / 2, 0.0)  # effective apical distance, observed apical distance is in (la-ln/2, la+ln/2)
ageLN = tap_root.calcAge(tap_root.getLength(True))  # theoretical age of root when lateral node is created
ageLG = tap_root.calcAge(tap_root.getLength(True) + effectiveLa)  # age of the root, when the lateral starts growing (i.e when the apical zone is developed)
growthDelay = ageLG - ageLN  # time the lateral has to wait #|\label{l2_1_3:analyticaldend}|

print("\n=== Part A: Distance- (for roots) and time-  (for shoot) based delay ===")
print(f"  Root delay : {growthDelay:.2f} days  (or {effectiveLa:.2f} cm)")
print(f"  Leaf delay : {lrp.ldelay:.2f} days  (std {lrp.ldelays:.2f})")


# # ---------------------------------------------------------------------------
# # Part B Time-based for all (self-defined fixed delay)
# # ---------------------------------------------------------------------------
# plant = pb.MappedPlant(2)
# plant.readParameters(path + "example2_1_3.xml")

# srp = plant.getOrganRandomParameter(pb.seed)[0]
# srp.delayDefinitionRoot = pb.dd_time_self  # roots: Apcial zone decides delay      |\label{l2_1_3:ddRootB}|
# srp.delayDefinitionShoot = pb.dd_time_self  # shoot: lateral decides its own delay       |\label{l2_1_3:ddShootB}|

# # Root subType 2 (laterals)
# rrp = plant.getOrganRandomParameter(pb.root)[2]
# rrp.ldelay = 2.5  # mean lateral delay [days]   |\label{l2_1_3:root_ldelayB}|
# rrp.ldelays = 0.0  # std dev of delay  [days]    |\label{l2_1_3:root_ldelaysB}|

# # Leaf subType 1
# lrp = plant.getOrganRandomParameter(pb.leaf)[1]
# lrp.ldelay = 1  # |\label{l2_1_3:leaf_ldelayB}|
# lrp.ldelays = 0.0  # |\label{l2_1_3:leaf_ldelaysB}|

# show_me(plant, simTime=15)

# print("\n=== Part B: Time-based delay ===")
# print(f"  Root delay : {rrp.ldelay:.2f} days (std {rrp.ldelays:.2f})")
# print(f"  Leaf delay : {lrp.ldelay:.2f} days  (std {lrp.ldelays:.2f})")


# # ---------------------------------------------------------------------------
# # Part C - Mixed distance-based delays for roots, time-based for shoots
# # ---------------------------------------------------------------------------
# plant = pb.MappedPlant(2)
# plant.readParameters(path + "example2_1_3.xml")

# srp = plant.getOrganRandomParameter(pb.seed)[0]
# srp.delayDefinitionRoot = pb.dd_distance  # roots: use distance-based delay       |\label{l2_1_3:ddRootC}|
# srp.delayDefinitionShoot = pb.dd_time_self  # shoot: lateral decides its own delay       |\label{l2_1_3:ddShootC}|

# # Leaf subType 1: shoot laterals wait 5 days.
# lrp = plant.getOrganRandomParameter(pb.leaf)[1]
# lrp.ldelay = 2.5  # |\label{l2_1_3:leaf_ldelayC}|
# lrp.ldelays = 0.0  # |\label{l2_1_3:leaf_ldelaysC}|

# # Stem subType 1 (main stem)
# strp = plant.getOrganRandomParameter(pb.stem)[1]
# strp.delayNGStart = 1  # |\label{l2_1_3:delayNGStart}|
# strp.delayNGEnd = 10  # |\label{l2_1_3:delayNGEnd}|

# show_me(plant)

# tap_root = plant.getOrgans(pb.root)[0]
# meanLn = tap_root.getParameter("lnMean")  # mean inter-lateral distance
# effectiveLa = max(tap_root.getParameter("la") - meanLn / 2, 0.0)  # effective apical distance, observed apical distance is in (la-ln/2, la+ln/2)
# ageLN = tap_root.calcAge(tap_root.getLength(True))  # theoretical age of root when lateral node is created
# ageLG = tap_root.calcAge(tap_root.getLength(True) + effectiveLa)  # age of the root, when the lateral starts growing (i.e when the apical zone is developed)
# growthDelay = ageLG - ageLN  # time the lateral has to wait

# print("\n=== Part C: Mixed distance (roots) + self-defined (shoot) ===")
# print(f"  Root delay : {growthDelay:.2f} days  (effectiveLa = {effectiveLa:.2f} cm)")
# print(f"  Leaf delay : {lrp.ldelay} days  (std {lrp.ldelays})")
