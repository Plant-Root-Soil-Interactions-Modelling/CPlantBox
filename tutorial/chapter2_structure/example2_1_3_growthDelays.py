"""Define the growth delays of organs"""
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
import plantbox.visualisation.vtk_animate as va
import numpy as np

path = "../../modelparameter/structural/plant/"
simTime = 10
dt = 0.1
N = int(simTime/dt)
min_ = np.array([-20, -20, -50])
max_ = np.array([20, 20, 30.])

# ---------------------------------------------------------------------------
# Part A - Distance-based delay
# ---------------------------------------------------------------------------
plant = pb.MappedPlant(2)  #|\label{l2_1_3:plantStart}|
plant.readParameters(path + "example2_1_3.xml")
srp = plant.getOrganRandomParameter(pb.seed)[0]

srp.delayDefinitionRoot  = 0  # roots: laterals wait for the apical zone to pass  |\label{l2_1_3:ddRoot}|
srp.delayDefinitionShoot = 2  # shoot: time-based delay     |\label{l2_1_3:ddShoot}|

# Leaf subType 1
lrp = plant.getOrganRandomParameter(pb.leaf)[1]
lrp.ldelay  = 2.5  # |\label{l2_1_3:leaf_ldelayC}|
lrp.ldelays = 0.   # |\label{l2_1_3:leaf_ldelaysC}|

plant.initialize(False)  #|\label{l2_1_3:init}|

anim = va.AnimateRoots(plant) #|\label{l2_1_3:growth}|
anim.min = min_
anim.max = max_
anim.res = [1, 1, 1]
anim.plant = True
anim.start()

for i in range(0, N):
    plant.simulate(dt, False)
    anim.update()
vp.plot_plant(plant, "creationTime")  #|\label{l2_1_3:growthEnd}|

tap_root = plant.getOrgans(pb.root)[0] #|\label{l2_1_3:analyticald}|
meanLn = tap_root.getParameter("lnMean") # mean inter-lateral distance
effectiveLa = max(tap_root.getParameter("la") - meanLn / 2, 0.) # effective apical distance, observed apical distance is in (la-ln/2, la+ln/2)
ageLN = tap_root.calcAge(tap_root.getLength(True)) # theoretical age of root when lateral node is created
ageLG = tap_root.calcAge(tap_root.getLength(True) + effectiveLa) # age of the root, when the lateral starts growing (i.e when the apical zone is developed)
growthDelay = ageLG - ageLN # time the lateral has to wait #|\label{l2_1_3:analyticaldend}|

print("\n=== Part A: Distance- (for roots) and time-  (for shoot) based delay ===")
print(f"  Root delay : {growthDelay:.2f} days  (or {effectiveLa:.2f} cm)")
print(f"  Leaf delay : {lrp.ldelay:.2f} days  (std {lrp.ldelays:.2f})")


# ---------------------------------------------------------------------------
# Part B - Mixed distance-based delays for roots, time-based for shoot (self-defined fixed delay)
# ---------------------------------------------------------------------------
plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_1_3.xml")

srp = plant.getOrganRandomParameter(pb.seed)[0]
srp.delayDefinitionRoot  = 2  # roots: time-based delay      |\label{l2_1_3:ddRootC}|
srp.delayDefinitionShoot = 2  # shoot: lateral decides its own delay       |\label{l2_1_3:ddShootC}|

# Root subType 2 (laterals)
rrp = plant.getOrganRandomParameter(pb.root)[2]
rrp.ldelay  = 2.5  # mean lateral delay [days]   |\label{l2_1_3:root_ldelayB}|
rrp.ldelays = 0.   # std dev of delay  [days]    |\label{l2_1_3:root_ldelaysB}|

# Stem subType 1 (main stem)
srp = plant.getOrganRandomParameter(pb.stem)[1]
srp.delayNGStart = 3  # |\label{l2_1_3:delayNGStart}|
srp.delayNGEnd   = 7  # |\label{l2_1_3:delayNGEnd}|

# Leaf subType 1
lrp = plant.getOrganRandomParameter(pb.leaf)[1]
lrp.ldelay  = 2.5  # |\label{l2_1_3:leaf_ldelayC}|
lrp.ldelays = 0.   # |\label{l2_1_3:leaf_ldelaysC}|

plant.initialize(False)

anim = va.AnimateRoots(plant)
anim.min = min_
anim.max = max_
anim.res = [1, 1, 1]
anim.plant = True
anim.start()

for i in range(0, N):
    plant.simulate(dt, False)
    anim.update()
vp.plot_plant(plant, "creationTime") #"creationTime")


print("\n=== Part B: Time-based delay ===")
print(f"  Root delay : {rrp.ldelay:.2f} days (std {rrp.ldelays:.2f})")
print(f"  Leaf delay : {lrp.ldelay:.2f} days  (std {lrp.ldelays:.2f})")


