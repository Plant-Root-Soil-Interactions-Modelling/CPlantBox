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

srp.delayDefinition      = 0  # laterals wait for the apical zone to pass  |\label{l2_1_3:ddRoot}|
srp.delayDefinitionShoot = 0  # same rule for shoot laterals                |\label{l2_1_3:ddShoot}|

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

print("\n=== Part A: Distance-based delay ===")
print(f"  Root delay : {growthDelay:.2f} days  (or {effectiveLa:.2f} cm)")

# ---------------------------------------------------------------------------
# Part B - Time-based delay (defined for the laterals)
# ---------------------------------------------------------------------------


plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_1_2.xml")

srp = plant.getOrganRandomParameter(pb.seed)[0]
srp.delayDefinition      = 1  # parent decides delay for root laterals    |\label{l2_1_3:ddRootB}|
srp.delayDefinitionShoot = 1  # parent decides delay for shoot laterals   |\label{l2_1_3:ddShootB}|


# Set ldelay on the parent organ random parameters.
# Root subType 1 (taproot): laterals wait 3 days on average.
rrp = plant.getOrganRandomParameter(pb.root)[1]
rrp.ldelay  = 2.5  # mean lateral delay [days]   |\label{l2_1_3:root_ldelayB}|
rrp.ldelays = 0.   # std dev of delay  [days]    |\label{l2_1_3:root_ldelaysB}|


# Stem subType 1 (main stem): shoot laterals wait 5 days.
srp = plant.getOrganRandomParameter(pb.stem)[1]
srp.ldelay  = 2.5   # |\label{l2_1_3:stem_ldelayB}|
srp.ldelays = 0.    # |\label{l2_1_3:stem_ldelaysB}|

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
vp.plot_plant(plant, "creationTime") 

print("\n=== Part B: Time-based delay (defined for the laterals) ===")
print(f"  Root    delay : {rrp.ldelay} days (std {rrp.ldelays})")
print(f"  Leaf    delay : {srp.ldelay} days (std {srp.ldelays})")

# ---------------------------------------------------------------------------
# Part C - Mixed distance-based delays for roots, time-based for shoot (self-defined fixed delay)
# ---------------------------------------------------------------------------


plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_x.xml")

srp = plant.getOrganRandomParameter(pb.seed)[0]
srp.delayDefinition      = 0  # roots: use distance-based delay       |\label{l2_1_3:ddRootC}|
srp.delayDefinitionShoot = 2  # shoot: lateral decides its own delay       |\label{l2_1_3:ddShootC}|



# Leaf subType 1: shoot laterals wait 5 days.
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


tap_root = plant.getOrgans(pb.root)[0]
meanLn = tap_root.getParameter("lnMean") # mean inter-lateral distance
effectiveLa = max(tap_root.getParameter("la") - meanLn / 2, 0.) # effective apical distance, observed apical distance is in (la-ln/2, la+ln/2)
ageLN = tap_root.calcAge(tap_root.getLength(True)) # theoretical age of root when lateral node is created
ageLG = tap_root.calcAge(tap_root.getLength(True) + effectiveLa) # age of the root, when the lateral starts growing (i.e when the apical zone is developed)
growthDelay = ageLG - ageLN # time the lateral has to wait

print("\n=== Part C: Mixed – distance (roots) + self-defined (shoot) ===")
print(f"  Root delay : {growthDelay:.2f} days  (effectiveLa = {effectiveLa:.2f} cm)")
print(f"  Leaf delay : {lrp.ldelay} days  (std {lrp.ldelays})")


