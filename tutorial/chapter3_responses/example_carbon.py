"""
revised 9.2024
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import numpy as np


def check_tfs(plant):
    """ checks if plant tropisms have a parent plant and plots parent plant id """
    for p in plant.getOrganRandomParameter(pb.root):
        if p.f_tf.isExpired():
            raise Exception("root with subType " + str(p.subType) + " has tropism without parent plant")
        print("root with subType {:g}, tropism {:s} is owned by plant {:g} of {:g}".format(p.subType, str(type(p.f_tf)), p.f_tf.getPlant().plantId, plant.plantId))

    for p in plant.getOrganRandomParameter(pb.stem):
        if p.f_tf.isExpired():
            raise Exception("stem with subType " + str(p.subType) + " has tropism without parent plant")
        print("stem with subType {:g}, tropism {:s} is owned by plant {:g} of {:g}".format(p.subType, str(type(p.f_tf)), p.f_tf.getPlant().plantId, plant.plantId))

    for p in plant.getOrganRandomParameter(pb.leaf):
        if p.f_tf.isExpired():
            raise Exception("leaf with subType " + str(p.subType) + "has tropism without parent plant")
        print("leaf with subType {:g}, tropism {:s} is owned by plant {:g} of {:g}".format(p.subType, str(type(p.f_tf)), p.f_tf.getPlant().plantId, plant.plantId))


def elongate(rs, dt, inc, se):

    accuracy = 0.001  # cm
    maxiter = 20

    ol = rs.getSummed("length")
    i = 0

    rs_ = rs.copy()

    se.setScale(1.)
    rs_.simulate(dt, True)
    inc_ = rs_.getSummed("length") - ol

    print("expected increase is ", inc_, " maximum is ", inc)
    if inc_ > inc and abs(inc_ - inc) > accuracy:  # check if we have to perform a binary search

        sl = 0.  # left
        sr = 1.  # right

        while abs(inc_ - inc) > accuracy and i < maxiter:  # binary search

            print("\nIteration", i)
            m = (sl + sr) / 2.  # mid
            rs_ = rs.copy()
            # print("copy ")
            # check_tfs(rs_)
            # print("original")
            # check_tfs(rs)
            se.setScale(m)
            rs_.simulate(dt, True)
            inc_ = rs_.getSummed("length") - ol
            print(f"{i}\tsl, mid, sr ", sl, m, sr, "inc", inc_, "err", abs(inc_ - inc))

            if inc_ > inc:  # concatenate
                sr = m
            else:
                sl = m

            i += 1

    rs.simulate(dt)


# Parameter
simtime = 50  # days
dt = 1  # day
N = round(simtime / dt)  # steps
maxinc = 20;  # 20  # maximal length increment (cm/day), TODO base this value on some fancy model
maxvol = 0.1;  # cm3

# Initialize root system
rs = pb.Plant()
name = "../../modelparameter/structural/rootsystem/Zea_mays_4_Leitner_2014"
rs.readParameters(name + ".xml")

# Set up depth dependent elongation scaling function
scale_elongation = pb.EquidistantGrid1D(0, -50, 100)  # for root elongation from 0 cm to -50 cm, 100 nodes = 99 layers
soil_strength = np.ones((99,)) * 0.5  # some data for the 99 layers, e.g. np.linspace(0.1, 1., 99)
scales = np.exp(-0.4 * soil_strength)  # scales from some equation (scale = function(soil_strength) ), where scale in (0,1)
scale_elongation.data = scales  # set proportionality factors

# Proportionally scale this function
se = pb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)

# Manually set scaling function
for p in rs.getOrganRandomParameter(pb.OrganTypes.root):  # rs.getOrganRandomParameter(pb.OrganTypes.root)
    p.f_se = se

rs.initialize()
# check_tfs(rs)

ol = 0

print("\n\n*****************")

# Simulation loop
for i in range(0, N):

    print("\nSimulation step", i)

    # if maxinc is dynamic: set maxinc (cm/day) according to some model

    # if soil_strength is dynamic: update soil_strength according to some model (update like in L58-L60)

    # elongate(rs, dt, maxinc, se)  #  for debugging
    # check_tfs(rs)

    rs.simulate(dt, maxinc, se, False)
    check_tfs(rs)

    # rs.simulateLimited(dt, maxvol, "volume", [0., 0., 0.1, 0.1, 0.1], se, True)  # { ot_organ = 0, ot_seed = 1, ot_root = 2, ot_stem = 3, ot_leaf = 4 };
    # "length", "lenghtTh"

    l = rs.getSummed("length")
    inc = l - ol
    ol = l
    print("elongated for", inc, " cm")

print("fin.")
