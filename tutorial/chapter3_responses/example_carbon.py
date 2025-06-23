"""
revised 9.2024
"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

import numpy as np

# Parameter
simtime = 50  # days
dt = 1
N = round(simtime / dt)  # steps
maxinc = 20;  # maximal length increment (cm/day), TODO base this value on some fancy model

# Initialize root system
rs = pb.RootSystem()
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
for p in rs.getRootRandomParameter():  # rs.getOrganRandomParameter(pb.OrganTypes.root)
    p.f_se = se

rs.initialize()

ol = 0

# Simulation loop
for i in range(0, N):

    print("\nSimulation step", i)

    # if maxinc is dynamic: set maxinc (cm/day) according to some model

    # if soil_strength is dynamic: update soil_strength according to some model (update like in L58-L60)

    rs.simulate(dt, maxinc, se, True)  # True = disable debug messages, False = enable debug messages

    l = rs.getSummed("length")
    inc = l - ol
    ol = l
    print("elongated for", inc, " cm")

rs.write("results/example_carbon.vtp")

# # DEPRICATED draft  was added to the c++ code
# accuracy = 0.1  # cm
# maxiter = 10
#
#
# #
# def elongate(rs, inc, dt, se):
#
#     ol = rs.getSummed("length")
#     i = 0
#
#     rs_ = rb.RootSystem(rs)  # copy ################# SHALLOW?!
#     se.setScale(1.)
#     rs_.simulate(dt, True)
#     inc_ = rs_.getSummed("length") - ol
#
#     if inc_ > inc and abs(inc_ - inc) > accuracy:  # check if we have to perform a binary search
#
#         sl = 0.  # left
#         sr = 1.  # right
#
#         while abs(inc_ - inc) > accuracy and i < maxiter:  # binary search
#
#             m = (sl + sr) / 2.  # mid
#             rs_ = rb.RootSystem(rs)    # copy ################# SHALLOW?!
#             se.setScale(m)
#             rs_.simulate(dt, True)
#             inc_ = rs_.getSummed("length") - ol
#             print("\tsl, mid, sr ", sl, m, sr, inc_)
#
#             if inc_ > inc:  # concatenate
#                 sr = m
#             else:
#                 sl = m
#
#             i += 1
#
#         return rs_
#
#     else:
#         return rs_

