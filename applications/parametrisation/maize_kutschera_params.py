"""everything from scratch (without parameter files)"""
import sys; sys.path.append("../..")

import matplotlib.pyplot as plt
import numpy as np

import rsml_reader as rsml
import estimate_root_params as es
import plantbox as pb


def insert_params(p0, p :np.array, order : int, successors :bool = True):
    """ inserts la, lb, ln, a, theta obtained by estimate_root_params into the RootRandomParameter object """
    p0.la, p0.las = p[0, 0], p[0, 1]  # [cm] apical zone
    p0.lb, p0.lbs = p[1, 0], p[1, 1]  # [cm] basal zone
    p0.ln, p0.lns = p[2, 0], p[2, 1]  # [cm] inter-lateral distance (16 branching nodes)
    p0.a, p0.a_s = p[3, 0], p[3, 1]  # [cm] radius
    p0.theta, p0.thetas = p[4, 0], p[4, 1]  # [rad]
    p0.subType = order + 1  # [-] index starts at 1
    if successors:
        p0.successor = [order + 1]  # add successors
        p0.successorP = [1]  # probability that successor emerges


print("reading")
time = [75]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

#
# see maize_kutschera.py
#
print("finding parameters la, lb, ln, a, theta")
max_order = np.max(np.array([r.order() for r in roots.values()], dtype = np.int64))
print("Number of root orders: ", max_order + 1)

roots_i = []
for i in range(0, max_order + 1):
    roots_i.append(es.get_order(i, roots))

params = []
for i in range(0, 3):
    params.append(es.get_params(roots_i[i], time[-1]))

params2 = []
for i in range(1, min(max_order + 1, 4)):
    roots_ = []
    for j in range(i, max_order + 1):
        roots_.extend(roots_i[j])
    params2.append(es.get_params(roots_, time[-1]))

#
# see maize_kutschera.py
#
print("finding production rate, growth rate, and maximal root length")
rate = 1

#
# build xml
#
rs = pb.RootSystem()
p0 = pb.RootRandomParameter(rs)  # with default values,
p1 = pb.RootRandomParameter(rs)  # all standard deviations are 0
p2 = pb.RootRandomParameter(rs)  # all standard deviations are 0
srp = pb.SeedRandomParameter(rs)  # with default values

srp.seedPos = pb.Vector3d(0., 0., -3.)  # [cm] seed position
srp.maxB = 100  # [-] number of basal roots (neglecting basal roots and shoot borne)
srp.firstB = rate  # [day] first emergence of a basal root
srp.delayB = rate  # [day] delay between the emergence of basal roots

p0.name = "base roots"
insert_params(p0, params[0], 0)
p0.r = 1  # [cm/day] initial growth rate
p0.lmax = 30  # [cm] maximal root length, number of lateral branching nodes = round((lmax-lb-la)/ln) + 1
p0.dx = 10  # [cm] axial resolution
p0.tropismT = pb.TropismType.gravi  #
p0.tropismN = 1.8  # [-] strength of tropism
p0.tropismS = 0.2  # [rad/cm] maximal bending

p1.name = "first order laterals"
insert_params(p1, params[1], 1)
p1.r = 2  # initial growth rate
p1.lmax = 15  # # [cm] apical zone
p1.lmaxs = 0.15  # [cm] standard deviation of the apical zone
p1.dx = 1  # [cm] axial resolution
p1.tropismT = pb.TropismType.gravi  # exo
p1.tropismN = 2  # [-] strength of tropism
p1.tropismS = 0.1  # [rad/cm] maximal bending

p2.name = "higher order laterals"
insert_params(p2, params2[1], 2, False)
p2.r = 2  # initial growth rate
p2.lmax = 15  # # [cm] apical zone
p2.lmaxs = 0.15  # [cm] standard deviation of the apical zone
p2.dx = 1  # [cm] axial resolution
p2.tropismT = pb.TropismType.gravi  # exo
p2.tropismN = 2  # [-] strength of tropism
p2.tropismS = 0.1  # [rad/cm] maximal bending

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)
rs.setOrganRandomParameter(p2)
rs.setRootSystemParameter(srp)
rs.writeParameters("kutschera.xml")
