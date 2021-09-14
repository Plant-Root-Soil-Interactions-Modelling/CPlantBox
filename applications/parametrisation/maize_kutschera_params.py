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
    p0.a, p0.a_s = p[3, 0] / 116.93, p[3, 1] / 116.93 / 10.  # [pixel] -> [cm] radius, reduce std
    p0.theta, p0.thetas = p[4, 0], p[4, 1]  # [rad]
    p0.subType = order + 1  # [-] index starts at 1
    if successors:
        p0.successor = [order + 2]  # add successors
        p0.successorP = [1]  # probability that successor emerges


print("reading")
lmax = 100.  # [cm] fixed
time = 75  # measurement times (not in the rsml)
name = "RSML/Maize_Kutschera.rsml"
roots = es.parse_rsml(name, time)

#
# finding parameters la, lb, ln, a, theta: see maize_kutschera.py
#
print("finding parameters la, lb, ln, a, theta")
max_order = np.max(np.array([r.order() for r in roots.values()], dtype = np.int64))
print("Number of root orders: ", max_order + 1)

roots_i = []  # sort roots per oreder
for i in range(0, max_order + 1):
    roots_i.append(es.get_order(i, roots))

params = []
for i in range(0, 2):
    params.append(es.get_params(roots_i[i], time))

roots_ = []
for j in range(2, max_order + 1):  # all roots > order 1
    roots_.extend(roots_i[j])
params.append(es.get_params(roots_, time))

#
# see maize_kutschera2.py
#
print()
print("\nfinding production rate and growth rate of order 0")
rate, r0, r0s = es.estimate_set_order0_rate(roots, lmax, time)
print("rate", rate, "r0", r0, "see", r0s)

print("\nfinding growth rate and maximal length of order 1")
order1 = es.get_order(1, roots)
lengths1 = np.array([r.length() for r in order1])  # for last measurement
ages1 = np.array([r.ages[time] for r in order1])

# res, f = es.estimate_rk(lengths1, ages1)
# r1, lmax1 = res.x[0], res.x[1]
# r1s, lmaxs1 = es.get_see_rk(lengths1, ages1, r1, lmax1)

lmax1 = 5
res, f = es.estimate_r(lengths1, ages1, lmax1)
r1 = res.x[0]
r1s, lmaxs1 = 0.1 * r1, 0.1 * lmax1  # es.get_see_rk(lengths1, ages1, r1, lmax1)

print("r1", r1, "see", r1s, "lmax1", lmax1, "see", lmaxs1)

#
# build xml
#
print("\nwriting kutschera.xml\n")
rs = pb.RootSystem()
p0 = pb.RootRandomParameter(rs)  # with default values,
p1 = pb.RootRandomParameter(rs)  # all standard deviations are 0
p2 = pb.RootRandomParameter(rs)  # all standard deviations are 0
srp = pb.SeedRandomParameter(rs)  # with default values

srp.name = "maize kutschera"
srp.seedPos = pb.Vector3d(0., 0., -3.)  # [cm] seed position
srp.maxB = 100  # [-] number of basal roots (neglecting basal roots and shoot borne)
srp.firstB = rate  # [day] first emergence of a basal root
srp.delayB = rate  # [day] delay between the emergence of basal roots

p0.name = "base roots"
insert_params(p0, params[0], 0)  # inserts la, lb, ln, a, theta
p0.r = r0  # [cm/day] initial growth rate
p0.rs = r0s
p0.lmax = lmax  # [cm] maximal root length, number of lateral branching nodes = round((lmax-lb-la)/ln) + 1
# not based on data
p0.dx = 0.5  # [cm] axial resolution
p0.tropismT = pb.TropismType.gravi  #
p0.tropismN = 1  # [-] strength of tropism
p0.tropismS = 0.2  # [rad/cm] maximal bending

p1.name = "first order laterals"
insert_params(p1, params[1], 1)  # inserts la, lb, ln, a, theta
p1.r = r1  # [cm/day] initial growth rate
p1.rs = r1s  # standard deviation of initial growth rate
p1.lmax = lmax1  # # [cm] apical zone
p1.lmaxs = lmaxs1  # # standard deviation of apical zone
# not based on data
p1.dx = 0.25  # [cm] axial resolution
p1.tropismT = pb.TropismType.gravi  # exo
p1.tropismN = 0.5  # [-] strength of tropism
p1.tropismS = 0.3  # [rad/cm] maximal bending

p2.name = "higher order laterals"
insert_params(p2, params[2], 2, False)  # inserts la, lb, ln, a, theta
p2.r = r1  # [cm/day]  initial growth rate
p2.rs = r1s  # standard deviation of initial growth rate
p2.lmax = lmax1  # # [cm] apical zone
p2.lmaxs = lmaxs1  # # standard deviation of apical zone
# not based on data
p2.dx = 0.25  # [cm] axial resolution
p2.tropismT = pb.TropismType.gravi  # exo
p2.tropismN = 0.5  # [-] strength of tropism
p2.tropismS = 0.3  # [rad/cm] maximal bending

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)
rs.setOrganRandomParameter(p2)
rs.setRootSystemParameter(srp)

print(p0)
print()
print(p1)
print()
print(p2)
print()
print(srp)
print()

rs.writeParameters("kutschera.xml")

