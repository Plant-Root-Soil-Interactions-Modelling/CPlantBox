""" determines the (more easy) parameters la, lb, ln, a, theta by order and prints mean and std"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import rsml.rsml_reader as rsml
import visualisation.vtk_plot as vp
import estimate_root_params as es
import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt


def insert_params(p0, p:np.array, order: int, successors:bool = True):
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


time = [1, 2, 3, 4, 5, 6, 8, 9]  # measurement times (not in the rsml)
name = ["RSML/m1/dicot/lupin/lupin_d{:g}.rsml".format(a) for a in time]

roots = es.parse_rsmls(name, np.array(time) + 2)
print("\nNumber of roots", len(roots.values()))

#
# finding parameters la, lb, ln, a, theta: see maize_kutschera.py
#
orders = np.array([r.order() for r in roots.values()], dtype = np.int64)
max_order = np.max(orders)
print("Number of root orders: ", max_order + 1, "\n")

roots_i = []
for i in range(0, max_order + 1):
    roots_i.append(es.get_order(i, roots))

params = []
for i in range(0, max_order + 1):
    print("order {:g}:".format(i), len(roots_i[i]), "roots")
    params.append(es.get_params(roots_i[i], time[-1]))
    print(params[i])
    print()

#
# find r and k for root order 0 and order 1
#
order0 = es.get_order(0, roots)
print(len(order0), "base_roots")
for r in order0:
    r.set_emergence_time(0.)
order0_lengths = np.array([r.length() for r in order0])
order0_ages = np.array([list(r.ages.values())[0] for r in order0])
res, f0 = es.estimate_rk(order0_lengths, order0_ages)
r0, lmax0 = res.x[0], res.x[1]
print("r0", r0, "lmax0", lmax0)
r0s, lmax0s = es.get_see_rk(order0_lengths, order0_ages, r0, lmax0)
print("r0s", r0s, "lmaxs0", lmax0s)

# set order 0 growth rate and maximal length
for r in order0:
    roots[r.id].calc_growth_rate(r0, lmax0)

# t_ = np.linspace(0, time[-1], 200)
# y1 = es.Root.negexp_length(t_, r0, lmax0)
# plt.plot(t_, y1)
# plt.scatter(order0_ages, order0_lengths)
# plt.show()

order1 = es.get_order(1, roots)
print(len(order1), "1st order roots")
order1_lengths = np.array([r.length() for r in order1])
order1_ages = np.array([list(r.ages.values())[0] for r in order1])
lmax1 = 15  # cm (TODO literature value)?
res, f1 = es.estimate_r(order1_lengths, order1_ages, lmax1)
r1 = res.x[0]
print("r1", r1, "lmax", lmax1)
r1s, lmaxs1 = es.get_see_rk(order1_lengths, order1_ages, r1, lmax1)
print("r1s", r1s, "lmaxs1", lmaxs1)

t_ = np.linspace(0, time[-1], 200)
y1 = es.Root.negexp_length(t_, r1, lmax1)
plt.plot(t_, y1)
plt.scatter(order1_ages, order1_lengths)
plt.xlabel("age (days)"); plt.ylabel("length (cm)")
plt.show()

#
# build xml
#
print("\nwriting m1_lupine.xml\n")
rs = pb.RootSystem()
p0 = pb.RootRandomParameter(rs)  # with default values,
p1 = pb.RootRandomParameter(rs)  # all standard deviations are 0
srp = pb.SeedRandomParameter(rs)  # with default values

srp.name = "m1 lupine"

p0.name = "base roots"
insert_params(p0, params[0], 0)  # inserts la, lb, ln, a, theta
p0.r = r0  # [cm/day] initial growth rate
p0.rs = r0s
p0.lmax = lmax0  # [cm] maximal root length, number of lateral branching nodes = round((lmax-lb-la)/ln) + 1
p0.lmaxs = lmax0s
# not based on data
p0.dx = 0.5  # [cm] axial resolution
p0.tropismT = pb.TropismType.gravi  #
p0.tropismN = 2  # [-] strength of tropism
p0.tropismS = 0.2  # [rad/cm] maximal bending

p1.name = "first order laterals"
insert_params(p1, params[1], 1, False)  # inserts la, lb, ln, a, theta
p1.r = r1  # [cm/day] initial growth rate
p1.rs = 0.1 * r1  # r1s  # standard deviation of initial growth rate
p1.lmax = lmax1  # # [cm] apical zone
p1.lmaxs = 0  # lmaxs1  # # standard deviation of apical zone
# not based on data
p1.dx = 0.25  # [cm] axial resolution
p1.tropismT = pb.TropismType.gravi  # exo
p1.tropismN = 0.5  # [-] strength of tropism
p1.tropismS = 0.3  # [rad/cm] maximal bending

rs.setOrganRandomParameter(p0)
rs.setOrganRandomParameter(p1)
rs.setRootSystemParameter(srp)

print(p0)
print()
print(p1)
print()
print(srp)
print()

rs.writeParameters("m1_lupine.xml")

print("done")

