import sys
sys.path.append("../..")
import plantbox as pb
import rsml_reader as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt

# from maize_basal
r_tap = 2.0707344055175803  # cm/day
lmax_tap = 100.0  # cm

times = [78]  # measurement times (not in the rsml)
names = ["Maize_Kutschera.rsml"]

""" open all files, and split into measurent times"""
polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
properties = [ [ [] for t in times] for i in range(0, len(names)) ]
functions = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    print(names[i])
    j = len(times) - 1
    p, properties[i][j], functions[i][j] = rsml.read_rsml(names[i])  # read file to final time step
    print(properties[i][j].keys())
    print(functions[i][j].keys())
#     print("diameter len", len(functions[i][j]["diameter"]))
#     print("root len", len(p))
    p = [[np.array([p[i][j][k] for k in range(0, 3)]) for j in range(0, len(p[i])) ] for i in range(0, len(p)) ]
    polylines[i][j] = p
    for k in range(0, j):  # truncate the others
        polylines[i][k], properties[i][k], functions[i][k] = es.measurement_time(polylines[i][j], properties[i][j], functions[i][j], times[k])

""" add diameters as property """
for i in range(0, len(names)):
    for j in range(0, len(times)):  # truncate the others
        es.create_diameters(polylines[i][j], properties[i][j], functions[i][j])  # add diameter property

""" add order as property """
for i in range(0, len(names)):
    for j in range(0, len(times)):  # truncate the others
        es.create_order(polylines[i][j], properties[i][j])  # add order property

# polylines[i][j] contains i-th plant, j-th measurement time

""" find all base roots """
base_polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
base_properties = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    for j in range(0, len(times)):
            base_polylines[i][j], base_properties[i][j] = es.base_roots(polylines[i][j], properties[i][j])

""" find all first order laterals """
new_polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
new_properties = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    for j in range(0, len(times)):
            new_polylines[i][j], new_properties[i][j] = es.laterals1(polylines[i][j], properties[i][j])
polylines = new_polylines; properties = new_properties;
print('first order laterals done')

# print()
# for i in range(0, len(names)):
#     print("Number of nodes", len(base_polylines[i][-1][0]))  # approx 1 node per mm

# for Faba this is exactly 1 for all 5 measurements

""" plots """
# prop = properties[0][0]["length"]
# rsml.plot_rsml(polylines[0][0],prop)
# plt.show()

""" recalculate length """
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_length(polylines[i][j], properties[i][j])  # add root length for all roots
            es.create_length(base_polylines[i][j], base_properties[i][j])

""" reconstruct laterals """
for i in range(0, len(names)):
    for j in range(0, len(times)):
            # print("\nMeasurement", i, ", time", j, "number of roots", len(polylines[i][j]))
            polylines[i][j], properties[i][j] = es.reconstruct_laterals(polylines[i][j], properties[i][j], base_polylines[i][j][0], snap_radius = 0.5 * 116.93)  # pixel?
            # es.create_order(polylines[i][j], properties[i][j])  # add root order

""" measure remaining tap root parameters """
print()
la_, la_delay_, lb_, ln_, theta_, a_, a0_ = [], [], [], [], [], [], []
for i in range(0, len(names)):

    lb0, ln0, la0 = es.analyze_zones(polylines[i][0], properties[i][0])  # only for the first measurement
    la_.append(la0)
    lwa = lb0 + np.sum(ln0)  # length without apical
    la_delay_.append(times[0] - es.negexp_age(lwa, r_tap, lmax_tap))  # calculate apical zone delay

    for j in range(0, len(times)):
        a0_.append(properties[i][j]["diameter"][0])  # tap root
        a_.extend(properties[i][j]["diameter"][1:])  # laterals
        lb, ln, la = es. analyze_zones(polylines[i][j], properties[i][j])
        lb_.append(lb)
        ln_.extend(ln)
        theta = es. analyze_theta(polylines[i][j], properties[i][j])
        theta_.extend(theta)
        k = la + lb + np.sum(np.array(ln))
        k2 = es.polyline_length(0, len(polylines[i][j][0]) - 1, polylines[i][j][0])  #
        print("Measurement", i, "time", j, "lb", lb, "la", la, "ln", np.sum(np.array(ln)), "k", k, k2, "length", properties[i][j]["length"][0])

print()
print("a0 ", np.mean(a0_) / 2., "cm ", np.std(a_) / 2., "cm")  # diameter -> radius
print("a  ", np.mean(a_) / 2., "cm ", np.std(a0_) / 2., "cm")  # diameter -> radius
print("la ", np.mean(la_), "cm  ", np.std(la_), "cm")
print("lad", np.mean(la_delay_), "days", np.std(la_delay_), "days")
print("lb ", np.mean(lb_), "cm  ", np.std(lb_), "cm")
print("ln ", np.mean(ln_), "cm", np.std(ln_), "cm")
print("theta", np.mean(theta_) / np.pi * 180, "  ", np.std(theta_) / np.pi * 180)
print("Maximal number of branches", round((lmax_tap - np.mean(la_) - np.mean(lb_)) / np.mean(ln)))
print()

# """ plots histogram of ln, its zero or ln, which is rather unexpected"""
# n, bins, patches = plt.hist(ln_, 1000)
# plt.ylim([0., 100.])
# plt.xlim([0., 0.3])
# plt.show()
# print("zeros...", np.sum(np.array(ln_) == np.zeros(np.array(ln_).shape)))

# """ plots histogram of theta (not normally distributed) """
# n, bins, patches = plt.hist(theta_, 50)
# plt.show()

""" estimate lateral age """
lad0 = np.mean(la_delay_)
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_age_delay(polylines[i][j], properties[i][j], times[j], r_tap, lmax_tap, la_delay_[i])

""" root age plot """
col = ["r*", "g*", "b*", "m*", "c*"]
for j in range(0, len(times)):  # Maybe drop measurement 3 ? TODO
    age_flat, l_flat = [], []
    for i in range(0, len(names)):  # len(names)
        lengths = properties[i][j]["length"]
        ages = properties[i][j]["age"]
        for k in range(1, len(polylines[i][j])):  # skip tap root
            age_flat.append(ages[k])
            l_flat.append(lengths[k])
    plt.plot(l_flat, age_flat, col[j])
plt.ylabel("(estimatd) root age [day]")
plt.xlabel("measured root length [cm]")

r, k, f = es.fit_taproot_rk(l_flat, age_flat)
print(r, "cm/day", k, "cm; error", f)
k1 = 50.
r1, f1 = es.fit_taproot_r(l_flat, age_flat, k1)
print(r1, "cm/day", k1, "cm; error", f1)
t_ = np.linspace(0., times[-1] - lad0, 200)
y0 = es.negexp_length(t_, r, k)
plt.plot(y0, t_, "k", label = "fit r and k")
y1 = es.negexp_length(t_, r1, k1)
plt.plot(y1, t_, "k:", label = "fit r")
plt.legend()
plt.show()

print("\nfin.")

