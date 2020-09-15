import sys;  sys.path.append("../..")
""" determines growth rate and production rate and ages of the basal roots """

import rsml_reader as rsml
import estimate_root_params as es
import numpy as np
import matplotlib.pyplot as plt

time = [75]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

# basals = es.get_order(0, roots)
# basal_ids = np.array([r.id for r in basals], dtype = np.int64)
# basal_lengths = np.array([r.length() for r in basals])
# basal_ages = np.array(time * len(basals))
#
# ii = np.argsort(basal_lengths)  # sort by ascending lengths
# basal_ids = basal_ids[ii]
# basal_lengths = basal_lengths[ii]
# basal_ages = basal_ages[ii]
#
# print("basal_ids", basal_ids)
#
# # Method 1 (predefine initial growth rate, iteration does not improve result)
# r = 3  # [cm/day] initial
# res, _, ages1 = es.estimate_order0_rate(basal_lengths, r, k, time[0])
# rate1 = res.x[0]
# res, f1 = es.estimate_r(basal_lengths, ages1, k)
# r1 = res.x[0]
# print("prodcution rate", rate1, "growth rate", r1, "err", f1(res.x))
#
# # Method 2 (fit both)
# res, f2, ages2 = es.estimate_order0_rrate(basal_lengths, k, time[0], 1.5)  # third is r0, unstable results
# rate2 = res.x[0]
# r2 = res.x[1]
# print("prodcution rate", rate2, "growth rate", r2, "err", f2(res.x))
#
# t_ = np.linspace(0, time[0], 200)
# y1 = es.Root.negexp_length(t_, r1, k)
# y2 = es.Root.negexp_length(t_, r2, k)
#
# plt.scatter(ages1, basal_lengths, label = "Method 1")
# plt.plot(t_, y1)
# plt.scatter(ages2, basal_lengths, label = "Method 2")
# plt.plot(t_, y2)
# plt.xlabel("estimated root age [day]")
# plt.ylabel("measured root length [cm]")
# plt.legend()
# plt.show()

#
# rate = 3.7, r =  2.93
# basal_ids = [10, 7, 1, 2, 9, 6, 8, 4, 11, 0, 12, 24, 28, 55, 142, 246, 348, 644, 709, 1042, 1298, 1442, 1807, 2203, 2567, 2666, 2768, 2794, 3018, 3383]
#

# from maize_kutschera.py
k0 = 100.  # [cm] fixed
rate0 = 3.7
r0 = 2.93
basal_ids = np.array([10, 7, 1, 2, 9, 6, 8, 4, 11, 0, 12, 24, 28, 55, 142, 246, 348, 644, 709,
                      1042, 1298, 1442, 1807, 2203, 2567, 2666, 2768, 2794, 3018, 3383])

# Set emergance time according to basal production model, see estimate_order0_rate, or estimate_order0_rrate
# need to simplify this ...
ages = np.zeros(basal_ids.shape)
for i, _ in enumerate(basal_ids):
    ages[basal_ids.shape[0] - i - 1] = max(time[0] - i * rate0, 0.)
for i, id in enumerate(basal_ids):
    et = time[0] - ages[i]
    roots[id].set_emergence_time(et)

# set order 0 growth rate and maximal length
for id in basal_ids:
    roots[id].calc_growth_rate(r0, k0)

order1 = es.get_order(1, roots)
print(len(order1), "order 1 roots")
lengths1 = np.array([r.length() for r in order1])  # for one measurement
ages1 = np.array([r.ages[time[-1]] for r in order1])

res, f = es.estimate_rk(lengths1, ages1)
r1 = res.x[0]
k1 = res.x[1]
print("growth rate", r1, "maximal length", k1, "err", f(res.x))
t_ = np.linspace(0, time[0], 200)
y1 = es.Root.negexp_length(t_, r1, k1)

plt.scatter(ages1, lengths1, label = "order 1")
plt.plot(t_, y1, "k")
plt.xlabel("estimated root age [day]")
plt.ylabel("measured root length [cm]")
plt.legend()
plt.show()

print("done")

