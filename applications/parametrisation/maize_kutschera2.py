import sys;  sys.path.append("../..")
""" determines growth rate and production rate and ages of the basal roots """

import rsml_reader as rsml
import estimate_root_params as es
import numpy as np
import matplotlib.pyplot as plt

lmax = 100.  # [cm] fixed
time = [75]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

basals = es.get_order(0, roots)
basal_ids = np.array([r.id for r in basals], dtype = np.int64)
basal_lengths = np.array([r.length() for r in basals])
basal_ages = np.array(time * len(basals))

ii = np.argsort(basal_lengths)  # sort by ascending lengths
basal_ids = basal_ids[ii]
basal_lengths = basal_lengths[ii]
basal_ages = basal_ages[ii]

print("basal_ids", basal_ids)

# Method 1 (predefine initial growth rate, iteration does not improve result)
r = 3  # [cm/day] initial
res, _, ages1 = es.estimate_order0_rate(basal_lengths, r, lmax, time[0])
rate1 = res.x[0]
res, f1 = es.estimate_r(basal_lengths, ages1, lmax)
r1 = res.x[0]
r1s, _ = es.get_see_rk(basal_lengths, ages1, r1, lmax)
print("prodcution rate", rate1, "growth rate", r1, "see", r1s, "err", f1(res.x))

# Method 2 (fit both)
res, f2, ages2 = es.estimate_order0_rrate(basal_lengths, lmax, time[0], 1.5)  # third is r0, unstable results
rate2, r2 = res.x[0], res.x[1]
r2s, _ = es.get_see_rk(basal_lengths, ages1, r2, lmax)
print("prodcution rate", rate2, "growth rate", r2, "see", r2s, "err", f2(res.x))

t_ = np.linspace(0, time[0], 200)
y1 = es.Root.negexp_length(t_, r1, lmax)
y2 = es.Root.negexp_length(t_, r2, lmax)

plt.scatter(ages1, basal_lengths, label = "Method 1")
plt.plot(t_, y1)
plt.scatter(ages2, basal_lengths, label = "Method 2")
plt.plot(t_, y2)
plt.xlabel("estimated root age [day]")
plt.ylabel("measured root length [cm]")
plt.legend()
plt.show()

rate0, r0 = rate2, r2  # or rate2, r2

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
    roots[id].calc_growth_rate(r0, lmax)

order1 = es.get_order(1, roots)
print(len(order1), "order 1 roots")
lengths1 = np.array([r.length() for r in order1])  # for one measurement
ages1 = np.array([r.ages[time[-1]] for r in order1])

res, f = es.estimate_rk(lengths1, ages1)
r1, lmax1 = res.x[0], res.x[1]
r1s, lmax1s = es.get_see_rk(lengths1, ages1, r1, lmax1)

# lmax1 = 5
# res, f = es.estimate_r(lengths1, ages1, lmax1)
# r1 = res.x[0]
# r1s, lmaxs1 = es.get_see_rk(lengths1, ages1, r1, lmax1)

print("growth rate", r1, "see", r1s, "maximal length", lmax1, "see", lmax1s, "err", f(res.x))
t_ = np.linspace(0, time[0], 200)
y1 = es.Root.negexp_length(t_, r1, lmax1)

plt.scatter(ages1, lengths1, label = "order 1")
plt.plot(t_, y1, "k")
plt.xlabel("estimated root age [day]")
plt.ylabel("measured root length [cm]")
plt.legend()
plt.show()

print("done")

