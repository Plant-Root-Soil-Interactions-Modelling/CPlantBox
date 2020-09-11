import sys;  sys.path.append("../..")
""" determines growth rate and production rate and ages of the basal roots """

import numpy as np
import matplotlib.pyplot as plt

import rsml_reader as rsml
import estimate_root_params as es

time = [75]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

basal_ids = np.array(es.get_order(0, roots), dtype = np.int64)
basal_lengths = np.array([roots[i].length() for i in basal_ids])
basal_ages = np.array(time * len(basal_ids))

ii = np.argsort(basal_lengths)  # sort by ascending lengths
basal_lengths = basal_lengths[ii]
basal_ids = basal_ids[ii]
basal_ages = basal_ages[ii]

print("basal_ids", basal_ids)
k = 100.  # [cm] fixed

# Method 1 (predefine initial growth rate, iteration does not improve result)
r = 3  # [cm/day] initial
res, _, ages1 = es.estimate_order0_rate(basal_lengths, r, k, time[0])
rate1 = res.x[0]
res, f1 = es.estimate_r(basal_lengths, ages1, k)
r1 = res.x[0]
print("prodcution rate", rate1, "growth rate", r1, "err", f1(res.x))

# Method 2 (fit both)
res, f2, ages2 = es.estimate_order0_rrate(basal_lengths, k, time[0], 1.)  # third is r0, unstable results
rate2 = res.x[0]
r2 = res.x[1]
print("prodcution rate", rate2, "growth rate", r2, "err", f2(res.x))

t_ = np.linspace(0, time[0], 200)
y1 = es.negexp_length(t_, r1, k)
y2 = es.negexp_length(t_, r2, k)

plt.scatter(ages1, basal_lengths, label = "Method 1")
plt.plot(t_, y1)
plt.scatter(ages2, basal_lengths, label = "Method 2")
plt.plot(t_, y2)
plt.xlabel("estimated root age [day]")
plt.ylabel("measured root length [cm]")
plt.legend()
plt.show()

#
# rate = 6.824382579932482, r =  1.5414324819206926
# basal_ids = [10, 7, 1, 2, 9, 6, 8, 4, 11, 0, 12, 24, 28, 55, 142, 246, 348, 644, 709, 1042, 1298, 1442, 1807, 2203, 2567, 2666, 2768, 2794, 3018, 3383]
#

print("done")

