import sys;  sys.path.append("../..")

import rsml_reader as rsml
import estimate_root_params as es
import numpy as np
import matplotlib.pyplot as plt

time = [140]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

# from maize_kutschera.py
k0 = 100.  # [cm] fixed
rate0 = 6.824382579932482  # [day-1]
r0 = 1.5414324819206926  # [cm/day]
basal_ids = np.array([10, 7, 1, 2, 9, 6, 8, 4, 11, 0, 12, 24, 28, 55, 142, 246, 348, 644, 709, 1042, 1298, 1442, 1807, 2203, 2567, 2666, 2768, 2794, 3018, 3383])

# calculate la, lb, ln, a, theta
for r in roots.items():
    r[1].calc_params()

# need to simplify this ...
ages = np.zeros(basal_ids.shape)
for i, _ in enumerate(basal_ids):
    ages[basal_ids.shape[0] - i - 1] = max(time[0] - i * rate0, 0.)
for i, id in enumerate(basal_ids):
    et = time[0] - ages[i]
    roots[id].set_emergence_time(et)

for id in basal_ids:
    roots[id].calc_growth_rate(r0, k0)

orders = np.array([r[1].order() for r in roots.items()])
print("Number of root orders: ", np.max(orders) + 1)

print("done")

