import sys;  sys.path.append("../..")
""" determines the (more easy) parameters la, lb, ln, a, theta by order and prints mean and std"""

import numpy as np
import matplotlib.pyplot as plt

import rsml_reader as rsml
import estimate_root_params as es

time = [75]  # measurement times (not in the rsml)
name = ["RSML/Maize_Kutschera.rsml"]
roots = es.parse_rsml(name, time)

orders = np.array([r.order() for r in roots.values()], dtype = np.int64)
max_order = np.max(orders)
print("Number of root orders: ", max_order + 1)

roots_i = []
for i in range(0, max_order + 1):
    roots_i.append(es.get_order(i, roots))

params = []
for i in range(0, 3):
    print("order {:g}".format(i))
    params.append(es.get_params(roots_i[i], time[-1]))
    print(params[i])
    print()

params2 = []
for i in range(1, min(max_order + 1, 4)):
    roots_ = []
    for j in range(i, max_order + 1):
        roots_.extend(roots_i[j])
    print("order >{:g}".format(i - 1))
    params2.append(es.get_params(roots_, time[-1]))
    print(params2[-1])
    print()

print("done")

