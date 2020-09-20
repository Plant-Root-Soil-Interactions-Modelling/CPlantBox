import sys;  sys.path.append("../..")
""" determines the (more easy) parameters la, lb, ln, a, theta by order and prints mean and std"""

import numpy as np
import matplotlib.pyplot as plt

import rsml_reader as rsml
import estimate_root_params as es
from xylem_flux import XylemFluxPython
import vtk_plot as vp

time = [1, 2, 3, 4, 5, 6, 8, 9]  # measurement times (not in the rsml)
name = ["RSML/m1/dicot/lupin/lupin_d{:g}.rsml".format(a) for a in time]

# time = [75]  # measurement times (not in the rsml)
# name = ["RSML/Maize_Kutschera.rsml"]

roots = es.parse_rsml(name, time)

orders = np.array([r.order() for r in roots.values()], dtype = np.int64)
max_order = np.max(orders)
print("Number of root orders: ", max_order + 1)

roots_i = []
for i in range(0, max_order + 1):
    roots_i.append(es.get_order(i, roots))

params = []
for i in range(0, max_order + 1):
    print("order {:g}".format(i))
    params.append(es.get_params(roots_i[i], time[-1]))
    print(params[i])
    print()

base_roots = es.get_order(0, roots)
for b in base_roots:
    b.set_emergence_time(0.)
base_lengths = np.array([r.length() for r in base_roots])
base_ages = np.array(time * len(base_roots))

print("done")

