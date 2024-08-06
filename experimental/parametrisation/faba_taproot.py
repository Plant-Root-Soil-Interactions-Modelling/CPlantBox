import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import rsml.rsml_reader as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt

fixed_k = 150.  # for the fits based on a literature value

times = [7, 11, 15]  # measurement times (not in the rsml)
names = ["Faba_12/DAP15.rsml", "Faba_14/DAP15.rsml", "Faba_16/DAP15.rsml", "Faba_20/DAP15.rsml", "Faba_24/DAP15.rsml" ]

""" open all files, and split into measurent times"""
polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
properties = [ [ [] for t in times] for i in range(0, len(names)) ]
functions = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    print(names[i])
    j = len(times) - 1
    p, properties[i][j], functions[i][j], _ = rsml.read_rsml("RSML/" + names[i])  # read file to final time step
    print(properties[i][j].keys())
    print(functions[i][j].keys())
    p = [[np.array([p[i][j][k] / 10 for k in range(0, 3)]) for j in range(0, len(p[i])) ] for i in range(0, len(p)) ]  # mm -> cm
    polylines[i][j] = p
    for k in range(0, j):  # truncate the others
        polylines[i][k], properties[i][k], functions[i][k] = es.measurement_time(polylines[i][j], properties[i][j], functions[i][j], times[k])
# polylines[i][j] contains i-th plant, j-th measurement time

""" find all base roots """
base_polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
base_properties = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    for j in range(0, len(times)):
            base_polylines[i][j], base_properties[i][j] = es.base_roots(polylines[i][j], properties[i][j])

# for Faba this is exactly 1 for all 5 measurements

""" recalculate length (should not alter anything)"""
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_length(base_polylines[i][j], base_properties[i][j])

""" fit data """
print()
l = np.array([[base_properties[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])

r0, f0 = es.fit_taproot_r(l[0:1,:], [times[0]], fixed_k)
k0 = fixed_k
print("k fixed, first ", r0, "cm/day", k0, "cm; error", f0)

r1, f1 = es.fit_taproot_r(l, times, fixed_k)
k1 = fixed_k
print("k fixed        ", r1, "cm/day", k1, "cm; error", f1)

r2, k2, f2 = es.fit_taproot_rk(l, times)
print("fit r, k       ", r2, "cm/day", k2, "cm; error", f2)

""" length plot """
t_ = np.linspace(0, times[-1], 200)
y0 = es.negexp_length(t_, r0, k0)
y1 = es.negexp_length(t_, r1, k1)
y2 = es.negexp_length(t_, r2, k2)
plt.plot([0.], [0.], "r*")  # we can add that point
for i in range(0, len(names)):
    for j in range(0, len(times)):
        l = base_properties[i][j]["length"]
        plt.plot([times[j] for k in l], l, "k*")
plt.plot(t_, y0, "b", label = "first only, k fixed")
plt.plot(t_, y1, "g", label = "k fixed")
plt.plot(t_, y2, "r", label = "fit r, k")

plt.legend()
plt.title("Faba")

plt.xlabel("Time [days]")
plt.ylabel("Length [cm]")
plt.show()

print("\nfin.")

