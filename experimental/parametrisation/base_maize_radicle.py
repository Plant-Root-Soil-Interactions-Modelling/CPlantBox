import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import rsml.rsml_reader as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt

# times = [7, 11, 15]  # not in the rsml
# names = ["Faba_12/DAP15.rsml", "Faba_14/DAP15.rsml", "Faba_16/DAP15.rsml", "Faba_20/DAP15.rsml", "Faba_24/DAP15.rsml" ]
times = [11, 15, 19]
names = ["Maize2/DAP19.rsml", "Maize4/DAP19.rsml", "Maize6/DAP19.rsml", "Maize8/DAP19.rsml", "Maize10/DAP19.rsml" ]

""" open all files, and split into measurent times"""
polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
properties = [ [ [] for t in times] for i in range(0, len(names)) ]
functions = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    print(names[i])
    j = len(times) - 1
    p, properties[i][j], functions[i][j], _ = rsml.read_rsml("RSML/" + names[i])  # read file to final time step
    p = [[np.array([p[i][j][k] / 10 for k in range(0, 3)]) for j in range(0, len(p[i])) ] for i in range(0, len(p)) ]  # mm -> cm
    polylines[i][j] = p
    es.create_order(polylines[i][j], properties[i][j])  # add root order
    # es.create_length(polylines[i][j], properties[i][j])  # add root order
    for k in range(0, j):  # truncate the others
        polylines[i][k], properties[i][k], functions[i][k] = es.measurement_time(polylines[i][j], properties[i][j], functions[i][j], times[k])
# polylines[i][j] contains i-th plant, j-th measurement time

""" find all base roots """
base_polylines = [ [ [] for t in times] for i in range(0, len(names)) ]
base_properties = [ [ [] for t in times] for i in range(0, len(names)) ]
for i in range(0, len(names)):
    for j in range(0, len(times)):
            base_polylines[i][j], base_properties[i][j] = es.base_roots(polylines[i][j], properties[i][j])

""" recalculate length """
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_length(base_polylines[i][j], base_properties[i][j])

""" plots """
vis_prop_name = "type"

for mi in range(0, len(names)):
    base_vis_p = [[] for t in times]
    for i in range(0, len(times)):
        base_vis_p[i] = base_properties[0][i][vis_prop_name]  # so many indices...
    rsml.plot_multiple_rsml([base_polylines[mi][0], base_polylines[mi][1], base_polylines[mi][2]], base_vis_p, times)  # looks good
plt.show()

l = np.array([[base_properties[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])
# l_tap = np.array([[l[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])
print(l)

k0 = 50.
r0, f0 = es.fit_taproot_r(l[0:1,:], [times[0]], k0)
print(r0, "cm/day", k0, "cm")

k1 = 50.
r1, f1 = es.fit_taproot_r(l, times, k1)
print(r1, "cm/day", k1, "cm")

r2, k2, f2 = es.fit_taproot_rk(l, times)
print(r2, "cm/day", k2, "cm")

t_ = np.linspace(0, times[-1], 200)
y0 = es.negexp_length(t_, r0, k0)
y1 = es.negexp_length(t_, r1, k1)
y2 = es.negexp_length(t_, r2, k2)

""" length plot """
c = ["r*", "g*", "b*", "m*", "c*"]
plt.plot([0.], [0.], "r*")  # we can add that point
for i in range(0, len(names)):
    for j in range(0, len(times)):
        # print(len(times), l.shape, i, j)
        plt.plot(times[j], l[j, i], c[i])

plt.plot(t_, y0, "b", label = "first only, k fixed")
plt.plot(t_, y1, "g", label = "fixed")
plt.plot(t_, y2, "r", label = "fit r, k")
plt.legend()
plt.title("Maize radicle")

plt.xlabel("Time [days]")
plt.ylabel("Length [cm]")
plt.show()

print("fin.")

