import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import rsml.rsml_reader as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt
import copy

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

col = ["r*", "g*", "b*", "m*", "c*"]
l_tap = np.array([[base_properties[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])
l = np.array([[base_properties[i][j]["length"][1:] for i in range(0, len(names))] for j in range(0, len(times)) ])  # [1:] without the primary root
number_of_l = np.array([[len(l[j, i]) for i in range(0, len(names))] for j in range(0, len(times)) ])

""" 
Plots 
"""

# """ length plot """
# plt.plot([0.], [0.], "r*")  # we can add that point
# c = 0
# for i in range(0, len(names)):
#     for j in range(0, len(times)):
#         l[j, i].sort()
#         plt.plot(list(range(c, c + len(l[j, i]))), l[j, i], col[i])
#     c += len(l[-1, i])
# plt.ylabel("Root length [cm]")
# plt.title("Maize seminals")

# """ number bar plot """
# width = 0.2  # the width of the bars
# x = np.arange(len(names))
# fig, ax = plt.subplots()
# for i, t in enumerate(times):
#     ax.bar(x - 2 * width + i * width, number_of_l[i], width)
# ax.set_xticks(x)
# ax.set_xticklabels(names)
# plt.ylabel("Number of seminal roots [1]")

""" estimate delay of seminals """
initial_number = 2.  # we start with to seminals (need to check literature)
delay = es.fit_number_of_roots(times, number_of_l, initial_number)
print("Seminal delay", delay, "days")
delay_ = []
for i in range(0, len(names)):
    delay_.append(es.fit_number_of_roots(times, number_of_l[:, i:i + 1], initial_number))
    print("individual delay", delay_[-1], "days")

# """ number plot """
# width = 0.2  # the width of the bars
# x = np.arange(len(names))
# plt.plot(0., initial_number, "k*")
# plt.ylabel("Number of seminal roots [1]")
# plt.xlabel("Time [days]")
# t_ = np.linspace(0, times[-1], 400)
# y = initial_number + np.round(t_ / delay)
# plt.plot(t_, y, "k")
# for i in range(0, len(names)):
#     y = initial_number + np.round(t_ / delay_[i])
#     plt.plot(t_, y, col[i][0] + ":")
# for j, n in enumerate(names):  # PLOT DATA
#     for i, t in enumerate(times):
#         plt.plot(t, number_of_l[i][j], col[j])
# plt.ylim(0, 20)

for i in range(0, len(names)):
    for j in range(0, len(times)):
        l[j, i].sort(reverse = True)

age = copy.deepcopy(l)
times0 = [0]
times0.extend(times)
for i in range(0, len(names)):
    for j in range(0, len(times)):

        dt = times0[j + 1] - times0[j]  # time span
        if j == 0:
            non = number_of_l[j, i] - initial_number  # number of new roots
        else:
            non = number_of_l[j, i] - number_of_l[j - 1, i]  # number of new roots
        c = 0
        for k, l_ in enumerate(l[j, i]):

#             """ method 1 using linear model """
#             et = max((k - initial_number) * delay_[i], 0)  # or delay instead of delay_
#             age[j, i][k] = max(times[j] - et, 0)

            """ data only, linearly interpolate between data """
            if j == 0:
                if k <= initial_number:
                    age[j, i][k] = times0[1]
                else:
                    age[j, i][k] = times0[1] - (k - initial_number) * (dt / non)
            else:
                if k < len(age[j - 1, i]):  # old root
                    age[j, i][k] = age[j - 1, i][k] + dt
                else:  # new root
                    age[j, i][k] = times0[j + 1] - c * (dt / non)
                    c += 1

# """ modeled length plot """
# for i in range(0, len(names)):
#     for j in [0]:  # range(0, len(times)) PICK time step 0, 1, 2 (all are too confusing)
#         for k, l_ in enumerate(l[j, i]):
#             et = max((k - initial_number) * delay_[i], 0)  # emergence time of the root
#             plt.plot(et, l[j, i][k], col[i])
# plt.xlabel("Modelled emergence time [day]")
# plt.ylabel("Root length at measurement time [cm]")
# plt.title("Maize seminals")

""" root age plot """
for j in range(0, len(times)):
    age_flat, l_flat, c_flat = [], [], []
    for i in range(0, len(names)):
        for k, l_ in enumerate(l[j, i]):
            age_flat.append(age[j, i][k])
            l_flat.append(l[j, i][k])
    plt.plot(age_flat, l_flat, col[j])

age_flat = [item for sublist in age for sublist2 in sublist for item in sublist2]
l_flat = [item for sublist in l for sublist2 in sublist for item in sublist2]
r, k, f = es.fit_taproot_rk(l_flat, age_flat)
k1 = 150.
r1, f1 = es.fit_taproot_r(l_flat, age_flat, k1)
print(r, "cm/day", k, "cm")
print(r1, "cm/day", k1, "cm")
t_ = np.linspace(0, times[-1], 200)
y0 = es.negexp_length(t_, r, k)
y1 = es.negexp_length(t_, r1, k1)
plt.plot(t_, y0, "k")
plt.plot(t_, y1, "k:")
plt.xlabel("root age [day]")
plt.ylabel("root length [cm]")

plt.show()
print("fin.")

