import sys
sys.path.append("../..")
import plantbox as pb
import rsml_reader as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt

times = [7, 11, 15]  # not in the rsml
names = ["Faba_12/DAP15.rsml", "Faba_14/DAP15.rsml", "Faba_16/DAP15.rsml", "Faba_20/DAP15.rsml", "Faba_24/DAP15.rsml" ]

""" open all files, and split into measurent times"""
polylines = [ [ [] for t in times] for i in range(0, len(names)) ] 
properties = [ [ [] for t in times] for i in range(0, len(names)) ] 
functions = [ [ [] for t in times] for i in range(0, len(names)) ]    
for i in range(0, len(names)):
    print(names[i])   
    j = len(times) - 1     
    p, properties[i][j], functions[i][j] = rsml.read_rsml("RSML/" + names[i])  # read file to final time step
    p = [[np.array([p[i][j][k] / 10 for k in range(0, 3)]) for j in range(0, len(p[i])) ] for i in range(0, len(p)) ]  # mm -> cm
    polylines[i][j] = p
    es.create_order(polylines[i][j], properties[i][j])  # add root order
    # es.create_length(polylines[i][j], properties[i][j])  # add root order    
    for k in range(0, j):  # truncate the others
        polylines[i][k], properties[i][k] = es.measurement_time(polylines[i][j], properties[i][j], functions[i][j], times[k])
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
mi = 3
base_vis_p = [[] for t in times]
for i in range(0, len(times)):
    base_vis_p[i] = base_properties[mi][i][vis_prop_name]  # so many indices...

vis_p = [[] for t in times]
for i in range(0, len(times)):
    vis_p[i] = properties[mi][i][vis_prop_name]  # so many indices...

for i in range(0, len(times)):
    print("number of roots in time", i, "is", len(polylines[mi][i]))
rsml.plot_multiple_rsml([polylines[mi][0], polylines[mi][1], polylines[mi][2]], vis_p, times)  # looks good

l = np.array([[base_properties[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])
# l_tap = np.array([[l[i][j]["length"][0] for i in range(0, len(names))] for j in range(0, len(times)) ])
print(l)

plt.show()  

print("fin.")

