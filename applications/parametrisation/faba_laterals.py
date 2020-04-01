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
    # es.create_order(polylines[i][j], properties[i][j])  # add root order   
    for k in range(0, j):  # truncate the others
        polylines[i][k], properties[i][k] = es.measurement_time(polylines[i][j], properties[i][j], functions[i][j], times[k])
# polylines[i][j] contains i-th plant, j-th measurement time

""" find all base roots """
base_polylines = [ [ [] for t in times] for i in range(0, len(names)) ] 
base_properties = [ [ [] for t in times] for i in range(0, len(names)) ] 
for i in range(0, len(names)):
    for j in range(0, len(times)):
            base_polylines[i][j], base_properties[i][j] = es.base_roots(polylines[i][j], properties[i][j])
            
for i in range(0, len(names)):
    print("Number of nodes", len(base_polylines[i][-1][0]))  # approx 1 node per mm
# for Faba this is exactly 1 for all 5 measurements

""" recalculate length """
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_length(polylines[i][j], properties[i][j])  # add root length
            es.create_length(base_polylines[i][j], base_properties[i][j])

""" reconstruct laterals """
for i in range(0, len(names)):
    for j in range(0, len(times)):            
            print("Measurement", i, "number of roots in time", j, "is", len(polylines[i][j]))
            polylines[i][j], properties[i][j] = es.reconstruct_laterals(polylines[i][j], properties[i][j], base_polylines[i][j][0]) 
            es.create_order(polylines[i][j], properties[i][j])  # add root order  

# """ plots """
mi = 2
vis_p = [[] for t in times]
for i in range(0, len(times)):
    vis_p[i] = properties[mi][i]["order"]  
rsml.plot_multiple_rsml([polylines[mi][0], polylines[mi][1], polylines[mi][2]], vis_p, times) 

""" measure tap root parameters """
la_ = []
lb_ = []
ln_ = []
for i in range(0, len(names)):
    lb, ln, la0 = es. analyze_zones(polylines[i][0], properties[i][0])    
    la_.append(la0)    
    for j in range(0, len(times)):            
            lb, ln, la = es. analyze_zones(polylines[i][j], properties[i][j])
            lb_.append(lb)
            ln_.extend(ln)
            k = la + lb + np.sum(np.array(ln))            
            k2 = es.polyline_length(0, len(polylines[i][j][0]) - 1, polylines[i][j][0])
            print("Measurement", i, "time", j, "lb", lb, "la", la,
                  "banching zone", np.sum(np.array(ln)), "k", k, k2, "length", properties[i][j]["length"][0])

r = 2.0707344055175803  # cm/day , see faba_taproot (blue fit)
lmax = 150  # need literature value             
print("r", r, "cm/day")
print("lmax", lmax, "cm")
print("la", np.mean(la_), "cm", np.std(la_), "cm")            
print("lb", np.mean(lb_), "cm", np.std(lb_), "cm")
print("ln", np.mean(ln), "cm", np.std(ln), "cm")
print("maximal number of branches", round((lmax - np.mean(la_) - np.mean(lb_)) / np.mean(ln)))

print("delay",)
# delay = (age(k) - age(k - la)) - final_lateral / r  

""" estimate lateral age """

# """ root age plot """ 
# for j in range(0, len(times)):             
#     age_flat, l_flat, c_flat = [], [], []
#     for i in range(0, len(names)):    
#         for k, l_ in enumerate(l[j, i]):   
#             age_flat.append(age[j, i][k])
#             l_flat.append(l[j, i][k])   
#     plt.plot(age_flat, l_flat, col[j])
#    
# age_flat = [item for sublist in age for sublist2 in sublist for item in sublist2]
# l_flat = [item for sublist in l for sublist2 in sublist for item in sublist2]       
# r, k = es.fit_taproot_rk(l_flat, age_flat)
# k1 = 150.
# r1 = es.fit_taproot_r(l_flat, age_flat, k1)
# print(r, "cm/day", k, "cm") 
# print(r1, "cm/day", k1, "cm")        
# t_ = np.linspace(0, times[-1], 200)
# y0 = es.negexp_growth(t_, r, k)
# y1 = es.negexp_growth(t_, r1, k1)
# plt.plot(t_, y0, "k")            
# plt.plot(t_, y1, "k:")            
# plt.xlabel("root age [day]")    
# plt.ylabel("root length [cm]")

print("fin.")

