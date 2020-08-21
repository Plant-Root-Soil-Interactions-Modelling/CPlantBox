import sys
sys.path.append("../..")
import plantbox as pb
import rsml_reader_px as rsml
import estimate_params as es

import math
import numpy as np
import matplotlib.pyplot as plt

fixed_k = 100.  # for the fits based on a literature value

times = [140]  # measurement times (not in the rsml)
names = ["Maize_Kutschera.rsml"]

""" open all files, and split into measurent times"""
polylines = [ [ [] for t in times] for i in range(0, len(names)) ] 
properties = [ [ [] for t in times] for i in range(0, len(names)) ] 
functions = [ [ [] for t in times] for i in range(0, len(names)) ]    
for i in range(0, len(names)):
    print(names[i])   
    j = len(times) - 1     
    p, properties[i][j], functions[i][j] = rsml.read_rsml(names[i])  # read file to final time step
    p = [[np.array([p[i][j][k] for k in range(0, 3)]) for j in range(0, len(p[i])) ] for i in range(0, len(p)) ]  
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


""" recalculate length (should not alter anything)"""
for i in range(0, len(names)):
    for j in range(0, len(times)):
            es.create_length(base_polylines[i][j], base_properties[i][j]) 

""" fit data """
for i in range(0, len(names)):
    for j in range(0, len(times)):   
        l = base_properties[i][j]["length"]
print("maximum root length: ",np.amax(l))




print("\nfin.")

