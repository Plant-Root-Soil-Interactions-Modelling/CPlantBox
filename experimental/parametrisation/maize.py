import sys
sys.path.append("../..")
import plantbox as pb
import rsml_reader as rsml
import estimate_params_fibrous as es

import math
import numpy as np
import matplotlib.pyplot as plt

fixed_k = 150.  # for the fits based on a literature value

name = "Maize_Kutschera.rsml"

polylines, properties, functions = rsml.read_rsml(name)  # read rsml file

""" recalculate length (should not alter anything) and add order as property"""
es.create_length(polylines, properties)
es.create_order(polylines, properties)

# for i,p in enumerate(polylines):
    # print(properties["parent-node"][i])

# plot root system
# pp=properties["parent-node"]
# rsml.plot_rsml(polylines,pp)
# print(properties["length"][0])

""" find all base roots """
base_polylines, base_properties = es.base_roots(polylines, properties)
pp = base_properties['parent-poly']
# rsml.plot_rsml(base_polylines,pp)

""" find all first order lateral roots """
lateral_polylines = []; lateral_properties = {};
for i, p in enumerate(polylines):
    if properties["order"][i] == 2:
        lateral_polylines.append(p);
        for k in properties.keys():
            try:
                dummy = properties[k][i]
            except IndexError:
                dummy = -99
            if dummy is -99:
                print('Roots with missing property found: ', k)
                lateral_properties.setdefault(k, []).append(dummy)
            else:
                lateral_properties.setdefault(k, []).append(properties[k][i])

# rsml.plot_rsml(lateral_polylines,lateral_properties["parent-poly"])

""" find all first order laterals that are attached to a given polyline """
nob = len(base_polylines);  # ## number of base roots
pp = lateral_properties["parent-poly"]  # parent-poly
upp = list(set(pp))  # unique parent-poly values, number of subsets of laterals with same parent-poly

lateral_polylines_2base = [ [] for k in range(0, len(upp)) ];
lateral_properties_2base = [ [] for k in range(0, len(upp)) ];
# print(len(lateral_polylines_2base),len(upp))   # 21 different base roots and thus clusters of laterals

for i, p in enumerate(lateral_polylines):
    dummy = {};
    ind = upp.index(lateral_properties["parent-poly"][i])  # index of this specific parent polyline in upp
    lateral_polylines_2base[ind].append(p);
    for k in lateral_properties.keys():
        dummy.setdefault(k, []).append(lateral_properties[k][i])
    lateral_properties_2base[ind].append(dummy)

ind = 3;  # Todo: replace with a loop over all clusters
# print(lateral_properties_2base[ind])

prop = lateral_properties_2base[ind]
prop_ = [d['parent-poly'] for d in prop]
flattened = [val for sublist in prop_ for val in sublist]
plt.figure(1)
rsml.plot_rsml(lateral_polylines_2base[ind], flattened)

ppn = flattened[0];  # common parent-poly number of all laterals in a subset with same base root
base = []; prop = {};
for i, p in enumerate(polylines):
    if i == ppn:
        base.append(p)
        for k in properties.keys():
            prop.setdefault(k, []).append(base_properties[k][i])
plt.figure(2)
rsml.plot_rsml(base, prop["length"])

"""this does not work yet
# insert base root at the beginning of the laterals --> this structure can then be used with the ones developed for tap roots with the assumption that first root is the tap root
lateral_polylines_2base[ind].insert(0,base);
lateral_properties_2base[ind].insert(0,prop);


#prop = lateral_properties_2base[ind]
#prop_=[d['parent-poly'] for d in prop]
#flattened = [val for sublist in prop_ for val in sublist]
#rsml.plot_rsml(lateral_polylines_2base[ind],flattened)
"""

print("\nfin.")

