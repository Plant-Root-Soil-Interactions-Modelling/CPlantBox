"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import rsml.rsml_reader as rsml

import numpy as np
import matplotlib.pyplot as plt


polylines, props, funcs, metadata = rsml.read_rsml("Wurzelatlas_Tomate.rsml")
print(len(polylines), "roots")
# print(props["parent-poly"])
# print(props["parent-node"])
print("Taproot number of nodes", len(polylines[0]))
firstOrderLengths = []
firstOrderIndices = []
parentNodeIndices = []
childNodeIndices = []
secondOrderLengths = []
    
for i in range(1, len(polylines)-1):
    pn = int(props["parent-node"][i])
    try:
        parentNodeIndices.append(pn)
        childNodeIndices.append(i)
        if pn == 1:
            firstOrderLengths.append(props["length"][i])
            firstOrderIndices.append(i)
        elif pn > 1 and pn in firstOrderIndices:
            secondOrderLengths.append(props["length"][i])
    except:
        print("root", i, "index out of range", "parent-node", pn, "length", props["length"][i])

plt.hist(parentNodeIndices, bins=200)
plt.show()

parentNodeIndices = np.unique(parentNodeIndices)
childNodeIndices = np.unique(childNodeIndices)


NotTapRoots = [int(x) for x in childNodeIndices if x in parentNodeIndices]
NotTapRoots = np.unique(NotTapRoots)

print("First order laterals:", firstOrderLengths)
print("Second order laterals:", secondOrderLengths)
print("parent node indices:", parentNodeIndices)
# print("child node indices:", childNodeIndices)
print("Not tap roots (both first and second order laterals):", NotTapRoots)


print("Number of parent nodes:", len(parentNodeIndices), "Number of child nodes that are also parent nodes:", len(NotTapRoots))