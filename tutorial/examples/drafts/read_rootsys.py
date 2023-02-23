import numpy as np
from io import StringIO  # StringIO behaves like a file object

fname = "RootSys1"

with open(fname) as f:
    content = f.readlines()

# read relevant data

table1 = ""  # holds: segID#    x          y          z      prev or  br#  length   surface  mass
table2 = ""  # holds: origination time

i = 28  # start row
while i < len(content):
    line = content[i]
    if len(line) < 40:  # lousy stopping criteria
        break
    table1 += (line + "\n")
    table2 += (content[i + 1] + "\n")
    i += 2

id, x, y, z, prev, order, bn, length, surface , mass = np.loadtxt(StringIO(table1), unpack = True)
origination, time = np.loadtxt(StringIO(table2), unpack = True)
