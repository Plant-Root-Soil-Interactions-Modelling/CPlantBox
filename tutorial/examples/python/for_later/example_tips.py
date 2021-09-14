import py_rootbox as rb

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt

rootsystem = rb.RootSystem()
name = "Zea_mays_2_Pagès_2014"

#
# Open plant and root parameter from a file
#
rootsystem.openFile(name, "modelparameter/")
# rootsystem.writeParameters() # not exposed to python yet

#
# Initialize
#
rootsystem.initialize(4, 5)  # TODO expose default values

#
# Simulate
#
simtime = 60;
rootsystem.simulate(simtime);

#
# Analyse root system
#
tips = rootsystem.getRootTips()
notips = len(tips)

ana = rb.AnalysisSDF(rootsystem)
totalLength = ana.getSummed(rb.ScalarType.length)
l = ana.getScalar(rb.ScalarType.length)
nos = len(l)

print('\nNumber of root tips is ' + str(notips))
print('Number of nodes is ' + str(rootsystem.getNumberOfNodes()))
print('Number of segmentes is ' + str(nos))
print('Length of root system is ' + str(totalLength) + " cm")

rp1 = rootsystem.getRootTypeParameter(1)
rp2 = rootsystem.getRootTypeParameter(2)
rp3 = rootsystem.getRootTypeParameter(3)
rp4 = rootsystem.getRootTypeParameter(4)
print('\nTropisms')
print('#1 Type ' + str(rp1.tropismT) + ', N ' + str(rp1.tropismN) + ', sigma ' + str(rp1.tropismS) + ', dx ' + str(rp1.dx))
print('#2 Type ' + str(rp2.tropismT) + ', N ' + str(rp2.tropismN) + ', sigma ' + str(rp2.tropismS) + ', dx ' + str(rp2.dx))
print('#3 Type ' + str(rp3.tropismT) + ', N ' + str(rp3.tropismN) + ', sigma ' + str(rp3.tropismS) + ', dx ' + str(rp3.dx))
print('#4 Type ' + str(rp4.tropismT) + ', N ' + str(rp4.tropismN) + ', sigma ' + str(rp4.tropismS) + ', dx ' + str(rp4.dx))

#
# Analyse tip distribution
#
tipcoords = np.zeros((notips, 3))
z_ = np.zeros(notips)
l_ = np.zeros(notips)

#
# Copy everything we need
#
c = 0;
for t in tips:
    tipcoords[c, 0] = t.x
    tipcoords[c, 1] = t.y
    tipcoords[c, 2] = t.z
    z_[c] = t.z
    l_[c] = sqrt(t.x * t.x + t.y * t.y)
    c += 1
#
# Figure params
#
lbins = 10
lrange = (0, 10)
zbins = 120
zrange = (-120, 0)

#
# Plot histograms
#
plt.figure(1)
plt.hist(l_, bins = lbins, range = lrange)
plt.title("Root tip radial distance")
plt.show(False)
hl, bins = np.histogram(l_, bins = lbins, range = lrange)

np.savetxt('radialdistribution.txt', hl);  # save results

plt.figure(2)
plt.hist(z_, bins = zbins, range = zrange)
plt.title("Root tip depth")
plt.show()
hz, bins = np.histogram(z_, bins = zbins, range = zrange)

np.savetxt('depthdistribution.txt', hz);  # save results

