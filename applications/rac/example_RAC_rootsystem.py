import sys; sys.path.append("../../..")
from cmath import pi

import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/rootsystem/"
name = "wheat_zero_std"  # "maize_p1_zero_std" #"wheat"
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.setSeed(1)
# Manually set scale elongation function
scale_elongation = pb.EquidistantGrid1D(0, -200, 100)  # for root elongation from 0 cm to -100 cm, 100 layers
scales = np.ones((99,)) * 1  # imedance
scale_elongation.data = scales  # set proportionality factors
for p in rs.getRootRandomParameter():
    p.f_se = scale_elongation
    # p.r = 1;

############################################################################

# Soil core analysis
depth, layers = 100., 20
interrow = 1  # inter-row spacing
row = 1  # row spacing
layerVolume = depth / layers * interrow * row
simtime = 200;
dt = 0.5; N = round(simtime / dt);
z_ = np.linspace(0, -1 * depth, layers)

rs.initialize()
time = np.linspace(0, simtime - 1, N); rld = np.zeros(N);
for i in range(0, N):
    rs.simulate(dt, True)
    ana = pb.SegmentAnalyser(rs)
    rl = ana.distribution("length", 0., -depth, layers, True);
    rl = np.array(rl);
    rld[i] = rl[2] / layerVolume;
rs.write("RAC_unimpeded.vtp")
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 8))
axes.plot(time, rld)
axes.set_title('RAC')

#################################################################################

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
rs.setSeed(1)
# Manually set scale elongation function
scale_elongation = pb.EquidistantGrid1D(0, -200, 100)  # for root elongation from 0 cm to -100 cm, 100 layers
scales = np.ones((99,)) * 0.5  # imedance
scale_elongation.data = scales  # set proportionality factors
for p in rs.getRootRandomParameter():
    p.f_se = scale_elongation
    # p.r = p.r/2;

rs.initialize()
time = np.linspace(0, simtime - 1, N); rld = np.zeros(N);
for i in range(0, round(N)):
    rs.simulate(dt, True)
    ana = pb.SegmentAnalyser(rs)
    rl = ana.distribution("length", 0., -depth, layers, True);
    rl = np.array(rl);
    rld[i] = rl[2] / layerVolume;
rs.write("RAC_impeded.vtp")
axes.plot(time, rld, 'r--')

axes.legend(["unimpeded", "impeded"])
plt.savefig('wheat_heterogeneous scaling.png')

plt.show()
