import sys
from cmath import pi
sys.path.append("../../..")
import plantbox as pb

import numpy as np
import matplotlib.pyplot as plt

path = "../../modelparameter/rootsystem/"

rs = pb.RootSystem()

################### set parameters ##########################################

# Root random parameter
p0 = pb.RootRandomParameter(rs)  # with default values,
p0.name = "taproot"
p0.a = 0.2  # cm radius
p0.subType = 1
p0.lb = 1
p0.la = 10
p0.lmax = 200
p0.ln = 20.
p0.theta = 0. 
p0.r = 2.0  # initial growth rate
p0.dx = 0.5
p0.gf = 2
p0.tropismT = pb.TropismType.gravi
p0.tropismN = 1.
p0.tropismS = 0.0

rs.setOrganRandomParameter(p0)

# Seed random parameter (neglecting shoot borne)
srp = pb.SeedRandomParameter(rs)
srp.seedPos = pb.Vector3d(0., 0., -3.)
srp.maxB = 0
srp.firstB = 100.
srp.delayB = 100.
rs.setRootSystemParameter(srp)

############################################################################

rs.initialize()
rs.simulate(120)
rs.write("RAC.vtp")


# Soil core analysis
depth, layers = 100., 10

z_ = np.linspace(0, -1 * depth, layers)
fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 8))


# Make a root length distribution for periodic soil domain
ana = pb.SegmentAnalyser(rs)
interrow=1 # inter-row spacing
row=1 # row spacing
layerVolume = depth / layers * interrow * row  
N=120*2;
time=np.linspace(0,119,N); rld = np.zeros(N);

for i in range(N-1,0,-1):
    ana.filter("creationTime", 0, time[i])
    rl = ana.distribution("length", 0., -depth, layers, True);
    rl = np.array(rl);
    rld[i] = rl[1]/ layerVolume;
axes.plot(time,rld)
axes.set_title('RAC')

#################################################################################

rs = pb.RootSystem()
############################### parameters

################### set parameters ##########################################

# Root random parameter
p0 = pb.RootRandomParameter(rs)  # with default values,
p0.name = "taproot"
p0.a = 0.2  # cm radius
p0.subType = 1
p0.lb = 1
p0.la = 2
p0.lmax = 200
p0.ln = 0.5
p0.theta = 0. 
p0.r = 2.0  # initial growth rate
p0.dx = 0.5
p0.gf = 2
p0.tropismT = pb.TropismType.gravi
p0.tropismN = 1.
p0.tropismS = 0.0
rs.setOrganRandomParameter(p0)

# Seed random parameter (neglecting shoot borne)
srp = pb.SeedRandomParameter(rs)
srp.seedPos = pb.Vector3d(0., 0., -3.)
srp.maxB = 0
srp.firstB = 100.
srp.delayB = 100.
rs.setRootSystemParameter(srp)

############################################################################

###############################
rs.setOrganRandomParameter(p0)
rs.setRootSystemParameter(srp)
#Manually set scale elongation function
scale_elongation = pb.EquidistantGrid1D(0, -100, 100)  # for root elongation from 0 cm to -100 cm, 100 layers
scales = np.ones((99,)) * 0.5  # imedance
scale_elongation.data = scales  # set proportionality factors
for p in rs.getRootRandomParameter():
    p.f_se =scale_elongation
    #p.r = 1;

rs.initialize()
for i in range(0, round(N)):
    rs.simulate(120/N, True)

# Make a root length distribution for periodic soil domain
ana = pb.SegmentAnalyser(rs)

time=np.linspace(0,119,N); rld = np.zeros(N);
for i in range(N-1,0,-1):
    ana.filter("creationTime", 0, time[i])
    rl = ana.distribution("length", 0., -depth, layers, True);
    rl = np.array(rl);
    rld[i] = rl[1]/ layerVolume;
axes.plot(time,rld,'r--')
axes.legend(["unimpeded", "impeded"])



plt.show()
