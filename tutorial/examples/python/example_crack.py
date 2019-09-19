import py_rootbox as rb
import sys
import numpy as np
import matplotlib.pyplot as plt

rs = rb.RootSystem()
# name = "Triticum_aestivum_a_Bingham_2011" # is this the same as your wheat, Shehan?
name = "Zea_mays_1_Leitner_2010"
rs.openFile(name)

# Pore Geometry
x_ = (-10, -5, 1, 15)  # not 0, otherwise we start in crack
y_ = (0, 0, 0, 0)
crack = rb.SDF_PlantBox(0.25, 100, 160)  # cm
cracks_ = rb.std_vector_SDF_()
py_cracks = []

for i in range(0, len(y_)):
    v = rb.Vector3d(x_[i], y_[i], 0)
    py_cracks.append(rb.SDF_RotateTranslate(crack, v))
    cracks_.append(py_cracks[-1])

cracks = rb.SDF_Union(cracks_)
rs.setPoreGeometry(cracks)

# Increased elongation within the pores
maxS = 2  # twice the elongation rate within the pore
minS = 1  # normal elongation rate
slope = 0
soil_prop = rb.SoilLookUpSDF(cracks, maxS, minS, slope)

# Adjust Tropism
sigma = [0.4] * 10
for i in range(0, 10):
    p = rs.getRootTypeParameter(i + 1)
    p.dx = 0.25  # adjust resolution
    p.tropismT = rb.TropismType.gravi
    p.tropismN = 1  # strength of tropism
    p.tropismS = sigma[i]
    p.se = soil_prop

# Pore Local Axes
v1 = rb.Vector3d(0, 0, -1)
v2 = rb.Vector3d(1, 0, 0)
v3 = rb.Vector3d(0, 1, 0)
rs.setPoreLocalAxes(rb.Matrix3d(v1, v2, v3));

# Pore Conductivity Tensor
t1 = rb.Vector3d(1.5, 0, 0)
t2 = rb.Vector3d(0, 0.5, 0)
t3 = rb.Vector3d(0, 0, 0.5)
rs.setPoreConductivity(rb.Matrix3d(t1, t2, t3));

# Initialize
rs.initialize()

# Simulate
simtime = 240  # e.g. 30 or 60 days
dt = 1
N = round(simtime / dt)

for _ in range(0, N):
    rs.simulate(dt)

# Export results (as vtp)
rs.write("../results/crack.vtp")

# Export cracks
rs.setGeometry(cracks)  # just for vizualisation
rs.write("../results/crack.py")

# Make a root length distribution
# ana = rb.SegmentAnalyser(rs)
# rl_ = ana.distribution(rb.ScalarType.length, 0, depth, layers, True)
# np.set_printoptions(precision = 4)
# np.savetxt("c_" + str(cc) + ".txt", rl_, fmt = "%.2f")

