import py_rootbox as rb
from rb_tools import *
import sys
import numpy as np
import matplotlib.pyplot as plt
import csv

rs = rb.RootSystem()
#name = "Triticum_aestivum_a_Bingham_2011" # is this the same as your wheat, Shehan?
name = "Zea_mays_1_Leitner_2010"
rs.openFile(name)

# Pore Geometry
x_ = (-10, -5, 1, 15)  # not 0, otherwise we start in crack
y_ = (0, 0, 0, 0)
x_=(-10, -5)
y_=(0,0)
crack = rb.SDF_PlantBox(1.0, 100, 160)  # cm
cracks_ = rb.std_vector_SDF_()
py_cracks = []

for i in range(0, len(y_)):
    v = rb.Vector3d(x_[i], y_[i], 0)
    py_cracks.append(rb.SDF_RotateTranslate(crack, v))
    cracks_.append(py_cracks[-1])

cracks = rb.SDF_Union(cracks_)
rs.setPoreGeometry(cracks)

# Increased elongation within the pores
maxS = 20  # twice the elongation rate within the pore
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
t1 = rb.Vector3d(2, 0, 0)
t2 = rb.Vector3d(0, 0.5, 0)
t3 = rb.Vector3d(0, 0, 0.5)
rs.setPoreConductivity(rb.Matrix3d(t1, t2, t3));

# Set up depth dependent elongation scaling function
scale_elongation = rb.EquidistantGrid1D(0,-100, 11) # todo: replace this by reading in data from CSV file      
scales = np.zeros(len(scale_elongation.grid))+0.1 # scales from some equation (scale = function(soil_strength) ), where scale in (0,1)
scale_elongation.data = a2v(scales) # set proportionality factors
  
# Proportionally scale this function
se2 = rb.ProportionalElongation()
se2.setBaseLookUp(scale_elongation)
  
# multiply the scale elongation functions
se3 = rb.MultiplySoilLookUps(se2,soil_prop)  
  
# Manually set scaling function 
for i in range(0,10):  
    p = rs.getRootTypeParameter(i+1)
    p.se = se3   

# Initialize
rs.initialize()

# Simulate
simtime = 14  # e.g. 30 or 60 days
dt = 1
N = round(simtime / dt)

for i in range(0, N):
    
    # time-dependent and depth-dependent scaling function    
    scales = np.loadtxt('CSV data.csv', delimiter='\t', usecols=i)   # reading in ith column from CSV file
    scale_elongation.data = a2v(scales) # set the data of scale elongation 
    rs.simulate(dt)

# Export results (as vtp)
rs.write("../results/crack.vtp")

# Export cracks
rs.setGeometry(cracks)  # just for vizualisation
rs.write("../results/crack.py")

#
# Compute vertical RLD in layers
#
x = np.linspace(0,100,101);    # one layer is 1 cm deep
nl = len(x)
depth = max(x)
analysis = rb.SegmentAnalyser(rs)
RLD = analysis.distribution(rb.ScalarType.length,0,depth,nl,False)
RLD = np.array(RLD)/(depth/nl)/(75.0*15.0);   # divide by surface area per plant (inter-row spacing times inter-plant spacing)
#
plt.plot(RLD,x*(-1),'r-',linewidth=2)
plt.xlabel('RLD (cm/cm³)')
plt.ylabel('Depth (cm)')
plt.show()

