"""analysis of results/ using signed distance functions"""
import sys
from cmath import pi
sys.path.append("../..")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt
#from vtk_tools import *

path = "../../modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"  # ""

rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")

# Create and set geometry
rs.setMinDx(1.e-3)
x0 = pb.Vector3d(0., 0., -1.)
nx = pb.Vector3d(1., 0., -1.)
ny = pb.Vector3d(0., 1., -1.)
soil_layer = pb.SDF_HalfPlane(x0, nx, ny)  # there was bug, with updated CPlantBox
rs.setGeometry(soil_layer)

rs.setSeed(2)
rs.initialize()

rs.simulate(7)
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.2, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/" + name + "/" + name +  "_7days.vtp")
ana.write("results/" + name + "/" + name + "_7days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(7)
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.2, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/" + name + "/" + name +  "_14days.vtp")
ana.write("results/" + name + "/" + name + "_14days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))

rs.simulate(140)
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.2, True)  # ct, radius, insert first
    vt = rs.getSummed("volume")
print("Volume (cm3)", vt)
ana.write("results/" + name + "/" + name +  "_154days.vtp")
ana.write("results/" + name + "/" + name + "_154days.dgf")

l = np.array(ana.getParameter("length"))
print("Min length", np.min(l))
a = np.array(ana.getParameter("radius"))
print("Min radius", np.min(a))
