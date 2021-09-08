""" analysis and plots a root system """
import sys; sys.path.append("../../.."); sys.path.append("../../../src/python_modules")
import numpy as np

import plantbox as pb
import vtk_plot as vp
from xylem_flux import XylemFluxPython

path = "../../../modelparameter/rootsystem/"  # to cplantbox parameter folder
name = "Moraes2018_dry_2008"  # of xml cplantbox parameter file

simtime = 154

print("\n1. Initialize")
rs = pb.MappedRootSystem()
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

print("\n2. Simulate")
rs.simulate(simtime, True)

print("\n3. SegmentAnalyser and artificial shoot")
ana = pb.SegmentAnalyser(rs)
aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)
ana.write("results/Moraes2018_{:s}days.dgf".format(str(simtime)))
# r = XylemFluxPython(rs)  # just for test
# r.test()

nodes = ana.nodes
segs = ana.segments
print("Number of nodes", len(nodes))
print("Number of segments", len(segs))
minx = np.min([n.x for n in nodes])
maxx = np.max([n.x for n in nodes])
miny = np.min([n.y for n in nodes])
maxy = np.max([n.y for n in nodes])
minz = np.min([n.z for n in nodes])
maxz = np.max([n.z for n in nodes])
print("Bounding box", minx, maxx, miny, maxy, minz, maxz)

# Plot, using vtk
print("\n4. Plot")
vp.plot_roots(rs, "subType")
# vp.plot_roots(ana, "subType")  # DOES NOT WORK (TODO)
