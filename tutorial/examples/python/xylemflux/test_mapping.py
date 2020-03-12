import sys
sys.path.append("../../../..")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt

min_ = np.array([-3, -3, -10])
max_ = np.array([6, 2, -5])
res_ = np.array([3, 2, 5])


def vector_3d(a):
    return pb.Vector3d(a[0], a[1], a[2])

""" 
Tests if MappedSegments work 
"""
fig, (ax1, ax2) = plt.subplots(1, 2)

""" Roots """
r = XylemFluxPython("RootSystem.rsml")  # returns a MappedSegments object
r.rs.setRectangularGrid(vector_3d(min_), vector_3d(max_), vector_3d(res_))  # cut and map segments

# Assign soil cell index (per hand todo)
segs = r.rs.segments
nodes = r.rs.nodes
x = np.zeros(len(segs))
for c, s in enumerate(segs):
    i, j = s.x, s.y
    n1 = np.array(nodes[i])
    n2 = np.array(nodes[j])
    mid = 0.5 * (n1 + n2)
    x[c] = r.rs.soil_index(mid[0], mid[1], mid[2])

y = np.zeros(x.shape[0])
for i in range(0, x.shape[0]):
   y[i] = 1 if x[i] >= 0 else -1  #  i % 2
y[0] = -2

for i in range(0, x.shape[0]):
   x[i] = x[i] if x[i] >= 0 else -1

ana = pb.SegmentAnalyser(r.rs)
ana.addData("linear_index", x)  # node data are converted to segment data
ana.addData("in", y)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime", "linear_index", "in"])
rootActor, scalarBar = vp.plot_roots(pd, "in", False)

""" Mesh """
width_ = max_ - min_
ind_ = res_ / width_
print("min    ", min_)
print("max    ", max_)
print("width  ", width_)
print("cuboids", 1 / ind_)
grid = vp.uniform_grid(min_, max_, res_)
meshActor = vp.plot_mesh_wireframe(grid, "", False)

vp.render_window([rootActor, meshActor], "mixed fun", scalarBar)

