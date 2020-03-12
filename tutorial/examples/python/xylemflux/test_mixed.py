import sys
sys.path.append("../../../..")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt


def vector_3d(a):
    return pb.Vector3d(a[0], a[1], a[2])

""" 
Tests if MappedSegments work 
"""
fig, (ax1, ax2) = plt.subplots(1, 2)

""" root problem """
r = XylemFluxPython("RootSystem.rsml")  # returns a MappedSegments object
segs = r.rs.segments

nodes = r.rs.nodes
for i in range(0, len(nodes)):
    nodes[i] = vector_3d(np.array(nodes[i]) / 100.)
r.rs.nodes = nodes

""" Mixed plot """
ana = pb.SegmentAnalyser(r.rs)
pd = vp.segs_to_polydata(ana, 1.e-2, ["radius", "subType", "creationTime"])
print("Root system bounds", pd.GetBounds())
rootActor = vp.plot_roots(pd, "creationTime", False)

ug = vp.read_vtu("benchmark3d_2-00001.vtu")
print("Mesh bounds", ug.GetBounds())
p_name = "water content"  # "pressure head"  # e.g. "S_liq" "water content"
meshActor = vp.plot_mesh_wireframe(ug, "", False)

vp.render_window([rootActor, meshActor], "mixed fun")
