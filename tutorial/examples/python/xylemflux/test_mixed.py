""" 
Tests if we can combine different plots in vtk
"""
import sys
sys.path.append("../../../..")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np


def vector_3d(a):
    return pb.Vector3d(a[0], a[1], a[2])


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
rootActor, rootCBar = vp.plot_roots(pd, "creationTime", False)

ug = vp.read_vtu("benchmark3d_2-00001.vtu")
print("Mesh bounds", ug.GetBounds())
meshActor, meshCBar = vp.plot_mesh(ug, "water content", "", False)  # "pressure head"  # e.g. "S_liq" "water content"

vp.render_window([rootActor, meshActor], "mixed fun", meshCBar).Start()


# # Plot, using vtk
# rootActor, cBar = vp.plot_roots(rs, "creationTime", False)
# rootActor.RotateX(90) # to look at it from top
# vp.render_window(rootActor,"top view", cBar).Start()

