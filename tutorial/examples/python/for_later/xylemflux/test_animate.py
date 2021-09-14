""" 
Tests if we can animate root growth (only roots, nos soil)
"""
import sys
sys.path.append("zzz../../../..")

import plantbox as pb
import rsml_reader as rsml
import vtk_plot as vp

from math import *
import numpy as np
import matplotlib.pyplot as plt

""" root problem """
rs = pb.RootSystem()
path = "../../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")
for p in rs.getRootRandomParameter():  # Modify axial resolution
    p.dx = 5  # [cm] adjust resolution
rs.initialize()

""" initial growth """
rs.simulate(1)
ana = pb.SegmentAnalyser(rs)
pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime"])
rootActor, rootCBar = vp.plot_roots(pd, "creationTime", False)
iren = vp.render_window(rootActor, "Animation", rootCBar)

c = 0
max_age = 60


def timer_callback_(obj, ev):
    """ animation call back function (called every 0.1 second) """
    global rootActor
    global c

    c += 1
    print("hello", c)
    rs.simulate(1)
    ana = pb.SegmentAnalyser(rs)
    pd = vp.segs_to_polydata(ana, 1., ["radius", "subType", "creationTime"])
    newRootActor, rootCBar = vp.plot_roots(pd, "creationTime", False)
    renWin = iren.GetRenderWindow()
    ren = renWin.GetRenderers().GetFirstRenderer ()
    ren.RemoveActor(rootActor)
    newRootActor.RotateX(-90)
    ren.AddActor(newRootActor)
    ren.ResetCamera()
    rootActor = newRootActor
    iren.Render()

    if c >= max_age:
        c = 0
        rs.initialize()


""" initial """
iren.AddObserver('TimerEvent', timer_callback_, 1.0)
iren.Start()  # Start the event loop.

print("the end")
