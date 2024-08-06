""" plots simulation of parameter set m1_lupine.xml and maps it to 2D """
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp


# sets all standard deviation to a percantage, i.e. value*s
def set_all_sd(rs, s):
    for p in rs.getRootRandomParameter():
        p.lmaxs = p.lmaxs * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.rs = p.r * s
        p.a_s = p.a * s


rs = pb.RootSystem()

# Open plant and root parameter from a file
path = ""
name = "m1_lupine"
rs.readParameters(path + name + ".xml")
# set_all_sd(rs, 0.)

# print(rs.getRootRandomParameter(0))
print(rs.getRootRandomParameter(1))
print(rs.getRootRandomParameter(2))

# Initialize
rs.initializeLB(1, 5)  # change basal to tap

# Simulate
rs.simulate(11, True)
ana = pb.SegmentAnalyser(rs)
ana.map2D()

# Export final result (as vtp)
ana.write("m1_lupine.vtp")

# Plot, using vtk
vp.plot_roots(ana, "creationTime")  # "creationTime", "subType"
