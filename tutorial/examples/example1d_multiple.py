"""multiple root systems"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

path = path = "../../modelparameter/structural/rootsystem/"
name = "Zea_mays_4_Leitner_2014"

simtime = 120
N = 3  # number of columns and rows
dist = 40  # distance between the root systems [cm]

# Initializes N*N root systems
allRS = []
for i in range(0, N):
    for j in range(0, N):
        rs = pb.RootSystem()
        rs.readParameters(path + name + ".xml")
        rs.getRootSystemParameter().seedPos = pb.Vector3d(dist * i, dist * j, -3.)  # cm
        rs.initialize(False)  # verbose = False
        allRS.append(rs)

# Simulate
for rs in allRS:
    rs.simulate(simtime, False)  # verbose = False

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()  # see example 3b
for i, rs in enumerate(allRS):
      vtpname = "results/example_1d_" + str(i) + ".vtp"
      rs.write(vtpname)
      ana.addSegments(rs)  # collect all

# Write all into single file (segments)
ana.write("results/example_1d_all.vtp")

# Plot, using vtk
vp.plot_roots(ana, "radius")
