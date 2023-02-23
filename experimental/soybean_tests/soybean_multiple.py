"""multiple soybean root systems"""
import sys
sys.path.append("../.."); sys.path.append("../../src/python_modules")
import plantbox as pb
import vtk_plot as vp

path = "../../modelparameter/rootsystem/"
name = "Glycine_max_Moraes2020_opt2"

simtime = 154
N = 17  # number of columns
M = 3  # number of rows
dist = 38  # inter-row distance [cm]
distp = 6  # inter-plant distance within the rows [cm]

# Initializes N*M root systems
allRS = []
for i in range(0, N):
    for j in range(0, M):
        rs = pb.RootSystem()
        rs.readParameters(path + name + ".xml")
        rs.getRootSystemParameter().seedPos = pb.Vector3d(distp * i, dist * j, -3.)  # cm
        rs.initialize(False)  # verbose = False
        allRS.append(rs)

# Simulate
rs.setSeed(2)
for rs in allRS:
    rs.simulate(simtime, False)  # verbose = False

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()  # see example 3b
for i, rs in enumerate(allRS):
      vtpname = "results/" + name + "/" + str(i) + ".vtp"
      rs.write(vtpname)
      ana.addSegments(rs)  # collect all

# Write all into single file (segments)
ana.write("results/" + name + "/" + "soybean_all.vtp")

# Plot, using vtk
vp.plot_roots(ana, "radius", True, 'oblique')
