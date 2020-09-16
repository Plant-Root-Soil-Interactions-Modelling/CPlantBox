"""small example"""
import sys; sys.path.append("../..")
import plantbox as pb
import vtk_plot as vp

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = ""
name = "kutschera"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeLB(1, 5)  # change basal to tap

# Simulate
rs.simulate(75, True)

# Export final result (as vtp)
rs.write("kutschera.vtp")

# Plot, using vtk
vp.plot_roots(rs, "creationTime")  # "creationTime"

