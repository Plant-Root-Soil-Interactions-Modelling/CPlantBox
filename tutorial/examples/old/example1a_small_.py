"""small example"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp

rs = pb.Plant()

# Open plant and root parameter from a file
path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# stem = pb.StemRandomParameter(rs)
# stem.subType = 1
# rs.setOrganRandomParameter(stem)

# Initialize
rs.initialize()

# Simulate
rs.simulate(30, True)

# for o in rs.getOrgans():
#     print(o.param())

# Export final result (as vtp)
rs.write("results/example_1a.vtp")

# Plot, using vtk
vp.plot_roots(rs, "subType")
