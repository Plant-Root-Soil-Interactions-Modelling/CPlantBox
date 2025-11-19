"""
small example for CI testing
original source: tutorial/examples/example1a_small.py

m.vianna@fz-juelich.de
Jun-2025
"""


import plantbox as pb
#import plantbox.visualisation.vtk_plot as vp

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

for p in rs.getOrganRandomParameter(pb.root):
    print(p.subType)
    print(p.name)

rs.writeParameters("test.xml")

# Simulate
rs.simulate(30, True)

# for o in rs.getOrgans():
#     print(o.param())

# Export final result (as vtp)
rs.write("results/example_1a.vtp")

# Plot, using vtk [switched-off as runners do not have x11]
# vp.plot_roots(rs, "subType")
