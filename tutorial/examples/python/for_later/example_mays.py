import py_rootbox as rb
from rb_tools import *

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Zea_mays_4_Leitner_2014"  # "Anagallis_femina_Leitner_2010"

rs.openFile(name)

# maxB = 5
# firstB = 7.
# delayB = 7.
# rsp = rb.RootSystemParameter()
# rsp.set(-3., firstB, delayB, maxB, 0, 1.e9, 1.e9,  1.e9, 0., 0.)
# rs.setRootSystemParameter(rsp)

# Initialize
rs.initialize()

# Simulate
rs.simulate(60, True)

ana = rb.SegmentAnalyser(rs)
l = ana.getSummed("length")
print("Total root length", l, "cm")

# Export final result (as vtp)
rs.write("results/example_P1.vtp")
