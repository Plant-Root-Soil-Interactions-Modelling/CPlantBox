"""small example"""
import sys
sys.path.append("../../..")
import plantbox as pb

rootsystem = pb.RootSystem()

# Open plant and root parameter from a file
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rootsystem.readParameters(path + name + ".xml")

# Initialize
rootsystem.initialize()

# Simulate
rootsystem.simulate(30, True)

# Export final result (as vtp)
rootsystem.write("results/example_1a.vtp")

aseg = rootsystem.getShootSegments()
for s in aseg:
    print("Shoot", s)

# or as DGF
ana = pb.SegmentAnalyser(rootsystem)
for s in aseg:  # its actually only one
    ana.addSegment(s, 0, 0.1)  # ct, radius

ana.write("results/example_1b.dgf")
