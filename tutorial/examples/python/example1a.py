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
rootsystem.write("../results/example_1a.vtp")

print("done.")
