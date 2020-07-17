"""small example"""
import sys
sys.path.append("../../..")
import plantbox as pb
import vtk_plot as vp

rs = pb.RootSystem()

# Open plant and root parameter from a file
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")

# Initialize
rs.initializeDB(1, 4)   # numbers indicate basal and shoot borne root types 

# Simulate
rs.simulate(30, True)

# Export final result (as vtp)
rs.write("results/example_1a.vtp")

# Plot, using vtk
vp.plot_roots(rs, "creationTime")
