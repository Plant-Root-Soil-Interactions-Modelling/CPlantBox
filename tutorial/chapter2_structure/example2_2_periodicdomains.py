"""increase axial resolution (e.g. for animation)"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
filename = "Anagallis_femina_Leitner_2010"
plant.readParameters(path + filename + ".xml")

# Simulate
plant.initialize()
plant.simulate(60)  # days

# Export results as segments
ana = pb.SegmentAnalyser(plant)
ana.mapPeriodic(15, 10)  # cm |\label{l2_2_3:periodicity}|
ana.write("results/periodic.vtp")

# Plot final (periodic) image, using vtk
vp.plot_roots(ana, "creationTime")
