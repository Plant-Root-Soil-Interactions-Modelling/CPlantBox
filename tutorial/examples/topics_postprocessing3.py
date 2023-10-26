"""dgf and vtp, and rsml export example TODO change to plant"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

plant = pb.Plant()
path = "../../modelparameter/structural/rootsystem/"
name = "wheat"  # "Brassica_napus_a_Leitner_2010"  # "Brassica_napus_a_Leitner_2010"  # "Anagallis_femina_Leitner_2010"  #
plant.readParameters(path + name + ".xml")

rhizotron = pb.SDF_PlantBox(7, 7, 14)
plant.setGeometry(rhizotron)  # soilcore, or rhizotron
plant.initialize()
plant.simulate(30, True)

ana = pb.SegmentAnalyser(plant)

# aseg = plant.getShootSegments()  # if there are no shoot borne roots, it is only one segment
# for s in aseg:
#     print("Shoot segment", s)
#     ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

ana.write("results/example_3c.vtp", ["radius", "surface"])
ana.write("results/example_3c.dgf")

# segment analyser cannot write rsml files becasue rsml is based on polylines, not segments
# use RootSystem::write to export a RSML
plant.write("results/example_3c.rsml")
