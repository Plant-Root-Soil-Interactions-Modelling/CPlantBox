"""dgf and vtp, and rsml export example"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb

rs = pb.RootSystem()
path = "../../modelparameter/structural/rootsystem/"
name = "wheat"  # "Brassica_napus_a_Leitner_2010"  # "Brassica_napus_a_Leitner_2010"  # "Anagallis_femina_Leitner_2010"  #
rs.readParameters(path + name + ".xml")

rhizotron = pb.SDF_PlantBox(7, 7, 14)
rs.setGeometry(rhizotron)  # soilcore, or rhizotron

rs.initialize()
rs.simulate(30, True)

ana = pb.SegmentAnalyser(rs)

aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

ana.write("results/example_3c.vtp", ["radius", "surface"])
ana.write("results/example_3c.dgf")

# segment analyser cannot write rsml files becasue rsml is based on polylines, not segments
# use RootSystem::write to export a RSML
rs.write("results/example_3c.rsml")
