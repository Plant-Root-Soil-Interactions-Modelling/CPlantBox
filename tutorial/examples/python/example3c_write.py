"""dgf and vtp export example"""
import sys; sys.path.append("../../..")
import plantbox as pb

rs = pb.RootSystem()
path = "../../../modelparameter/rootsystem/"
name = "Anagallis_femina_Leitner_2010"
rs.readParameters(path + name + ".xml")
rs.initialize()
rs.simulate(15, True)

ana = pb.SegmentAnalyser(rs)

aseg = rs.getShootSegments()  # if there are no shoot borne roots, it is only one segment
for s in aseg:
    print("Shoot segment", s)
    ana.addSegment(s, 0., 0.1, True)  # ct, radius, insert first

ana.write("results/example_3c.vtp", ["radius", "surface"])
ana.write("results/example_3c.dgf")
