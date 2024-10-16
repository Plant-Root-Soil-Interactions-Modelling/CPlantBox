import sys; sys.path.append(".."); sys.path.append("../src/")
import plantbox as pb
import visualisation.vtk_plot as vp


plant_path = "structural/plant"
root_path = "structural/rootsystem"
filename = "Anagallis_femina_Leitner_2010"

plant1 = pb.Plant()
plant2 = pb.Plant()
simtime = 120


twoPlants = []
plant1.readParameters(root_path +"/"+ filename +".xml",verbose = True)
seed1 = plant1.getOrganRandomParameter(pb.seed)[0]
seed1.seedPos = pb.Vector3d(-20, 0, -3.)
plant1.initialize()
twoPlants.append(plant1)
plant2.readParameters(plant_path +"/"+ filename +".xml",verbose = True)
seed2 = plant2.getOrganRandomParameter(pb.seed)[0]
seed2.seedPos = pb.Vector3d(20, 0, -3.)
plant2.initialize()
twoPlants.append(plant2)
for rs in twoPlants:
    rs.simulate(simtime, False)  # verbose = False

# Export results as single vtp files (as polylines)
ana = pb.SegmentAnalyser()  # see example 3b
for i, rs in enumerate(twoPlants):
      ana.addSegments(rs)  # collect all


# Plot, using vtk
vp.plot_roots(ana, "radius")