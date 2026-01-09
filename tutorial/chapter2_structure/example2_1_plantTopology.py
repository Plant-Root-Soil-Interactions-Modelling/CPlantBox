"""change parameters in the input script"""

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

path = "../../modelparameter/structural/plant/"

# Define a simple plant topology, part A
plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_1_2.xml", verbose = True)

rrp = plant.getOrganRandomParameter(pb.root)[1]  # laterals of taproot #|\label{l2_1:arrayStart1}|
rrp.successorOT = [[pb.root]]
rrp.successorST = [[2]]
rrp.successorP = [[1.0]]
rrp.successorNo = [1]
rrp.successorWhere = [[]]  # |\label{l2_1:emptywhere}|

srp = plant.getOrganRandomParameter(pb.stem)[1]  # laterals of stem
srp.successorOT = [[pb.stem]]
srp.successorST = [[2]]
srp.successorP = [[1.0]]
srp.successorNo = [1]
rrp.successorWhere = [[]]  # |\label{l2_1:arrayEnd1}|

plant.initialize(False)  # |\label{l2_1:simulateStart1}|
plant.simulate(100, False)
vp.plot_plant(plant, "organType")
plant.write("results/example2_1_2a.vtp")  # |\label{l2_1:simulateEnd1}|

# Several successor types, specific locations, part B
plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_1_2.xml")

rrp = plant.getOrganRandomParameter(pb.root)[1]  # |\label{l2_1:arrayStart2}|
rrp.successorOT = [[pb.root], [pb.root]]
rrp.successorST = [[2], [3]]
rrp.successorNo = [1, 20]
rrp.successorP = [[1.0], [1.0]]
rrp.successorWhere = [[-1, -2, -3, -4, -5, -7], [7]]  # |\label{l2_1:arrayEndroot2}|

srp = plant.getOrganRandomParameter(pb.stem)[1]
srp.successorOT = [[pb.root], [pb.leaf]]
srp.successorST = [[4], [1]]
srp.successorNo = [1, 4]
srp.successorP = [[1.0], [1.0]]
srp.successorWhere = [[0.0], [-0.0]]  # |\label{l2_1:arrayEndstem2}|

plant.initialize(False)
plant.simulate(100, False)
vp.plot_plant(plant, "organType")
plant.write("results/example2_1_2b.vtp")

# Probabilistic branching, part C
plant = pb.MappedPlant(2)
plant.readParameters(path + "example2_1_2.xml")

rrp = plant.getOrganRandomParameter(pb.root)[1]  # |\label{l2_1:arrayStart3}|
rrp.successorOT = [[pb.root]]
rrp.successorST = [[2]]
rrp.successorP = [[0.7]]  # |\label{l2_1:probabilityroot3}|
rrp.successorNo = [1]
rrp.successorWhere = []

srp = plant.getOrganRandomParameter(pb.stem)[1]
srp.successorOT = [[pb.stem, pb.leaf], [pb.stem]]
srp.successorST = [[2, 1], [2]]
srp.successorP = [[0.2, 0.8], [1.0]]  # |\label{l2_1:probabilitystem3}|
srp.successorNo = [4, 10]
srp.successorWhere = [[-3], [3]]  # |\label{l2_1:arrayEnd3}|

plant.initialize(False)
plant.simulate(100, False)

vp.plot_plant(plant, "organType")
plant.write("results/example2_1_2c.vtp")

