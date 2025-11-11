import sys; sys.path.append("../..")
sys.path.append("../../src/") # |\label{l2_1:importStart}|
path = "../../modelparameter/structural/plant/"
sys.path.append(path)
import plantbox as pb
import visualisation.vtk_plot as vp # |\label{l2_1:importEnd}|


### Simple example
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1_2.xml", verbose = True)

rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot #|\label{l2_1:arrayStart1}|
rrp.successorOT = [[pb.root]] 
rrp.successorST = [[2]] 
rrp.successorP  = [[1.0]] 
rrp.successorNo = [1] 
rrp.successorWhere  = [[]] #|\label{l2_1:emptywhere}|

srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
srp.successorOT = [[pb.stem]] 
srp.successorST = [[2]] 
srp.successorP  = [[1.0]] 
srp.successorNo = [1] 
rrp.successorWhere  = [[]] #|\label{l2_1:arrayEnd1}|

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1_2a.vtp")

### Several successor types, specific locations    
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1_2.xml")

rrp = p.getOrganRandomParameter(pb.root)[1] #|\label{l2_1:arrayStart2}|
rrp.successorOT     = [[pb.root], [pb.root]] 
rrp.successorST     = [[2], [3]] 
rrp.successorNo     = [1, 20] 
rrp.successorP      = [[1.0], [1.0]] 
rrp.successorWhere  = [[-1, -2, -3,-4, -5, -7], [7]]

srp = p.getOrganRandomParameter(pb.stem)[1]
srp.successorOT     = [[pb.root], [pb.leaf]] 
srp.successorST     = [[4], [1]] 
srp.successorNo     = [1, 4] 
srp.successorP      = [[1.0], [1.0]] 
srp.successorWhere  = [[0.], [-0.]] #|\label{l2_1:arrayEnd2}|


p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1_2b.vtp")

### Probabilistic branching    
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1_2.xml")

rrp = p.getOrganRandomParameter(pb.root)[1] #|\label{l2_1:arrayStart3}|
rrp.successorOT     = [[pb.root]] 
rrp.successorST     = [[2]] 
rrp.successorP      = [[0.7]] 
rrp.successorNo     = [1] 
rrp.successorWhere  = []

srp = p.getOrganRandomParameter(pb.stem)[1] 
srp.successorOT     = [[pb.stem, pb.leaf], [pb.stem]] 
srp.successorST     = [[2, 1], [2]] 
srp.successorP      = [[0.2, 0.8], [1.0]] 
srp.successorNo     = [4, 10] 
srp.successorWhere  = [[-3], [3]] #|\label{l2_1:arrayEnd3}|

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1_2c.vtp")
