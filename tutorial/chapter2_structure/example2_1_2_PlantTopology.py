# <2_1_AStart>
import sys; sys.path.append("../.."); sys.path.append("../../src/") # |\label{l2_1:importStart}|
path = "../../modelparameter/structural/plant/"
sys.path.append(path)
import plantbox as pb
import visualisation.vtk_plot as vp
from example1f import template_text # |\label{l2_1:importEnd}|


### Simple example
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1.xml", verbose = True)

rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot #|\label{l2_1:arrayStart1}|
rrp.successorOT = [[2]] 
rrp.successorST = [[2]] 
rrp.successorP  = [[1.0]] 
rrp.successorNo = [1] 
rrp.successorWhere  = [[]] #|\label{l2_1:emptywhere}|

srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
srp.successorOT = [[3]] 
srp.successorST = [[2]] 
srp.successorP  = [[1.0]] 
srp.successorNo = [1] 
rrp.successorWhere  = [[]] #|\label{l2_1:arrayEnd1}|

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1a.vtp")
# <2_1_AEnd>

# ## Several successor types, specific locations    
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1.xml")

rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot #|\label{l2_1:arrayStart2}|
rrp.successorOT     = [[2], [4], [3]] 
rrp.successorST     = [[2], [1], [2]] 
rrp.successorNo     = [4, 1, 6] 
rrp.successorP      = [[1.0], [1.0], [1.0]] 
rrp.successorWhere  = [[-1.0, -3.0, -5.0, -7.0], [], []]

srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
srp.successorOT     = [[2], [4]] 
srp.successorST     = [[1], [1]] 
srp.successorNo     = [1, 4] 
srp.successorP      = [[1.0], [1.0]] 
srp.successorWhere  = [[4.0, 6.0, 8.0], [-3.0]] #|\label{l2_1:arrayEnd2}|


p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1b.vtp")

### Several successor types, specific locations and probabilistic branching    
p = pb.MappedPlant(2)
p.readParameters(path + "example2_1.xml")

rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot #|\label{l2_1:arrayStart3}|
rrp.successorOT     = [[2]] 
rrp.successorST     = [[2]] 
rrp.successorP      = [[1.0]] 
rrp.successorNo     = [1] 
rrp.successorWhere  = [[3.0, 5.0]]

srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
srp.successorOT     = [[3, 4, 2], [3]] 
srp.successorST     = [[2, 1, 3], [2]] 
srp.successorP      = [[0.2, 0.3, 0.3], [1.0]] 
srp.successorNo     = [4, 1] 
srp.successorWhere  = [[-3.0], [3.0]] #|\label{l2_1:arrayEnd3}|

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example2_1c.vtp")
