import sys
sys.path.append("..")
import plantbox as pb
from rb_tools import *

# set up simulation
rs = pb.RootSystem()
name = "Zea_mays_4_Leitner_2014"
rs.openFile(name)
soilcore = pb.SDF_PlantContainer(5,5,40,False)
rs.setGeometry(soilcore)  # soilcore, or rhizotron
rs.initialize()

# test push and pop
rs.simulate(20)
nodes0a = vv2a(rs.getNodes())
rs.push() # push simualtion at time 20

rs.simulate(100)
nodes1a = vv2a(rs.getNodes())

rs.pop() # should be at state 20 again
nodes0b = vv2a(rs.getNodes())
stime0 = rs.getSimTime()

rs.simulate(100)
nodes1b = vv2a(rs.getNodes())
stime1 = rs.getSimTime()

# check
print("\nResults ")

uneq0 = np.sum(nodes0a!=nodes0b)/3
print("Time", stime0,"has", uneq0, "unequal nodes")
print(nodes0a.shape, nodes0b.shape)

uneq1 = np.sum(nodes1a!=nodes1b)/3
print("Time", stime1,"has", uneq1, "unequal nodes")
print(nodes1a.shape, nodes1b.shape)
