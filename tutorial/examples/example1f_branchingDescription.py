import sys; sys.path.append("../.."); sys.path.append("../../src/")
path = "../../modelparameter/structural/plant/"
sys.path.append(path)
import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from example1f import template_text


'''
Definition of organ successors.
Arguments:
    ruleId: id [int] of the rule, in case the parameters are divided between several lines
    organType: one or several organ type(s) of the successor [int or array]. 
    subType: one or several subtype(s) of the successor [int or array]
    numLat: number of lateral to create at each branching point [int]
    probability: one value for each successor given. Use for probabilistic branching
    where: at which linking node is the rule applyed. The number and position of the linking nodes are defined by lmax, la, lb and ln.

Below, we set several parameter file describing how different branching patterns can be defined
The plant grows one tap root and one main shoot. The dictionnary below deine the successor of those two organs
'''

changeParamFile = False # re-define the topology via the parameter file [True] or change directly the random parameters [False]

# ## Simple example
'''
example with one type of successor
'''
if changeParamFile:
    successors = {
            "root":"""<parameter name="successor" ruleId="0" organType="2" subType="2" probability="1"/>""",
            "stem":"""<parameter name="successor" ruleId="0" organType="3" subType="2" probability="1"/>"""
            }

    filledText = template_text.format(**successors)
    f = open(path + "example1f.xml", "w")
    f.write(filledText)
    f.close()
    
p = pb.MappedPlant(2)
p.readParameters(path + "example1f.xml", verbose = True)

if not changeParamFile:
    rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot
    rrp.successorOT = [[2]] 
    rrp.successorST = [[2]] 
    rrp.successorP  = [[1.0]] 
    rrp.successorNo = [1] 

    srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
    srp.successorOT = [[3]] 
    srp.successorST = [[2]] 
    srp.successorP  = [[1.0]] 
    srp.successorNo = [1] 

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example_1fa.vtp")

# ## Several successor types, specific locations
'''
giving positive integers for "where" define the linking node to include (here "root" successor). Giving negative integers to "where" define linking nodes to exclude  (here 1st "stem" successor).
By default, the rule is applyed at all linking nodes   (here 2nd "stem" successor).
via "numlat" several laterals can grow at each selected linking nodes.
Roots can grow out of shoot organs and vice-versa
 
'''
if changeParamFile:
    successors = {
            "root":
    """<parameter name="successor" ruleId="0" numLat="4" where="-1,-3,-5,-7" organType="2" subtype="2" probability="1"/>
        <parameter name="successor" ruleId="1" numLat="1" organType="4" subtype="1" probability="1"/>
        <parameter name="successor" ruleId="2" numLat="6" organType="3" subtype="2" probability="1"/>""",
            "stem":"""<parameter name="successor" ruleId="0" where="4,6,8" organType="2" subtype="1" probability="1"/>
            <parameter name="successor" ruleId="1" numLat="4" where="-3" organType="4" subtype="1" probability="1"/>"""
            }

    filledText = template_text.format(**successors)
    f = open(path + "example1f.xml", "w")
    f.write(filledText)
    f.close()
p = pb.Plant(2)
p.readParameters(path + "example1f.xml")

if not changeParamFile:
    rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot
    rrp.successorOT     = [[2], [4], [3]] 
    rrp.successorNo     = [4, 1, 6] 
    rrp.successorST     = [[2], [1], [2]] 
    rrp.successorP      = [[1.0], [1.0], [1.0]] 
    rrp.successorWhere  = [[-1.0, -3.0, -5.0, -7.0], [], []]

    srp = p.getOrganRandomParameter(pb.stem)[1] # define laterals of stem
    srp.successorOT     = [[2], [4]] 
    srp.successorST     = [[1], [1]] 
    srp.successorP      = [[1.0], [1.0]] 
    srp.successorNo     = [1, 4] 
    srp.successorWhere  = [[4.0, 6.0, 8.0], [-3.0]]

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")

# Export final result (as vtp)
p.write("results/example_1fb.vtp")

# ## Several successor types, specific locations and probabilistic branching
'''
We can have probalistic branching for the selection of the successor via the "probability" parameter (here 1st "stem" successor).
The probability for one successor rule must be < or = 1.
'''
if changeParamFile:
    successors = {
            "root":"""<parameter name="successor" ruleId="0" where="3, 5" organType="2" subtype="2" probability="1"/>""",
            "stem":"""<parameter name="successor" ruleId="0" where="-3" numLat="4" organType="3,4,2" subtype="2,1,3" probability="0.2,0.3,0.3"/>
            <parameter name="successor" ruleId="1" where="3" numLat="1" organType="3" subtype="2" probability="1"/>"""
            }

    filledText = template_text.format(**successors)
    f = open(path + "example1f.xml", "w")
    f.write(filledText)
    f.close()
    
p = pb.MappedPlant(2)
p.readParameters(path + "example1f.xml")

if not changeParamFile:
    rrp = p.getOrganRandomParameter(pb.root)[1] # define laterals of taproot
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
    srp.successorWhere  = [[-3.0], [3.0]]

p.initialize(False)
time = 100
p.simulate(time, False)
vp.plot_plant(p, "organType")
p.write("results/example_1fc.vtp")
