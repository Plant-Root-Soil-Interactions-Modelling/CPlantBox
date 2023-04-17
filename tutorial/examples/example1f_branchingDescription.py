
import sys; sys.path.append("../.."); sys.path.append("../../src/")
path =  "../../modelparameter/structural/plant/"
sys.path.append( path)
import plantbox as pb
import visualisation.vtk_plot as vp
from example1f import template_text

'''
Definition of organ successors.
Arguments:
    ruleId: id [int] of the rule, in case the parameters are divided between several lines
    organType: one or several organ type(s) of the successor [int or array]. ATT: roots can grow out of shoot organs but not vice-versa
    subType: one or several subtype(s) of the successor [int or array]
    numLat: number of lateral to create at each branching point [int]
    probability: one value for each successor given. Use for probabilistic branching
    where: at which linking node is the rule applyed. The number and position of th elinking nodes are defined by lmax, la, lb and ln.

Below, we set several parameter file describing how different branching patterns can be defined
The plant grows one tap root and one main shoot. The dictionnary below deine the successor of those two organs
'''

### Simple example
'''
example with one type of successor
'''
successors= {
        "root":"""<parameter name="successor" ruleId="0" organType="2" subType="2" probability="1"/>""",
        "stem":"""<parameter name="successor" ruleId="0" organType="3" subType="2" probability="1"/>"""
        }
        
filledText = template_text.format(**successors)
f = open(path+"example1f.rsml", "w")
f.write(filledText)
f.close()
p = pb.MappedPlant(2)
p.readParameters(path+"example1f.rsml") #the parameters are given directly as string and not via a text file        
p.initialize(False)
time = 100
p.simulate(time, False)
ana = pb.SegmentAnalyser()  # see example 3b

vp.plot_plant(p, "organType")

### Several successor types, specific locations and probabilistic branching
'''
giving positive integers for "where" define the linking node to include (here "root" successor). Giving negative integers to "where" define linking nodes to exclude  (here 1st "stem" successor).
By default, the rule is applyed at all linking nodes   (here 2nd "stem" successor).

We can have probalistic branching for the selection of the successor via the "probability" parameter (here 1st "stem" successor).
The probability for one successor rule must be < or = 1. 
'''
successors= {
        "root":"""<parameter name="successor" ruleId="0" where="3, 5" organType="2" type="2" probability="1"/>""",
        "stem":"""<parameter name="successor" ruleId="0" where="-3" numLat="4" organType="3,4,2" type="2,1,3" probability="0.2,0.3,0.3"/>
        <parameter name="successor" ruleId="1" numLat="1" organType="3" type="2" probability="1"/>"""
        }
        
filledText = template_text.format(**successors)
f = open(path+"example1f.rsml", "w")
f.write(filledText)
f.close()
p = pb.Plant(2)
p.readParameters(path+"example1f.rsml") #the parameters are given directly as string and not via a text file        
p.initialize(False)
time = 100
p.simulate(time, True)
ana = pb.SegmentAnalyser()  # see example 3b

vp.plot_plant(p, "organType")

### Several successor types, specific locations and probabilistic branching
'''
via "numlat" several laterals can grow at each selected linking nodes.
Roots can grow shoot organs and vice-versa
'''

successors= {
        "root":
"""<parameter name="successor" ruleId="0" numLat="4" where="-1,-3,-5,-7" organType="2" type="2" probability="1"/>
    <parameter name="successor" ruleId="1" numLat="1" organType="4" type="1" probability="1"/>
    <parameter name="successor" ruleId="2" numLat="6" organType="3" type="2" probability="1"/>""",
        "stem":"""<parameter name="successor" ruleId="0" where="4,6,8" organType="2" type="1" probability="1"/>
        <parameter name="successor" ruleId="1" numLat="4" where="-3" organType="4" type="1" probability="1"/>"""
        }
        
filledText = template_text.format(**successors)
f = open(path+"example1f.rsml", "w")
f.write(filledText)
f.close()
p = pb.MappedPlant(2)
p.readParameters(path+"example1f.rsml") #the parameters are given directly as string and not via a text file        
p.initialize(False)
time = 100
p.simulate(time, False)
ana = pb.SegmentAnalyser()  # see example 3b

vp.plot_plant(p, "organType")