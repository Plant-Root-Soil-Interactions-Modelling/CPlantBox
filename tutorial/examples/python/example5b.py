import py_rootbox as rb
from rb_tools import *

# Simulate a root system
name = "Zea_mays_4_Leitner_2014"  # "Zea_mays_4_Leitner_2014"  # "Anagallis_femina_Leitner_2010" #"Sorghum_bicolor_NA_NA" "Zea_mays_4_Leitner_2014" #
rs = rb.RootSystem()
rs.openFile(name)
rs.initialize()

# example
simtime = 60  # days
dt = 0.1
N = round(simtime / dt)

# Incrementally build nodes and segments
nodes = vv2a(rs.getNodes())  # contains the initial nodes of tap, basal and shootborne roots
seg = np.array([], dtype = np.int64).reshape(0, 2)

print(nodes)

print()
for i in range(0, N):

     rs.simulate(dt, False)
     print("Number of nodes", rs.getNumberOfNodes())  # equals the number of new segments
     print("Number of roots", len(rs.getRootTips()))

     print("Number of new nodes", rs.getNumberOfNewNodes())  # equals the number of new segments
     print("Number of new roots", rs.getNumberOfNewOrgans())

     uni = v2ai(rs.getUpdatedNodeIndices())
     unodes = vv2a(rs.getUpdatedNodes())
     print("Number of node updates", len(unodes), len(uni))

     nodes[uni] = unodes  # do the update

     newnodes2 = rs.getNewNodes()
     newnodes = vv2a(newnodes2)
     newsegs = seg2a(rs.getNewSegments())
     if len(newnodes) != 0:
         print(nodes)
         print(newnodes)
         nodes = np.vstack((nodes, newnodes))
     if len(newsegs) != 0:
         seg = np.vstack((seg, newsegs))

     print()

# test is everything right?
nodes_ = vv2a(rs.getNodes());
seg_ = seg2a(rs.getSegments());

uneq = np.sum(nodes_ != nodes) / 3
print("unequal nodes: ", uneq)
if uneq > 0:
    i = np.nonzero(nodes_[:, 0] != nodes[:, 0])
    print(i)
    print()
    print(nodes[i, :])
    print(nodes_[i, :])
    print()

# node indices have meaning, the ordering should be the same

# segment indices have no special meaning, and the ordering is different
seg = np.sort(seg, axis = 0)  # per default along the last axis
seg_ = np.sort(seg_, axis = 0)
print("unequal segs: ", np.sum(seg_ != seg) / 2)

# rs.write("results/example_5b.vtp")

#      if rs.getNumberOfNewNodes()!=newnodes.shape[0]: # stuff for debugging...
#          print("oh noooooooo")
#      for i in range(0,newnodes.shape[0]):
#          if np.sum(newnodes[i,:]) == 0:
#              i_=i+nodes.shape[0]
#              print("New node ", i, "is empty =",i_)
#
#              nodes_ = vv2a(rs.getNodes())
#              seg_ = seg2a(rs.getSegments())
#              origins = rs.getSegmentsOrigin()
#              # print(type(origins))
#
#              print("should be ", nodes_[i_])
#
#              arg = np.argmax(seg_[:,1] == i_)
#              print("missing segment ", arg, seg_[arg,:])
#              r = origins[int(arg)]
#              print(r)
#              print(r.param)
#              arg = np.argmax(seg_[:,1] == seg_[arg,0])
#              print("predecessor ", arg, seg_[arg,:])
#
#              arg = np.argmax(newsegs[:,1] == i_)
#              print(arg)
#              print("segment ", newsegs[arg,:])
#
#              print("New segments")
#              print(newsegs)

