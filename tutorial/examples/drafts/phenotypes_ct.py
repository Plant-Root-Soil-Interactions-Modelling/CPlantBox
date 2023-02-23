#
# mimic root system part simulation from rb_swbot, rb_swtop
#

import py_rootbox as rb
from rb_tools import *
import matplotlib.pyplot as plt
import math


def simulate(name, simtime):
    return rs


names = ["maize_p1_zero_std", "maize_p2_zero_std", "maize_p3_zero_std"]

rs = rb.RootSystem()
rs.openFile(names[0], "params/")
rs.initialize()
rs.simulate(2, True)  # initial growth

nodes = vv2a(rs.getNodes())  # contains the initial nodes of tap, basal and shootborne roots
seg = np.array([], dtype = np.int64).reshape(0, 2)
cts = v2a(rs.getSegmentCTs())
nonm = 0
simtime = 60  # days
dt = 0.1
N = round(simtime / dt)

for i in range(0, N):
    rs.simulate(dt, False)
    uni = v2ai(rs.getUpdatedNodeIndices())  # MOVED NODES
    unodes = vv2a(rs.getUpdatedNodes())
    ucts = vv2a(rs.getUpdatedNodes())
    nodes[uni] = unodes  # do the update
    nonm += uni.shape[0]
    newnodes2 = rs.getNewNodes()  # NEW NODES
    newnodes = vv2a(newnodes2)
    newsegs = seg2a(rs.getNewSegments())  # NEW SEGS
    newcts = v2a(rs.getNewSegmentCTs())

    if len(newnodes) != 0:
        nodes = np.vstack((nodes, newnodes))
    if len(newsegs) != 0:
        seg = np.vstack((seg, newsegs))
        cts = np.vstack((cts, newcts))
        simtime = rs.getSimTime()
        print("Simtime", simtime, " max cts", np.max(newcts))

nodes_ = vv2a(rs.getNodes())
nodeCTs_ = v2a(rs.getNodeCTs())
seg_ = seg2a(rs.getSegments())
sct_ = v2a(rs.getSegmentCTs())

rs.write("rootsystem.vtp")

ana = rb.SegmentAnalyser(rs)
ana.write("rootsystem_segs.vtp")

print("---")
print("Creation times range from ", np.min(cts), " to ", np.max(cts), "len", cts.shape)
print("Creation times range from ", np.min(sct_), " to ", np.max(sct_), "len", sct_.shape)
print("seg:", seg_[0])
print(nodes_[seg[0][0]], nodes_[seg[0][1]])
print(nodeCTs_[seg[0][0]], nodeCTs_[seg[0][1]])
print("---")

# assertEqual(nodes_.shape, nodes.shape, "incremental growth: node lists are not equal")
# ind = np.argwhere(nodes_[:, 1] != nodes[:, 1])
# for i in ind:
#     print(i, nodes_[i], "!=", nodes[i])
#
# uneq = np.sum(nodes_ != nodes) / 3
# self.assertEqual(uneq, 0, "incremental growth: node lists are not equal")
# self.assertEqual(seg_.shape, seg.shape, "incremental growth: segment lists are not equal")
# seg = np.sort(seg, axis = 0)  # per default along the last axis
# seg_ = np.sort(seg_, axis = 0)
# uneq = np.sum(seg_ != seg) / 2
# self.assertEqual(uneq, 0, "incremental growth: segment lists are not equal")
