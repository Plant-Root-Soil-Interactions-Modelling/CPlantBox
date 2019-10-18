""" analysis of nodes and segments from measurements """
import sys
sys.path.append("../../..")
import plantbox as pb

nodes = [ [0, 1, 0], [0, 1, -1], [0, 1, -2], [0, 1, -3], ]
segs = [ [0, 1], [1, 2], [2, 3] ]
cts = [0., 0., 0.]
radii = [ 0.1, 0.1, 0.1]

# convert from python to c++ types
nodes = [pb.Vector3d(n[0], n[1], n[2]) for n in nodes]
segs = [pb.Vector2i(s[0], s[1]) for s in segs]

# use C++ without RootSystem
ana = pb.SegmentAnalyser(nodes, segs, cts, radii)

print("length", ana.getSummed("length"))

ana.write("example_3e.vtp", ["radius"])  # working
# ana.write("test.vtp") # not working, but with a meaningful exception
