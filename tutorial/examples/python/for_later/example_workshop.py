import py_rootbox as rb
from rb_tools import *
import matplotlib.pyplot as plt

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "workshop"  # "Anagallis_femina_Leitner_2010"  Zea_mays_4_Leitner_2014
rs.openFile(name)

# Initialize
rs.initialize()

# Simulate
simtime = 8
rs.simulate(simtime, True)

nodes = rs.getNodes()
seg = rs.getSegments()
times = rs.getNETimes(False)

# Remove artificial shoot
nodes = nodes[1:]
times = times[1:]
for s in seg:
    s.x -= 1
    s.y -= 1

# Type
segO = rs.getSegmentsOrigin()
types = list(range(0, len(seg)))
for i, s in enumerate(segO):
    types[i] = s.param.type

print()
print("Number of nodes", len(nodes))
print("Number of segments", len(seg))
print("Number of ages", len(times))
print("Number of types", len(types))
print()

str_n = "nodes = [ "
for i in range(0, len(nodes)):
    n = nodes[i]
    str_n += "["
    str_n += '{:1.2f}'.format(n.x)
    str_n += ","
    str_n += '{:1.2f}'.format(n.y)
    str_n += ","
    str_n += '{:1.2f}'.format(n.z)
    str_n += "], "
str_n = str_n[0:-2] + " ]"

str_s = "seg = [ "
for i in range(0, len(seg)):
    s = seg[i]
    str_s += "["
    str_s += '{:d}'.format(s.x)
    str_s += ","
    str_s += '{:d}'.format(s.y)
    str_s += "], "
str_s = str_s[0:-2] + " ]"

str_t = "age = [ "
for i in range(0, len(times)):
    t = times[i]
    str_t += '{:1.2f}'.format(simtime - t)
    str_t += ","
str_t = str_t[0:-1] + " ]"

str_tt = "types = [ "
for i in range(0, len(types)):
    t = types[i]
    str_tt += '{:d}'.format(t)
    str_tt += ","
str_tt = str_tt[0:-1] + " ]"

print(str_n)
print(str_s)
print(str_t)
print(str_tt)

cols = ['r', 'b', 'm', 'm']
for i, s in enumerate(seg):
    n1 = nodes[s.x]
    n2 = nodes[s.y]
    t = types[i] - 1
    plt.plot([n1.x, n2.x], [n1.z, n2.z], cols[t])

plt.axis('equal')
plt.show()

