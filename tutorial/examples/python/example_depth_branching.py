import py_rootbox as rb
import math


#
#
#
class My_Soil(rb.SoilLookUp):

    def getValue(self, pos, root):
        print(pos.z)
        return 0.1 + abs(pos.z / 100.)

    def __str__(self):
        return "mysoil"


rs = rb.RootSystem()
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)

soil = My_Soil()

# Manually set scaling function and tropism parameters
sigma = [0.4, 1., 1., 1., 1. ] * 2
for i in range(0, 10):
    p = rs.getRootTypeParameter(i + 1)
    p.dx = 0.25  # adjust resolution
    p.tropismS = sigma[i]
    p.ln = p.ln / 10
    p.nob = p.nob * 10
    p.lns = 0  # variation comes with branching propability anyway
    p.sbp = soil

# Simulation
rs.initialize()
simtime = 120.
dt = 1.
N = simtime / dt
for i in range(0, round(N)):
    rs.simulate(dt, True)

rs.write("results/example_depth_branching.vtp")  # root system

