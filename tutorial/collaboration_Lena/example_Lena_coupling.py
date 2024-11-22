""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.root_conductivities import *  # hard coded conductivities
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import math
import numpy as np
import matplotlib.pyplot as plt
import timeit
import csv

def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" Parameters """
min_b = [-4., -4., -25.] 
max_b = [4., 4., 0.]
cell_number = [16, 16, 30] #[4,4,12] #[8, 8, 25]  # [16, 16, 30]  # [32, 32, 60]
cellvol = ((max_b[0]-min_b[0])/cell_number[0])*((max_b[1]-min_b[1])/cell_number[1])*((max_b[1]-min_b[1])/cell_number[1])
periodic = False

path = "../../modelparameter/structural/rootsystem/"
name = "Faba_synMRI" #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -1000  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm
wcroot = 0.85 #water content within the root 

sim_time = 1  # [day] for task b
rs_age = 10  # root system initial age
age_dependent = False  # conductivities
dt = 120. / (24 * 3600)  # [days] Time step must be very small

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
s.setHomogeneousIC(initial, True)  # cm pressure head, equilibrium
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.initializeProblem()
s.setCriticalPressure(wilting_point)

""" Initialize xylem model """
rs = pb.MappedRootSystem()
rs.readParameters(path + name + ".xml")
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
else:
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, max_b[2] - min_b[2])
rs.setGeometry(sdf)
rs.setSeed(1)
r = XylemFluxPython(rs)
init_conductivities(r, age_dependent)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
r.rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)
rs.initialize()
rs.simulate(rs_age, False)
r.test()  # sanity checks
nodes = r.get_nodes()
cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution """
start_time = timeit.default_timer()
x_, y_ = [], []
sx = s.getSolutionHead()  # inital condition, solverbase.py
N = round(sim_time / dt)
t = 0.

for i in range(0, N):

    rs.simulate(dt)

    rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point)  # xylem_flux.py
    x_.append(t)
    y_.append(float(r.collar_flux(rs_age + t, rx, sx)))  # exact root collar flux

    fluxes = r.soilFluxes(rs_age + t, rx, sx, False)
    s.setSource(fluxes)  # richards.py
    s.solve(dt)
    sx = s.getSolutionHead()  # richards.py

    min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
    n = round(float(i) / float(N) * 100.)
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g} days {:g}"
            .format(min_sx, max_sx, min_rx, max_rx, s.simTime, rx[0]))
    t += dt

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" VTK visualisation """
#vp.plot_roots_and_soil(r.rs, "water content", rx, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)


#soil
wc = np.array(s.getWaterContent())
wc = np.reshape(wc, (cell_number[2],cell_number[0],cell_number[1]))
#wc = np.swapaxes(wc, 0, 2)
wc = wc[::-1]

#root
ana = pb.SegmentAnalyser(rs.mappedSegments())
segs = ana.segments
nodes = ana.nodes
radius = ana.getParameter("radius")
seglen = ana.getParameter("length") 

C = np.zeros((cell_number[0]*cell_number[1]*cell_number[2]))
x = np.zeros(len(segs))
for i, s in enumerate(segs):
    try:
        x[i] = rs.seg2cell[i]
    except:  # in case the segment is not within the domain
        x[i] = -1

sc = np.unique(x).astype(int) #all voxels with a segment 
for i in range(0,len(sc)):
    sg = rs.cell2seg[sc[i]] #numbers of all segments within that voxel

    rootvol = 0
    for j in range(0,len(sg)):
        rootvol = rootvol + radius[sg[j]] ** 2 * math.pi * seglen[sg[j]]
    frac = rootvol / cellvol
    C[sc[i]] = frac * wcroot

C[C>1] = 1
C = np.reshape(C, (cell_number[2],cell_number[0],cell_number[1]))
C = C[::-1]

#save everything for plotting
np.savez('sims_root_soil_fine', wc, C, cell_number)

with open('tmp_file.txt', 'w') as f:
    csv.writer(f, delimiter=' ').writerows(data)

print('simulation finished')
