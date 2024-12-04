""" coupling with DuMux as solver for the soil, dumux-rosi must be installed & compiled"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.xylem_flux import XylemFluxPython
from functional.root_conductivities import *  # hard coded conductivities

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" Parameters """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 25]  # [16, 16, 30]  # [32, 32, 60]
periodic = True

path = "../../modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -659.8 + 12.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm

sim_time = 7  # [day] for task b
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
    sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2])
rs.setGeometry(sdf)
rs.initialize()
rs.simulate(rs_age, False)
r = XylemFluxPython(rs)
init_conductivities(r, age_dependent)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
r.rs.setSoilGrid(picker)  # maps segments
r.rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)
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

    rx = r.solve(rs_age + t, -trans * sinusoidal(t), sx[cci], sx, True, wilting_point)  # xylem_flux.py
    x_.append(t)
    y_.append(float(r.collar_flux(rs_age + t, rx, sx)))

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
vp.plot_roots_and_soil(r.rs.mappedSegments(), "pressure head", rx, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)

""" transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential
ax1.plot(x_, -np.array(y_), 'g')  # actual
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
plt.show()
