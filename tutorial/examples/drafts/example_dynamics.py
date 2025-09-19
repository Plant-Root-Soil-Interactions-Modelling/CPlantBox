import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

""" SRA example draft """
import plantbox as pb

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.perirhizal import PerirhizalPython  # Steady rate helper
from functional.root_conductivities import *  # hard coded conductivities
import functional.van_genuchten as vg
from functional.xylem_flux import sinusoidal2

import visualisation.vtk_plot as vp

from rosi_richards import RichardsSP  # RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
import matplotlib.pyplot as plt
import timeit


def set_all_sd(rs, s):
    for p in rs.getOrganRandomParameter(2):
        p.lmaxs = p.lmax * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.rs = p.r * s
        p.a_s = p.a * s


""" Parameters """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 25]  # [16, 16, 30]  # [32, 32, 60]
periodic = True

path = "../../../../CPlantBox/modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -659.8 + 12.5  # -659.8

trans = 6.4  # cm3 /day (sinusoidal)
wilting_point = -15000  # cm
trans_f = lambda t, dt:-trans * sinusoidal2(t, dt)

sim_time = 1  # [day] for task b
rs_age = 10  # root system initial age
age_dependent = False  # conductivities
dt = 360. / (24 * 3600)  # [days] Time step must be very small

N = round(sim_time / dt)
t = 0.
wilting_point = -15000  # cm
skip = 1  # for output and results, skip iteration
max_iter = 10  # maximum for fix point iteration

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
rs = pb.Plant()

rs.readParameters(path + name + ".xml")
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
else:
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, max_b[2] - min_b[2])
rs.setGeometry(sdf)

# set_all_sd(rs, 0.2)

rs.initialize(True)

# picker = lambda x, y, z: s.pick([x, y, z])
# rs.setSoilGrid(picker)  # maps segments
# rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)

rs.simulate(rs_age, True)

# r = XylemFluxPython(rs)
# init_conductivities(r, age_dependent)
# r.test()

for i in range(0, N):

    print("\nIteration", i)
    rs.simulate(dt, True)
    # r.test()

vp.plot_roots(rs, "subType")

print("fin")
