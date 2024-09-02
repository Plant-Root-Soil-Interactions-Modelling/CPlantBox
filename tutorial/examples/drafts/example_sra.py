import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

""" 
SRA example draft  

uses:
Plant, MappedPlant (instead of RootSystem, MappedRootSystem) 
PerirhizalPython class (for Schröder et al. 2008 related things)

for a growing root system
"""
import plantbox as pb

from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.Perirhizal import PerirhizalPython  # Steady rate helper
from functional.root_conductivities import *  # hard coded conductivities
import functional.van_genuchten as vg

import visualisation.vtk_plot as vp

from rosi_richards import RichardsSP  # RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
import matplotlib.pyplot as plt
import timeit

""" Parameters """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [1, 1, 25]  # [16, 16, 30]  # [32, 32, 60]
periodic = True

path = "../../../../CPlantBox/modelparameter/structural/rootsystem/"
name = "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010
loam = [0.08, 0.43, 0.04, 1.6, 50]  # make sure these vg parameters correspond to the look up table
initial = -659.8  # total potential [cm]

trans = 6.4  # [cm3/day] sinusoidal
wilting_point = -15000  # cm
trans_f = lambda t, dt:-trans * HydraulicModel_Doussan.sinusoidal2(t, dt)  # rename

sim_time = 1  # simulation time [day]
initial_age = 10  # initial age [day]
age_dependent = False  # conductivities
dt = 360. / (24 * 3600)  # [days] Time step must be very small

N = round(sim_time / dt)  # number of iterations
t = 0.  # current simulation time [day]
skip = 1  # for output and results, skip iteration
max_iter = 1000  # maximum for fix point iteration

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

""" Initialize plant model """
rs = pb.MappedPlant()
rs.readParameters(path + name + ".xml")
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), -min_b[2])
else:
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, -min_b[2])
rs.setGeometry(sdf)
rs.initialize()
rs.simulate(initial_age, True)

""" Initialize plant hydraulic and perirhizal model"""
params = PlantHydraulicParameters()
init_conductivities(params, age_dependent)

r = HydraulicModel_Doussan(rs, params, cached = True)  # or HydraulicModel_Doussan, HydraulicModel_Meunier
rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False)
picker = lambda x, y, z: s.pick([x, y, z])
r.ms.setSoilGrid(picker)  # maps segment

peri = PerirhizalPython(rs)
peri.open_lookup("../table_loam")

""" Numerical solution """
start_time = timeit.default_timer()

psi_x_, psi_s_, sink_ , t_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
vol_ = [[], [], [], [], [], []]
surf_ = [[], [], [], [], [], []]
krs_ = []
depth_ = []

# initialize loop
collar_ind = r.collar_index()  # segment index of root collar (collar node index == 0)

psi_s_cell = s.getSolutionHead_()  # richards.py
psi_s = np.array(r.ms.getHs(psi_s_cell))

psi_rs = psi_s.copy()  # initial values for fix point iteration

N = int(np.ceil(sim_time / dt))  # number of iterations
print("Starting simulation loop", N, "iterations")

""" simulation loop """
for i in range(0, N):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ growth (and update everything)"""

    rs.simulate(dt, verbose = False)

    outer_r = peri.get_outer_radii("length")  # 0: segment volume, 1: segment surface, 2: segment length
    inner_r = rs.radii
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)
    kr_ = r.get_kr(initial_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

    mapping = rs.getSegmentMapper()

    psi_s_cell = s.getSolutionHead_()  # richards.py
    psi_s = np.array([psi_s_cell[j] for j in mapping])  # soil bulk matric potential per segment
    psi_s = np.maximum(psi_s, np.ones(psi_s.shape) * -15000.)  ############################################ (too keep within table)
    psi_s = np.minimum(psi_s, np.zeros(psi_s.shape))  ############################################ (too keep within table)

    psi_rs = np.hstack((psi_rs, psi_s[psi_rs.shape[0]:]))  # initial values for fix point iteration

    """ fix point iteration """
    wall_fixpoint = timeit.default_timer()

    psi_x = r.solve(initial_age + t, trans_f(initial_age + t, dt), psi_rs, cells = False)
    psi_x_old = psi_x.copy()

    # print("\n")
    # r.test()
    # print("psi_x", np.min(psi_x), np.max(psi_x))
    # print("psi_s", np.min(psi_s), np.max(psi_s))
    # print("\n")

    err_ = 1.e6  # cm
    c = 0
    while err_ > 1 and c < max_iter:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        psi_rs = peri.soil_root_interface_potentials(psi_x[1:], psi_s, inner_kr_, rho_)
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        psi_x = r.solve(initial_age + t, trans_f(initial_age + t, dt), psi_rs, cells = False)  # xylem_flux.py, cells = False
        err_ = np.linalg.norm(psi_x - psi_x_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        c += 1

    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    """ macroscopic soil """
    wall_soil = timeit.default_timer()
    fluxes = r.radial_fluxes(initial_age + t, psi_x, psi_rs)
    collar_flux = r.get_transpiration(initial_age + t, psi_x.copy(), psi_rs.copy())
    err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    if err > 1.e-6:
        print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
    err2 = np.linalg.norm(trans_f(initial_age + t, dt) - collar_flux)
    if r.last == "neumann":
        if err2 > 1.e-6:
            print("error: potential transpiration differs root collar flux in Neumann case" , err2)
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)

    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    sink = np.zeros(psi_s_cell.shape)
    for k, v in soil_fluxes.items():
        sink[k] += v
    t_.append(initial_age + t)  # day
    y_.append(np.sum(sink))  # cm3/day

    if i % skip == 0:

        # if i % (24 * skip) == 0:
        print("time {:g}".format(initial_age + t), "{:g}/{:g} fix point iterations {:g}, {:g}".format(i, N, c, err_),
              "wall times: fix point {:g}:{:g}; soil vs iter {:g}:{:g}".format(wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
                                                                                 wall_fixpoint / wall_iteration, wall_soil / wall_iteration),
              "\nnumber of segments", rs.getNumberOfSegments(), "root collar", psi_x[collar_ind])
        # print("wall_interpolation", wall_interpolation)

        # sink_.append(sink)  # cm3/day (per soil cell)
        # psi_s2_.append(psi_s_cell.copy())  # cm (per soil cell)
        # ana = pb.SegmentAnalyser(r.rs.mappedSegments())  # VOLUME and SURFACE
        # for j in range(0, 6):  # root types
        #     anac = pb.SegmentAnalyser(ana)psi_s
        #     anac.filter("subType", j)
        #     vol_[j].append(anac.getSummed("volume"))
        #     surf_[j].append(anac.getSummed("surface"))
                # depth_.append(ana.getMinBounds().z)
        krs, _ = r.get_krs(initial_age + t)
        krs_.append(krs)  # KRS

        """ direct vtp output """
        # psi_x_.append(psi_x.copy())  # cm (per root node)
        # psi_s_.append(psi_rs.copy())  # cm (per root segment)
        # ana.addData("psi_x", psi_x[1:])
        # ana.addData("psi_rs", psi_rs)
        # ana.addAge(initial_age + t)  # "age"
        # ana.addConductivities(r, initial_age + t)  # "kr", "kx"
        # ana.addFluxes(r, psi_x, psi_rs, initial_age + t)  # "axial_flux", "radial_flux"
        # ana.write("results/rs{0:05d}.vtp".format(int(i / skip)), ["radius", "subType", "creationTime", "organType", "psi_x", "psi_rs", "age", "kr", "kx", "axial_flux", "radial_flux"])

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" save final results """

""" transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(t_, [trans * sinusoidal2(t, dt) for t in t_], 'k')  # potential
ax1.plot(t_, -np.array(y_), 'g')  # actual
ax2 = ax1.twinx()
cum_y = np.cumsum(-np.array(y_) * dt)
ax2.plot(t_, cum_y, 'c--')  # cumulative
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((t_, -np.array(y_))), delimiter = ';')
print("total uptake of", cum_y[-1], "[cm3]")
plt.show()

np.savez("example_sra", time = t_, actual = y_, cumulative = np.cumsum(-np.array(y_) * dt), krs = krs_)

""" VTK visualisation """
vp.plot_roots_and_soil(r.rs, "pressure head", psi_x, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)

