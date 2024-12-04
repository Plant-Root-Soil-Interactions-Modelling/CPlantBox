import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

""" 
SRA example draft  

uses:

Plant, MappedPlant (instead of RootSystem, MappedRootSystem) 
PerirhizalPython class (for Schröder et al. 2008 related things)

PlandHydraulicModel and PlantHydraulicParameters (instead of xylem_flux)
"""

import plantbox as pb

# from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # Doussan solver
from functional.PlantHydraulicModel import PlantHydraulicModel  # Doussan solver

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

sim_time = 7  # final simulation time [day]
initial_age = 10  # initial age [day]
age_dependent = True  # conductivities
dt = 360. / (24 * 3600)  # [days]

t = 0.  # current simulation time [day]
N = round(sim_time / dt)
wilting_point = -15000  # cm
skip = 1  # for output and results, skip iteration
max_iter = 1000  # maximum for fix point iteration

# SLOW
soil_vg = vg.Parameters(loam)
vg.create_mfp_lookup(soil_vg, wilting_point = -15000, n = 15001)
peri = PerirhizalPython()
peri.open_lookup("table_loam")

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
rs = pb.MappedPlant()  # pb.MappedRootSystem()  #
rs.readParameters(path + name + ".xml")
if not periodic:
    sdf = pb.SDF_PlantBox(0.99 * (max_b[0] - min_b[0]), 0.99 * (max_b[1] - min_b[1]), max_b[2] - min_b[2])
else:
    sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2])
rs.setGeometry(sdf)
rs.initialize()
rs.simulate(initial_age, True)

params = PlantHydraulicParameters()
init_conductivities(params, age_dependent)
r = PlantHydraulicModel("Doussan", rs, params)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
rs.setSoilGrid(picker)  # maps segments
rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False)  # DISABLE cutting....
r.test()  # sanity checks

""" Numerical solution """
start_time = timeit.default_timer()

psi_x_, psi_s_, sink_ , t_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
vol_ = [[], [], [], [], [], []]
surf_ = [[], [], [], [], [], []]
krs_ = []
depth_ = []

collar_ind = r.collar_index()  # segment index of root collar

# initialize loop
mapping = rs.getSegmentMapper()  # because seg2cell is a dict
sx = s.getSolutionHead_()  # richards.py
hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
rsx = hsb.copy()  # initial values for fix point iteration
r.update(initial_age)
rx = r.solve(rsx, trans_f(0, dt), wilting_point)
rx_old = rx.copy()

N = int(np.ceil(sim_time / dt))  # number of iterations
print("Starting simulation loop", N, "iterations")

""" 
    Simulation loop 
"""
for i in range(0, N):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ growth (and update everything)"""

    rs.simulate(dt, False)

    outer_r = PerirhizalPython(rs).seg_outer_radii("length")
    inner_r = rs.radii
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

    mapping = rs.getSegmentMapper()

    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
    hsb_ = hsb  # - cell_centers_z  # from total potential to matric potential
    hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
    hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)

    rsx = np.hstack((rsx, hsb[rsx.shape[0]:]))  # initial values for fix point iteration, last value and hsb for new segs

    # cell_centers = s.getCellCenters_()
    # cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
    seg_centers_z = rs.getSegmentZ()

    kr_ = r.get_kr(initial_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

    """ fix point iteration """
    wall_fixpoint = timeit.default_timer()

    err_ = 1.e6  # cm
    c = 0

    # r.init_solve_static(initial_age + t, rsx, False, wilting_point, soil_k = [])  # LU factorisation for speed up
    r.update(initial_age + t)
    rx = r.solve(rsx, trans_f(initial_age + t, dt), wilting_point)
    rx_old = rx.copy()

    while err_ > 1 and c < max_iter:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        rx_ = rx - seg_centers_z  # from total potential to matric potential
        rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)
        rsx = peri.soil_root_interface_potentials(rx_ , hsb_, inner_kr_, rho_, soil_vg)
        rsx = rsx + seg_centers_z  # from matric potential to total potential
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        # print("Segment size from Python ", len(r.rs.segments), ns)
        rx = r.solve(rsx, trans_f(initial_age + t, dt), wilting_point)  # xylem_flux.py, cells = False
        # print(r.last)
        err_ = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        rx_old = rx.copy()
        c += 1

    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    """ macroscopic soil """
    wall_soil = timeit.default_timer()
    fluxes = r.radial_fluxes(rx, rsx)
    # collar_flux = r.collar_flux(initial_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
    # err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    # if err > 1.e-6:
    #     print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
    # err2 = np.linalg.norm(trans_f(initial_age + t, dt) - collar_flux)
    # if r.last == "neumann":
    #     if err2 > 1.e-6:
    #         print("error: potential transpiration differs root collar flux in Neumann case" , err2)
    print("fluxes", fluxes.shape)
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)

    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    sink = np.zeros(sx.shape)
    for k, v in soil_fluxes.items():
        sink[k] += v
    t_.append(initial_age + t)  # day
    y_.append(np.sum(sink))  # cm3/day

    if i % skip == 0:

        # if i % (24 * skip) == 0:
        print("time {:g}".format(initial_age + t), "{:g}/{:g} fix point iterations {:g}, {:g}".format(i, N, c, err_),
              "wall times: fzero vs. xylem {:g}:{:g}; fixpoint vs. soil {:g}:{:g}".format(wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
                                                                                 wall_fixpoint / wall_iteration, wall_soil / wall_iteration),
              "\nnumber of segments", rs.getNumberOfSegments(), "root collar {:g} cm\n".format(rx[collar_ind]))
        # print("wall_interpolation", wall_interpolation)

        # sink_.append(sink)  # cm3/day (per soil cell)
        # psi_s2_.append(sx.copy())  # cm (per soil cell)
        # ana = pb.SegmentAnalyser(r.rs.mappedSegments())  # VOLUME and SURFACE
        # for j in range(0, 6):  # root types
        #     anac = pb.SegmentAnalyser(ana)
        #     anac.filter("subType", j)
        #     vol_[j].append(anac.getSummed("volume"))
        #     surf_[j].append(anac.getSummed("surface"))
        # krs, _ = r.get_krs(initial_age + t, [collar_ind])
        # krs_.append(krs)  # KRS
        # depth_.append(ana.getMinBounds().z)

        """ direct vtp output """
        # psi_x_.append(rx.copy())  # cm (per root node)
        # psi_s_.append(rsx.copy())  # cm (per root segment)
        # ana.addData("rx", rx[1:])
        # ana.addData("rsx", rsx)
        # ana.addAge(initial_age + t)  # "age"
        # ana.addConductivities(r, initial_age + t)  # "kr", "kx"
        # ana.addFluxes(r, rx, rsx, initial_age + t)  # "axial_flux", "radial_flux"
        # ana.write("results/rs{0:05d}.vtp".format(int(i / skip)), ["radius", "subType", "creationTime", "organType", "rx", "rsx", "age", "kr", "kx", "axial_flux", "radial_flux"])

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" save final results """

""" transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(t_, [trans * sinusoidal2(t, dt) for t in t_], 'k')  # potential
ax1.plot(t_, -np.array(y_), 'g')  # actual
ax2 = ax1.twinx()
ax2.plot(t_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((t_, -np.array(y_))), delimiter = ';')
plt.show()

np.savez("example_sra", time = t_, actual = y_, cumulative = np.cumsum(-np.array(y_) * dt), krs = krs_)

""" VTK visualisation """
print(rx.shape)
vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)

