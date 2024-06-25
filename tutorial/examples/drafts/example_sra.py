import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

""" SRA example draft """
import plantbox as pb

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.root_conductivities import *  # hard coded conductivities
import functional.van_genuchten as vg
from functional.xylem_flux import sinusoidal2

import visualisation.vtk_plot as vp

from rosi_richards import RichardsSP  # RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
import matplotlib.pyplot as plt
import timeit
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
from scipy.optimize import ridder
from scipy.optimize import brenth
from scipy.optimize import minimize_scalar


def open_sra_lookup(filename):
    """ opens the look-up table from a file, to quickly find soil root interface potential """
    sra_table = np.load(filename + ".npy")
    x = np.load(filename + "_.npy", allow_pickle = True)
    rx_ = x[0]
    sx_ = x[1]
    inner_ = x[2]
    outer_ = x[3]
    return RegularGridInterpolator((rx_, sx_, inner_, outer_), sra_table)  # method = "nearest" fill_value = None , bounds_error=False


def soil_root_interface(rx, sx, inner_kr, rho, sp):
    """
    finds potential at the soil root interface
    
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    sp             soil van Genuchten parameters (type vg.Parameters)
    """
    k_soilfun = lambda hsoil, hint: (vg.fast_mfp[sp](hsoil) - vg.fast_mfp[sp](hint)) / (hsoil - hint)
    # rho = outer_r / inner_r  # Eqn [5]
    rho2 = rho * rho  # rho squared
    # b = 2 * (rho2 - 1) / (1 + 2 * rho2 * (np.log(rho) - 0.5))  # Eqn [4]
    b = 2 * (rho2 - 1) / (1 - 0.53 * 0.53 * rho2 + 2 * rho2 * (np.log(rho) + np.log(0.53)))  # Eqn [7]
    fun = lambda x, rx, sx, inner_kr, b: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    """ currently best """
    rsx = [fsolve(fun, (rx[i] + sx[i]) / 2, args = (rx[i], sx[i], inner_kr[i], b[i])) for i in range(0, len(rx))]
    """ """
    # rsx = rx
    # for i in range(0, len(rx)):
    #     print(rx[i], sx[i], fun(rx[i], rx[i], sx[i], inner_kr[i], b[i]), fun((9 * sx[i] + rx[i]), rx[i], sx[i], inner_kr[i], b[i]))
    #     x0 = ridder(fun, rx[i], (9 * sx[i] + rx[i]) / 10, args = (rx[i], sx[i], inner_kr[i], b[i]))
    #     print(x0)
    #     rsx[i] = x0
    # rsx = rx

    # for i in range(0, len(rx)):
    #     # print(rx[i], sx[i], fun(rx[i], rx[i], sx[i], inner_kr[i], b[i]), fun((9 * sx[i] + rx[i]), rx[i], sx[i], inner_kr[i], b[i]))
    #     a_ = rx[i]
    #     b_ = (99 * sx[i] + rx[i]) / 100
    #     # x0 = minimize_scalar(fun, bracket = None, bounds = (min(a_, b_), max(a_, b_)), args = (rx[i], sx[i], inner_kr[i], b[i]))
    #     # rsx[i] = x0.x
    #     x0 = brenth(fun, min(a_, b_), max(a_, b_), args = (rx[i], sx[i], inner_kr[i], b[i]))
    #     rsx[i] = x0

    # print("e.g. ind=100 ", rx[100], sx[100], rsx[100])
    return rsx


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
65
66
67
68
69
70
    inner_kr       root radius times hydraulic conductivity [cm/day] 
    rho            geometry factor [1]
    f              function to look up the potentials
    """
    try:
        rsx = f((rx, sx, inner_kr_ , rho_))
    except:
        print("rx", np.min(rx), np.max(rx))  # 0, -16000
        print("sx", np.min(sx), np.max(sx))  # 0, -16000
        print("inner_kr", np.min(inner_kr_), np.max(inner_kr_))  # 1.e-7 - 1.e-4
        print("rho", np.min(rho_), np.max(rho_))  # 1. - 200.
    return rsx


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

# SLOW
soil_vg = vg.Parameters(loam)
vg.create_mfp_lookup(soil_vg, wilting_point = -15000, n = 1501)
root_interface = soil_root_interface  # function defined above

# # FAST
# soil_vg = open_sra_lookup("table_laom")
# root_interface = soil_root_interface_table  # function defined above

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
rs.initialize()
rs.simulate(rs_age, True)
r = XylemFluxPython(rs)
init_conductivities(r, age_dependent)

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])
rs.setSoilGrid(picker)  # maps segments
rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), True)
r.test()  # sanity checks
# cci = picker(nodes[0, 0], nodes[0, 1], nodes[0, 2])  # collar cell index

""" Numerical solution """
start_time = timeit.default_timer()

psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
soil_c_, c_ = [], []
vol_ = [[], [], [], [], [], []]
surf_ = [[], [], [], [], [], []]
krs_ = []
depth_ = []

# rs = r.rs
# nodes = rs.nodes
segs = rs.segments
# ns = len(segs)
mapping = rs.getSegmentMapper()  # because seg2cell is a dict

for i in range(0, len(segs)):
    if segs[i].x == 0:
        collar_ind = i  # segment index of root collar
        print("Collar segment index is ", collar_ind)
        break

sx = s.getSolutionHead_()  # richards.py
hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
rsx = hsb.copy()  # initial values for fix point iteration
rx = r.solve(rs_age, trans_f(0, dt), 0., rsx, False, wilting_point, soil_k = [])
rx_old = rx.copy()

N = int(np.ceil(sim_time / dt))  # number of iterations
print("Starting simulation loop", N, "iterations")

""" simulation loop """
for i in range(0, N):

    wall_iteration = timeit.default_timer()

    t = i * dt  # current simulation time

    """ growth (and update everything)"""
    # rs.simulate(dt, False)

    # cell2seg = rs.cell2seg  # for debugging
    mapping = rs.getSegmentMapper()
    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration

    cell_centers = s.getCellCenters_()
    cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
    seg_centers_z = rs.getSegmentZ()

    outer_r = pb.Perirhizal(rs).segOuterRadii(2)  # 0: segment volume, 1: segment surface, 2: segment length
    inner_r = rs.radii
    types = rs.subTypes
    rho_ = np.divide(outer_r, np.array(inner_r))
    rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

    kr_ = r.getKr(rs_age + t)
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
    inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
    inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

    """ fix point iteration """
    wall_fixpoint = timeit.default_timer()

    err_ = 1.e6  # cm
    c = 0

    # r.init_solve_static(rs_age + t, rsx, False, wilting_point, soil_k = [])  # LU factorisation for speed up
    rx = r.solve(rs_age + t, trans_f(t, dt), 0., rsx, False, wilting_point, soil_k = [])
    rx_old = rx.copy()

    hsb_ = hsb - cell_centers_z  # from total potential to matric potential
    hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
    hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)

    while err_ > 1 and c < max_iter:

        """ interpolation """
        wall_interpolation = timeit.default_timer()
        rx_ = rx[1:] - seg_centers_z  # from total potential to matric potential
        rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)
        rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, soil_vg)
        rsx = rsx + seg_centers_z  # from matric potential to total potential
        wall_interpolation = timeit.default_timer() - wall_interpolation

        """ xylem matric potential """
        wall_xylem = timeit.default_timer()
        # print("Segment size from Python ", len(r.rs.segments), ns)

        rx = r.solve(rs_age + t, trans_f(rs_age + t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
        err_ = np.linalg.norm(rx - rx_old)
        wall_xylem = timeit.default_timer() - wall_xylem

        rx_old = rx.copy()
        c += 1

    wall_fixpoint = timeit.default_timer() - wall_fixpoint

    """ macroscopic soil """
    wall_soil = timeit.default_timer()
    fluxes = r.segFluxes(rs_age + t, rx, rsx, approx = False, cells = False)
    collar_flux = r.collar_flux(rs_age + t, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
    err = np.linalg.norm(np.sum(fluxes) - collar_flux)
    if err > 1.e-6:
        print("error: summed root surface fluxes and root collar flux differ" , err, r.neumann_ind, collar_flux, np.sum(fluxes))
    err2 = np.linalg.norm(trans_f(rs_age + t, dt) - collar_flux)
    if r.last == "neumann":
        if err2 > 1.e-6:
            print("error: potential transpiration differs root collar flux in Neumann case" , err2)
    soil_fluxes = r.sumSegFluxes(fluxes)
    s.setSource(soil_fluxes.copy())  # richards.py
    s.solve(dt)

    wall_soil = timeit.default_timer() - wall_soil

    wall_iteration = timeit.default_timer() - wall_iteration

    """ remember results ... """
    sink = np.zeros(sx.shape)
    for k, v in soil_fluxes.items():
        sink[k] += v
    x_.append(rs_age + t)  # day
    y_.append(np.sum(sink))  # cm3/day

    if i % skip == 0:

        # if i % (24 * skip) == 0:
        print("time {:g}".format(rs_age + t), "{:g}/{:g} fix point iterations {:g}, {:g}".format(i, N, c, err_),
              "\nwall times: fix point {:g}:{:g}; soil vs iter {:g}:{:g}".format(wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
                                                                                 wall_fixpoint / wall_iteration, wall_soil / wall_iteration),
              "number of segments", rs.getNumberOfSegments(), "root collar", rx[0])

        print("wall_interpolation", wall_interpolation)

        sink_.append(sink)  # cm3/day (per soil cell)
        psi_s2_.append(sx.copy())  # cm (per soil cell)
        ana = pb.SegmentAnalyser(r.rs.mappedSegments())  # VOLUME and SURFACE
        for j in range(0, 6):  # root types
            anac = pb.SegmentAnalyser(ana)
            anac.filter("subType", j)
            vol_[j].append(anac.getSummed("volume"))
            surf_[j].append(anac.getSummed("surface"))
        krs, _ = r.get_krs(rs_age + t, [collar_ind])
        krs_.append(krs)  # KRS
        depth_.append(ana.getMinBounds().z)

        """ direct vtp output """
        # psi_x_.append(rx.copy())  # cm (per root node)
        # psi_s_.append(rsx.copy())  # cm (per root segment)
        # ana.addData("rx", rx[1:])
        # ana.addData("rsx", rsx)
        # ana.addAge(rs_age + t)  # "age"
        # ana.addConductivities(r, rs_age + t)  # "kr", "kx"
        # ana.addFluxes(r, rx, rsx, rs_age + t)  # "axial_flux", "radial_flux"
        # ana.write("results/rs{0:05d}.vtp".format(int(i / skip)), ["radius", "subType", "creationTime", "organType", "rx", "rsx", "age", "kr", "kx", "axial_flux", "radial_flux"])

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

""" transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, -trans * sinusoidal2(t, dt), 'k')  # potential
ax1.plot(x_, -np.array(y_), 'g')  # actual
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[cm^3 d^{-1}]$")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
plt.show()

""" VTK visualisation """
print(rx.shape)
vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)

