import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../../dumux-rosi/python/modules/")  # python wrappers

""" CPLantBox tutorial example 6c (see CPlantBox/tutorial/latex/PlantBox_RootSystem/tutorial.tex) """
""" coupling with DuMux as solver for the soil, run in dumux-rosi """

from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
from functional.root_conductivities import *  # hard coded conductivities
import plantbox as pb
import visualisation.vtk_plot as vp
from rosi_richards import RichardsSP  # RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part

import numpy as np
import matplotlib.pyplot as plt
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


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
    fun = lambda x: (inner_kr * rx + b * sx * k_soilfun(sx, x)) / (b * k_soilfun(sx, x) + inner_kr) - x
    rsx = fsolve(fun, rx)
    return rsx


def soil_root_interface_table(rx, sx, inner_kr_, rho_, f):
    """
    finds potential at the soil root interface
        
    rx             xylem matric potential [cm]
    sx             bulk soil matric potential [cm]
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
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, max_b[2] - min_b[2])
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

""" VTK visualisation """
print(rx.shape)
vp.plot_roots_and_soil(r.rs, "pressure head", rx, s, periodic, np.array(min_b), np.array(max_b), cell_number, name)

"""
functions for the steady rate approach
"""
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import fsolve
import timeit

import vtk_plot as vtk
import plantbox as pb
import van_genuchten as vg
import evapotranspiration as evap


def simulate_dynamic(s, r, sra_table_lookup, sim_time, dt, trans_f, rs_age = 1., type_ = 1):
    """     
    simulates the coupled scenario       
        root architecture is not gowing  
        conductivities are not changing over time
        
    s                            soil model (RichardsWrapper(RichardsSP()))
    r                            xylem flux model (XylemFluxPython wrapping MappedSegments mapped to soil @param s)
    sra_table_lookup             potentials a root soil interface    
    sim_time                     simulation time
    dt                           time step
    trans_f                      potential transpiration function 
    rs_age                       initial root system age  
    type_                        1 = water only, 2 = water and nitrate
    """

    wilting_point = -15000  # cm
    skip = 10  # for output and results, skip iteration
    max_iter = 10  # maximum for fix point iteration

    """ tabularized values for finding the zeros """
    if isinstance(sra_table_lookup, RegularGridInterpolator):
        root_interface = soil_root_interface_table  # function defined above
    else:
        raise
        root_interface = soil_root_interface  # function defined above

    start_time = timeit.default_timer()

    psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    soil_c_, c_ = [], []
    vol_ = [[], [], [], [], [], []]
    surf_ = [[], [], [], [], [], []]
    krs_ = []
    depth_ = []

    rs = r.rs
    nodes = rs.nodes
    segs = rs.segments
    ns = len(segs)
    mapping = rs.getSegmentMapper()  # because seg2cell is a dict

    for i in range(0, len(segs)):
        if segs[i].x == 0:
            collar_ind = i  # segment index of root collar
            break

    sx = s.getSolutionHead_()  # richards.py
    hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
    rsx = hsb.copy()  # initial values for fix point iteration
    rx = r.solve(rs_age, trans_f(0, dt), 0., rsx, False, wilting_point, soil_k = [])
    rx_old = rx.copy()

    N = int(np.ceil(sim_time / dt))  # number of iterations

    print("Starting simulation loop")

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ grow root system and update everything"""
        rs.simulate(dt, False)

        cell2seg = rs.cell2seg  # for debugging
        mapping = rs.getSegmentMapper()
        sx = s.getSolutionHead_()  # richards.py
        hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
        rsx = hsb.copy()  # initial values for fix point iteration

        cell_centers = s.getCellCenters_()
        cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
        seg_centers_z = rs.getSegmentZ()

        outer_r = r.rs.segOuterRadii()
        inner_r = r.rs.radii
        types = r.rs.subTypes
        rho_ = np.divide(outer_r, np.array(inner_r))
        rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)

        kr_ = r.getKr(rs_age + t)
        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
        inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
        inner_kr_ = np.minimum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-4)  ############################################ (too keep within table)

        wall_iteration = timeit.default_timer()
        wall_fixpoint = timeit.default_timer()

        err = 1.e6  # cm
        c = 0

        r.init_solve_static(rs_age + t, rsx, False, wilting_point, soil_k = [])  # LU factorisation for speed up
        rx = r.solve(rs_age + t, trans_f(t, dt), 0., rsx, False, wilting_point, soil_k = [])
        rx_old = rx.copy()

        hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
        hsb_ = np.maximum(hsb_, np.ones(hsb_.shape) * -15000.)  ############################################ (too keep within table)
        hsb_ = np.minimum(hsb_, np.zeros(hsb_.shape))  ############################################ (too keep within table)

        while err > 1 and c < max_iter:

            """ interpolation """
            wall_interpolation = timeit.default_timer()
            rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
            rx_ = np.maximum(rx_, np.ones(rx_.shape) * -15000.)  ############################################ (too keep within table)
            rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)
            rsx = rsx + seg_centers_z  # from matric potential to total matric potential
            wall_interpolation = timeit.default_timer() - wall_interpolation

            """ xylem matric potential """
            wall_xylem = timeit.default_timer()
            # print("Segment size from Python ", len(r.rs.segments), ns)
            rx = r.solve(rs_age + t, trans_f(rs_age + t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
            err = np.linalg.norm(rx - rx_old)
            wall_xylem = timeit.default_timer() - wall_xylem

            rx_old = rx.copy()
            c += 1

        wall_fixpoint = timeit.default_timer() - wall_fixpoint

        if type_ == 2:
            cc = s.getSolution_(1)  # kg/m3
            rsc = np.array([cc[i] for i in mapping])  # kg/m3
            seg_sol_fluxes = np.array(r.solute_fluxes(rsc))  # [g/day]
            soil_sol_fluxes = r.sumSegFluxes(seg_sol_fluxes)  # [g/day]
            # evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.)  # = 1.14e-4 g/day
            evap.add_nitrificatin_source(s, soil_sol_fluxes, nit_flux = 0.)  # = 1.14e-4 g/day 1.e-7 * (75 * 16 * 1)
            s.setSource(soil_sol_fluxes.copy(), eq_idx = 1)  # [g/day], in moduels/richards.py

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

        # for key, value in cell2seg.items():  # check cell2seg
        #     if key < 0:
        #         nodes = r.rs.nodes
        #         print("key is negative", key)
        #         print("segments", cell2seg[key])
        #         print("coresponding nodes")
        #         segs = r.rs.segments
        #         for s in cell2seg[key]:
        #             print(segs[s])
        #             print(nodes[segs[s].x], nodes[segs[s].y])
        #         ana = pb.SegmentAnalyser(r.rs.mappedSegments())
        #         ana.addCellIds(r.rs.mappedSegments())
        #         vtk.plot_roots(ana, "cell_id")

        wall_soil = timeit.default_timer() - wall_soil

        wall_iteration = timeit.default_timer() - wall_iteration

        # if rs_age + t > 24.5:
        #     pass
        #     min_b = [-19, -2.5, -200.]  # for soybean
        #     max_b = [19, 2.5, 0.]
        #     cell_number = [1, 1, 200]
        #     vtk.plot_roots_and_soil(rs, "fluxes", fluxes.copy()[1:], s, True, min_b, max_b, cell_number, "nice_plot")
        #     # vtk.plot_roots(pd, p_name:str, win_title:str = "", render:bool = True):
        #     # ind0 = s.pick([0, 0, -3.5])
        #     # ind1 = s.pick([0, 0, -15.])
        #     # ind2 = s.pick([0, 0, -25.])
        #     # print("cell0", ind0)
        #     # print("cell1", ind1)
        #     # print("cell2", ind2)
        #     # cell2seg = r.rs.cell2seg
        #     # segs0 = cell2seg[ind0]
        #     # # segs1 = cell2seg[ind1]
        #     # # segs2 = cell2seg[ind2]
        #     # for i in segs:
        #     #     rs.plot_cylinder(i)
        #     dd

        """ remember results ... """
        sink = np.zeros(sx.shape)
        for k, v in soil_fluxes.items():
            sink[k] += v
        x_.append(rs_age + t)  # day
        y_.append(np.sum(sink))  # cm3/day
        if type_ == 2:
            c_.append(-np.sum(seg_sol_fluxes))  # [g/day]

        if i % skip == 0:

            # if i % (24 * skip) == 0:
            print("time", rs_age + t, "{:g}/{:g} {:g} iterations".format(i, N, c), "wall times",
                  wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
                  "number of segments", rs.getNumberOfSegments(), "root collar", rx[0])

            sink_.append(sink)  # cm3/day (per soil cell)

            psi_s2_.append(sx.copy())  # cm (per soil cell)

            if type_ == 2:
                soil_c_.append(cc)  # [kg/m3]

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

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_, vol_, surf_, krs_, depth_, soil_c_, c_

# def simulate_const(s, r, sra_table_lookup, trans, sim_time, dt, trans_f = None):
#     """
#     simulates the coupled scenario
#         root architecture is not gowing
#         conductivities are not changing over time
#
#     s                            soil model (RichardsWrapper(RichardsSP()))
#     r                            xylem flux model (XylemFluxPython wrapping MappedSegments mapped to soil @param s)
#     sra_table_lookup             potentials a root soil interface
#     trans                        daily transpiration
#     sim_time                     simulation time
#     dt                           time step
#
#     TODO recyle factorisation of left hand side ...
#     """
#     wilting_point = -15000  # cm
#     skip = 10  # for output and results, skip iteration
#     rs_age = 0.  # day
#     max_iter = 100  # maximum for fix point iteration
#
#     """ tabularized values for finding the zeros """
#     if isinstance(sra_table_lookup, RegularGridInterpolator):
#         root_interface = soil_root_interface_table
#     else:
#         root_interface = soil_root_interface
#
#     """ set defalut potential transpiration """
#     if not trans_f:
#         trans_f = lambda age, dt:-trans * sinusoidal2(age, dt)
#
#     start_time = timeit.default_timer()
#
#     nodes = r.rs.nodes
#     segs = r.rs.segments
#     ns = len(segs)
#     map = r.rs.seg2cell
#     mapping = np.array([map[j] for j in range(0, ns)])  # because seg2cell is a map
#     cell_centers = s.getCellCenters()
#     cell_centers_z = np.array([cell_centers[j][2] for j in mapping])
#     seg_centers_z = np.array([0.5 * (nodes[s.x].z + nodes[s.y].z) for s in segs])
#
#     outer_r = r.rs.segOuterRadii()
#     inner_r = np.array(r.rs.radii)
#     types = r.rs.subTypes
#     rho_ = np.divide(outer_r, inner_r)
#     print("rho", np.min(rho_), np.max(rho_))
#     rho_ = np.minimum(rho_, np.ones(rho_.shape) * 200)  ############################################ (too keep within table)
#
#     psi_x_, psi_s_, sink_ , x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
#
#     sx = s.getSolutionHead()  # inital condition, solverbase.py
#     hsb = np.array([sx[j][0] for j in mapping])  # soil bulk matric potential per segment
#     rsx = hsb.copy()  # initial values for fix point iteration
#
#     r.init_solve_static(rs_age, rsx, False, wilting_point, soil_k = [])  # speed up & and forever static... ########################
#
#     rx = r.solve(rs_age, trans_f(0, dt), 0., rsx, False, wilting_point, soil_k = [])
#     rx_old = rx.copy()
#
#     kr_ = r.getKr(rs_age + t)
#     inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const
#     inner_kr_ = np.maximum(inner_kr_, np.ones(inner_kr_.shape) * 1.e-7)  ############################################ (too keep within table)
#
#     N = int(np.ceil(sim_time / dt))  # number of iterations
#
#     """ simulation loop """
#     for i in range(0, N):
#
#         t = i * dt  # current simulation time
#
#         wall_iteration = timeit.default_timer()
#         wall_fixpoint = timeit.default_timer()
#
#         err = 1.e6  # cm
#         c = 0
#         while err > 1 and c < max_iter:
#
#             """ interpolation """
#             wall_interpolation = timeit.default_timer()
#             rx_ = rx[1:] - seg_centers_z  # from total matric potential to matric potential
#             hsb_ = hsb - cell_centers_z  # from total matric potential to matric potential
#             rsx = root_interface(rx_ , hsb_, inner_kr_, rho_, sra_table_lookup)
#             rsx = rsx + seg_centers_z  # from matric potential to total matric potential
#             wall_interpolation = timeit.default_timer() - wall_interpolation
#
#             """ xylem matric potential """
#             wall_xylem = timeit.default_timer()
#             rx = r.solve(rs_age, trans_f(t, dt), 0., rsx, False, wilting_point, soil_k = [])  # xylem_flux.py, cells = False
#             err = np.linalg.norm(rx - rx_old)
#             wall_xylem = timeit.default_timer() - wall_xylem
#
#             rx_old = rx.copy()
#             c += 1
#
#         wall_fixpoint = timeit.default_timer() - wall_fixpoint
#
#         wall_soil = timeit.default_timer()
#         fluxes = r.segFluxes(rs_age, rx, rsx, approx = False, cells = False)
#         collar_flux = r.collar_flux(rs_age, rx.copy(), rsx.copy(), k_soil = [], cells = False)  # validity checks
#         err = np.linalg.norm(np.sum(fluxes) - collar_flux)
#         if err > 1.e-6:
#             print("error: summed root surface fluxes and root collar flux differ" , err)
#         err2 = np.linalg.norm(trans_f(t, dt) - collar_flux)
#         if r.last == "neumann":
#             if err2 > 1.e-6:
#                 print("error: potential transpiration differs root collar flux in Neumann case" , err2)
#         soil_fluxes = r.sumSegFluxes(fluxes)
#         s.setSource(soil_fluxes.copy())  # richards.py
#         s.solve(dt)
#         sx = s.getSolutionHead()[:, 0]  # richards.py
#         hsb = np.array([sx[j] for j in mapping])  # soil bulk matric potential per segment
#         wall_soil = timeit.default_timer() - wall_soil
#
#         wall_iteration = timeit.default_timer() - wall_iteration
#
#         """ remember results ... """
#         print("{:g}/{:g} {:g} iterations".format(i, N, c), "wall times",
#               wall_interpolation / (wall_interpolation + wall_xylem), wall_xylem / (wall_interpolation + wall_xylem),
#               "number of segments", rs.getNumberOfSegments())
#         sink = np.zeros(sx.shape)
#         for k, v in soil_fluxes.items():
#             sink[k] += v
#         sink_.append(sink)  # cm3/day (per soil cell)
#         x_.append(t)  # day
#         y_.append(np.sum(sink))  # cm3/day
#         psi_s2_.append(sx.copy())  # cm (per soil cell)
#         if i % skip == 0:
#             psi_x_.append(rx.copy())  # cm (per root node)
#             psi_s_.append(rsx.copy())  # cm (per root segment)
#
#     print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")
#
#     return psi_x_, psi_s_, sink_, x_, y_, psi_s2_

