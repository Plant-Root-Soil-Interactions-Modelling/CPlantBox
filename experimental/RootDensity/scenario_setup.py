"""
Functions to simplify setup of the scenarios for the INARI project
"""

import timeit
from datetime import *

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plantbox as pb  # CPlantBox
import plantbox.functional.van_genuchten as vg
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
from rosi.richards import RichardsWrapper  # Python part, macroscopic soil model
from rosi.rosi_richards import (
    RichardsSP,  # C++ part (Dumux binding), macroscopic soil model
)
from rosi.rosi_richardsnc import (
    RichardsNCSP,  # C++ part (Dumux binding), macroscopic soil model
)

SMALL_SIZE = 16
MEDIUM_SIZE = 16
BIGGER_SIZE = 16
plt.rc("font", size=SMALL_SIZE)  # controls default text sizes
plt.rc("axes", titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc("axes", labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("ytick", labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc("legend", fontsize=SMALL_SIZE)  # legend fontsize
plt.rc("figure", titlesize=BIGGER_SIZE)  # fontsize of the figure title
prop_cycle = plt.rcParams["axes.prop_cycle"]
colors = prop_cycle.by_key()["color"]


def vg_enviro_type(i: int):
    """Van Genuchten parameter for enviro-types, called by maize() and soybean()"""
    soil = {}
    soil[0] = [0.0639, 0.3698, 0.0096, 1.4646, 4.47]
    soil[1] = [0.0619, 0.3417, 0.0132, 1.3258, 2.03]
    soil[36] = [0.0760, 0.3863, 0.0091, 1.4430, 2.99]
    soil[5] = [0.0451, 0.3721, 0.0325, 1.4393, 34.03]
    soil[59] = [0.0534, 0.3744, 0.0171, 1.4138, 13.09]
    table_name = "envirotype{:s}".format(str(i))
    return soil[i], table_name


def springbarley(i: int):
    """parameters for maize simulation"""
    soil, table_name = vg_enviro_type(i)
    min_b = np.array([-6.5, -1.5, -150.0])
    max_b = np.array([6.5, 1.5, 0.0])
    cell_number = np.array([1, 1, 150])
    width = max_b - min_b
    area = width[0] * width[1]  # cm2
    Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998
    return soil, table_name, min_b, max_b, cell_number, area, Kc_maize


def maize(i: int):
    """parameters for maize simulation"""
    soil, table_name = vg_enviro_type(i)
    min_b = np.array([-38.0, -8.0, -200.0])  # data from INARI
    max_b = np.array([38.0, 8.0, 0.0])
    cell_number = np.array([1, 1, 200])
    width = max_b - min_b
    area = width[0] * width[1]  # cm2
    Kc_maize = 1.2  # book "crop evapotranspiration" Allen, et al 1998
    return soil, table_name, min_b, max_b, cell_number, area, Kc_maize


def soybean(i: int):
    """parameters for soybean simulation"""
    soil, table_name = vg_enviro_type(i)
    min_b = np.array([-38, -2, -200.0])  # data from INARI
    max_b = np.array([38, 2, 0.0])
    cell_number = np.array([1, 1, 200])
    width = max_b - min_b
    area = width[0] * width[1]  # cm2
    Kc_soybean = 1.15  # book "crop evapotranspiration" Allen, et al 1998
    return soil, table_name, min_b, max_b, cell_number, area, Kc_soybean


def create_soil_model(soil_, min_b, max_b, cell_number, type, times=None, net_inf=None, wet=False):
    """
    Creates a soil domain from @param min_b to @param max_b with resolution @param cell_number
    soil type is fixed and homogeneous
    domain is periodic (if 2d or 3d)
    initial potentials are linear from @param p_top to @param p_bot

    returns soil_model (RichardsWrapper(RichardsSP())) and soil parameter (vg.Parameters)
    """
    soil = vg.Parameters(soil_)
    vg.create_mfp_lookup(soil, -1.0e5, 1000)
    if type == 1:
        s = RichardsWrapper(RichardsSP())  # water only
    elif type == 2:
        s = RichardsWrapper(RichardsNCSP())  # water and one solute
    else:
        print("choose type, 1 = Richards, 2 = RichardsNCSP")
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, True)  # [cm]

    # Initial conditions for fertilization
    if type == 2:  # solute IC
        z_ = [0.0, -30.0, -30.0, -200.0]
        v_ = np.array([2.6e-4, 2.6e-4, 0.75 * 2.6e-4, 0.75 * 2.6e-4])  # kg / m3 (~4.e-4)
        # [1.5 * 2.6e-4, 1.5 * 2.6e-4, 2.6e-4, 2.6e-4]  # kg/m3 [2.e-4, 2.e-4, 1.e-4, 1.e-4]  # TODO [0., 0., 0., 0.]  #
        # -> Fletcher et al. 2021 initial solution concentration = 0.43 mol/m3 (2.6e-4 = 0.43*62*1e-3) (nitrate 62 g/mol)
        s.setICZ_solute(v_[::-1], z_[::-1])  # ascending order...

    # BC
    if times is not None:
        if wet:
            net_inf[net_inf > 0] = net_inf[net_inf > 0] * 1.2  # increase precipitation for 20%
            net_inf[net_inf < 0] = net_inf[net_inf < 0] * 0.8  # decrease evaporation for 20%
        s.setTopBC("atmospheric", 0.5, [times, net_inf])  # 0.5 is dummy value
    else:
        s.setTopBC("noFlux")
    s.setBotBC("noFlux")

    if type == 2:  # solute BC
        #  90 lb/ha = 40.8 kg/ha -> 4.08 g /m2 *1.e-4 -> 4.08e-4 g/cm2
        f0 = 4.08e-4  # g/cm2 (90)
        f1 = 2.27e-5  # g/cm2 (5)
        f2 = 5.43e-4  # g/cm2 (120)
        sol_times = np.array([0.0, 1.0, 1.0, 50.0, 50.0, 51.0, 51.0, 1.0e3])
        sol_influx = -np.array([f1, f1, 0.0, 0.0, f2, f2, 0.0, 0.0])  # g/(cm2 day)
        s.setTopBC_solute("managed", 0.5, [sol_times, sol_influx])
        s.setBotBC_solute("outflow", 0.0)

    # Paramters
    s.setVGParameters([soil_])
    s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    if type == 2:
        s.setParameter("Component.MolarMass", "6.2e-2")  # nitrate 62,0049 g/mol
        s.setParameter("Component.LiquidDiffusionCoefficient", "1.7e-9")  # m2 s-1 # nitrate = 1700 um^2/sec
    s.initializeProblem()
    wilting_point = -15000
    s.setCriticalPressure(wilting_point)  # for boundary conditions constantFlow, constantFlowCyl, and atmospheric
    s.ddt = 1.0e-5  # [day] initial Dumux time step
    # print("A")
    # # IC
    # h = np.load("data/initial_potential.npy")
    # s.setInitialConditionHead(h)  # cm
    # print("B")
    # if type == 2:
    #     c = np.load("data/initial_concentration.npy")  # kg/m3
    #     s.setInitialCondition(c, 1)  # kg/m3

    # plt.plot(h, np.linspace(-200., 0., h.shape[0]))
    # plt.xlabel("soil matric potential [cm]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()
    # plt.plot(c, np.linspace(-200, 0., c.shape[0]))
    # plt.xlabel("nitrate concentration [g/cm3]")
    # plt.ylabel("depth (cm)")
    # plt.tight_layout()
    # plt.show()

    return s, soil


def create_mapped_singleroot(min_b, max_b, cell_number, soil_model, ns=100, l=50, a=0.05):
    """creates a single root mapped to a soil with @param ns segments, length l, and radius a"""
    global picker  # make sure it is not garbage collected away...
    r = create_singleroot(ns, l, a)
    r.rs.setRectangularGrid(
        pb.Vector3d(min_b[0], min_b[1], min_b[2]),
        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]),
        cut=False,
    )
    picker = lambda x, y, z: soil_model.pick(
        [x, y, z]
    )  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.rs.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    return r


def create_singleroot(ns=100, l=50, a=0.05):
    """creates a single root with @param ns segments, length l, and radius a"""
    radii = np.array([a] * ns)
    nodes = [pb.Vector3d(0, 0, 0)]
    segs = []
    dx = l / ns
    z_ = np.linspace(-dx, -l, ns)
    for i in range(0, ns):
        nodes.append(pb.Vector3d(0, 0, z_[i]))
        segs.append(pb.Vector2i(i, i + 1))
    rs = pb.MappedSegments(nodes, segs, radii)
    return XylemFluxPython(rs)


def set_all_sd(rs, s):
    """# sets all standard deviation to a percantage, i.e. value*s"""
    for p in rs.getOrganRandomParameter(pb.OrganTypes.root):
        p.a_s = p.a * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.lmaxs = p.lmax * s
        p.rs = p.r * s
        p.thetas = p.theta * s
        p.rlts = p.rlt * s  # no used
        p.ldelays = p.ldelay * s
    seed = rs.getOrganRandomParameter(pb.OrganTypes.seed, 0)
    seed.firstBs = seed.firstB * s
    seed.delayBs = seed.delayB * s
    seed.maxBs = seed.maxB * s
    seed.firstSBs = seed.firstSB * s
    seed.delaySBs = seed.delaySB * s
    seed.delayRCs = seed.delayRC * s
    seed.nCs = seed.nCs * s
    seed.nzs = seed.nzs * s
    # todo seed position s


def create_mapped_rootsystem(min_b, max_b, cell_number, soil_model, fname, stochastic=False, mods=None):
    """loads a rmsl file, or creates a rootsystem opening an xml parameter set,
    and maps it to the soil_model"""
    global picker  # make sure it is not garbage collected away...

    if fname.endswith(".rsml"):
        r = XylemFluxPython(fname)
    elif fname.endswith(".xml"):
        # if rank == 0:
        #     if stochastic:
        #         seed = np.random.randint(0, 1e6)
        #     else:
        #         seed = 1  # always the same random seed
        # else:
        #     seed = None
        # seed = comm.bcast(seed, root = 0)  # random seed must be the same for each process
        seed = 1

        rs = pb.MappedPlant()
        rs.enableExtraNode()
        rs.setSeed(seed)
        rs.readParameters(fname)

    if fname == "data/Zeamays_synMRI_modified.xml":
        print("\nMaize\n")
        params = rs.getOrganRandomParameter(pb.OrganTypes.root)
        for p in params:
            p.a = 2.0 * p.a  # at least ..... TODO parameterisation

    if fname == "data/spring_barley_CF12.xml":
        print("\nSpring barley\n")
        params = rs.getOrganRandomParameter(pb.OrganTypes.root)
        params[2].lmax *= 2
        params[1].theta = 1.31  # why is the tap root not always 0?

    if not stochastic:
        set_all_sd(rs, 0.0)

    if mods is not None:  # apply modifications
        rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
        srp = rs.getOrganRandomParameter(pb.OrganTypes.seed)
        if "lmax" in mods:  # all types
            for i in range(0, len(rrp)):
                rrp[i].lmax *= mods["lmax"]
            mods.pop("lmax")
        if "lmax145" in mods:
            rrp[1].lmax *= mods["lmax145"]
            rrp[4].lmax *= mods["lmax145"]
            if len(rrp) > 5:
                rrp[5].lmax *= mods["lmax145"]
            mods.pop("lmax145")
        if "lmax2" in mods:
            rrp[2].lmax *= mods["lmax2"]
            mods.pop("lmax2")
        if "lmax3" in mods:
            rrp[3].lmax *= mods["lmax3"]
            mods.pop("lmax3")
        if "lmax4" in mods:
            rrp[4].lmax *= mods["lmax4"]
            mods.pop("lmax4")
        if "theta45" in mods:
            if len(rrp) > 5:
                print("shootbore (theta45)")
                rrp[5].theta = mods["theta45"]
            else:
                print("seminal (theta45)")
                rrp[4].theta = mods["theta45"]
            mods.pop("theta45")
        if "r145" in mods:
            rrp[1].r *= mods["r145"]
            rrp[4].r *= mods["r145"]
            if len(rrp) > 5:
                rrp[5].r *= mods["r145"]
            mods.pop("r145")
        if "r2" in mods:
            rrp[2].r *= mods["r2"]
            mods.pop("r2")
        if "r3" in mods:
            rrp[3].r *= mods["r3"]
            mods.pop("r3")
        if "r" in mods:  # all types
            for i in range(0, len(rrp)):
                rrp[i].r *= mods["r"]
            mods.pop("r")
        if "a" in mods:  # all types
            for i in range(0, len(rrp)):
                rrp[i].a *= mods["a"]
            mods.pop("a")
        if "src" in mods:
            srp[0].maxB = mods["src"]
            mods.pop("src")
        if "delaySB" in mods:
            srp[0].delaySB = mods["delaySB"]
            mods.pop("delaySB")
        if mods:  # something unused in mods
            print("\nscenario_setup.create_mapped_rootsystem() WARNING mods have unused parameters:")
            for k, v in mods.items():
                print("key:", k)
            print()

    # rrp = rs.getOrganRandomParameter(pb.OrganTypes.root)
    # for p in rrp:
    #     print(p.dx, p.dxMin)

    rs.setGeometry(pb.SDF_PlantBox(max_b[0] - min_b[0], max_b[1] - min_b[1], np.abs(min_b[2])))
    # rs.setGeometry(pb.SDF_PlantBox(1.0e6, 1.0e6, np.abs(min_b[2])))
    # rs.initializeDB(4, 5)
    rs.initialize()
    rs.simulate(1.0, True)
    r = HydraulicModel_Doussan(rs, PlantHydraulicParameters())

    # print("HERE***********************************")
    # print([s.x for s in r.rs.segments])
    # print([s.y for s in r.rs.segments])
    # # for i in range(0, len(r.rs.segments)): # ????????
    # #     print(r.rs.seg2cell[i])

    r.ms.setRectangularGrid(
        pb.Vector3d(min_b[0], min_b[1], min_b[2]),
        pb.Vector3d(max_b[0], max_b[1], max_b[2]),
        pb.Vector3d(cell_number[0], cell_number[1], cell_number[2]),
        cut=False,
    )

    # print("HERE***********************************")
    # print([s.x for s in r.rs.segments])
    # print([s.y for s in r.rs.segments])
    # ss
    # comm.barrier()
    # print("survived setRectangularGrid", rank)

    picker = lambda x, y, z: soil_model.pick(
        [x, y, z]
    )  #  function that return the index of a given position in the soil grid (should work for any grid - needs testing)
    r.ms.setSoilGrid(picker)  # maps segments, maps root segements and soil grid indices to each other in both directions
    # comm.barrier()
    # print("survived setSoilGrid", rank)

    # if rank == 0:
    # init_conductivities_const(r)

    return r


def write_files(file_name, psi_x, psi_i, sink, times, trans, psi_s, vol_, surf_, krs_, depth_, conc=None, c_=None):
    """saves numpy arrays ass npy files"""

    np.save("results/psix_" + file_name, np.array(psi_x))  # xylem pressure head per segment [cm]
    np.save("results/psiinterface_" + file_name, np.array(psi_i))  # pressure head at interface per segment [cm]
    np.save("results/sink_" + file_name, -np.array(sink))  # sink per segment [cm3/day]
    np.save("results/transpiration_" + file_name, np.vstack((times, -np.array(trans))))  # time [day], transpiration [cm3/day]
    np.save("results/soil_" + file_name, np.array(psi_s))  # soil potential per cell [cm]

    np.save("results/vol_" + file_name, np.array(vol_))  # volume per subType [cm3]
    np.save("results/surf_" + file_name, np.array(surf_))  # surface per subType [cm2]
    np.save("results/krs_" + file_name, np.array(krs_))  # soil potential per cell [cm2/day]
    np.save("results/depth_" + file_name, np.array(depth_))  # root system depth [cm]

    if conc is not None:
        np.save("results/soilc_" + file_name, np.array(conc))  # soil potential per cell [cm]
    if c_ is not None:
        np.save("results/nitrate_" + file_name, np.array(c_))  # soil potential per cell [cm]


def simulate_const(s, r, trans, sim_time, dt):
    """
    classic model:
    potential at root soil interface equals mean matric potential of surrounding finite volume
    """
    wilting_point = -15000  # cm
    skip = 6  # for output and results, skip iteration
    rs_age = 0.0  # day

    start_time = timeit.default_timer()
    psi_x_, psi_s_, sink_, x_, y_, psi_s2_ = [], [], [], [], [], []  # for post processing
    sx = s.getSolutionHead()  # inital condition, solverbase.py
    ns = len(r.rs.segments)
    if rank == 0:
        map_ = r.rs.seg2cell  # because seg2cell is a map
        mapping = np.array([map_[j] for j in range(0, ns)], dtype=np.int64)  # convert to a list

    N = int(np.ceil(sim_time / dt))

    """ simulation loop """
    for i in range(0, N):

        t = i * dt  # current simulation time

        """ 1. xylem model """
        if rank == 0:  # Root part is not parallel
            rx = r.solve(rs_age, -trans * sinusoidal2(t, dt), 0.0, sx, cells=True, wilting_point=wilting_point)  # xylem_flux.py
            fluxes = r.soilFluxes(rs_age, rx, sx, False)  # class XylemFlux is defined in MappedOrganism.h, approx = False
        else:
            fluxes = None
        fluxes = comm.bcast(fluxes, root=0)  # Soil part runs parallel

        """ 2. soil model """
        s.setSource(fluxes.copy())  # richards.py
        s.solve(dt)
        sx = s.getSolutionHead()  # richards.py

        """ validity check """

        """ remember results ... """
        if rank == 0 and i % skip == 0:

            sx_ = sx[:, 0]
            psi_x_.append(rx.copy())  # cm (per root node)
            psi_s_.append(np.array([sx_[ci] for ci in mapping]))  # cm (per root segment)
            sink = np.zeros(sx_.shape)
            for k, v in fluxes.items():
                sink[k] += v
            sink_.append(sink)  # cm3/day (per soil cell)
            x_.append(t)  # day
            y_.append(np.sum(sink))  # cm3/day
            psi_s2_.append(sx_)  # cm (per soil cell)

            min_sx, min_rx, max_sx, max_rx = np.min(sx), np.min(rx), np.max(sx), np.max(rx)
            n = round(float(i) / float(N) * 100.0)
            print(
                "\n["
                + "".join(["*"]) * n
                + "".join([" "]) * (100 - n)
                + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm root at {:g}, {:g}".format(min_sx, max_sx, min_rx, max_rx, np.sum(sink), -trans * sinusoidal2(t, dt))
            )

    if rank == 0:
        print("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

    return psi_x_, psi_s_, sink_, x_, y_, psi_s2_
