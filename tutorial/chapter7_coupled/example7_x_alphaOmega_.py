""" Vanderborght et al. 2023 """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.Perirhizal import PerirhizalPython as Perirhizal
import functional.van_genuchten as vg

from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import figure_style
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def make_source(q, area):
    s = {}
    for i in range(0, len(q)):
        if not np.isnan(q[i]):
            s[i] = -q[i] * area

    return s


""" Parameters """  # |\label{l7xa:param}|
depth = -100
N = 100
min_b = [-10., -10., depth]  # [cm]
max_b = [10., 10., 0.]  # [cm]
cell_number = [1, 1, N]  # ~[4*4*1] cm3

kx = 10 * 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]

path = "../../modelparameter/structural/rootsystem/"
name = "Zeamays_synMRI_modified"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Zeamays_synMRI.xml  <<<<-------
trans = 1.5 * 400  # cm3 /day (sinusoidal) = mL/day
wilting_point = -15000  # cm
rs_age = 60  # root system initial age [day]

loam = [0.01, 0.43, 0.0083, 1.2539, 2.272]  # jan paper
sp = vg.Parameters(loam)  # needed for Perirhizal class
vg.create_mfp_lookup(sp, wilting_point = -16000, n = 1501)  # needed for Perirhizal class
initial = -500  # cm (-330 )

sim_time = 7  # [day]
dt = 3600. / (24 * 3600)  # [days]  # |\label{l7xa:param_end}|

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())  # |\label{l7xa:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, False)  # [cm] False = matrix, True, = total potential
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.setParameter("Soil.SourceSlope", "100")  # |\label{l7xa:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l7xa:soil_end}|

""" Initialize xylem model """
plant = pb.MappedPlant()  # |\label{l7xa:soil_plant}|
plant.enableExtraNode()
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.5)  # |\label{l7xa:domain}|
plant.setGeometry(sdf)  # |\label{l7xa:soil_plant_end}|
plant.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False, False)  # needed for Perirhizal class

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l7xa:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
# params.set_kr_const(kr)
# params.set_kx_const(kx)
# params.plot_conductivities(True) # |\label{l7xa:plot_conductivities}|
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = wilting_point  # |\label{l7xa:hydraulic_end}|

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([0, 0, z])  # |\label{l7xa:coupling}|
plant.setSoilGrid(picker)
plant.initialize(True)
plant.simulate(rs_age, True)
hm.test()  # |\label{l7xa:test}|

peri = Perirhizal(plant)
h_bs = s.getSolutionHead()
h_sr = np.ones(h_bs.shape) * wilting_point

# k_prhiz = peri.perirhizal_conductance_per_layer(h_bs, h_sr, sp)  # test 1
# print("k_prhiz", np.nanmin(k_prhiz), np.nanmax(k_prhiz))
# plt.plot(k_prhiz, np.linspace(-50, 0, 50))
# plt.show()

# suf_ = hm.get_suf(rs_age + sim_time)  # test 2
# suf = peri.aggregate(suf_)
# print("suf", np.min(suf), np.max(suf), np.sum(suf))
# plt.plot(suf, np.linspace(depth, 0, N))
# plt.show()

# print("\nK_srs")
# k_srs = hm.get_soil_rootsystem_conductance(sim_time, h_bs, h_sr, sp)  # test 2
# print("k_srs", np.nanmin(k_srs), np.nanmax(k_srs))
# plt.plot(k_srs, np.linspace(-50, 0, 50))
# plt.show()
# # dd

""" Numerical solution """
start_time = timeit.default_timer()
t = 0.
x_, y_ = [], []
N = round(sim_time / dt)
area = (plant.maxBound.x - plant.minBound.x) * (plant.maxBound.y - plant.minBound.y)  # [cm2]

for i in range(0, N):  # |\label{l7xa:loop}|

    h_bs = s.getSolutionHead()
    h_bs = np.array(plant.matric2total(h_bs))

    start_time_ao = timeit.default_timer()

    hm.update(rs_age + sim_time)

    # Alpha: root system averaged stress factor
    # krs, _ = hm.get_krs(rs_age + sim_time)  # [cm2/day] (could be precomputed for static case)
    krs = hm.krs
    krs = krs / area
    print("krs", krs)

    k_srs = hm.get_soil_rootsystem_conductance(rs_age + sim_time, h_bs, wilting_point, sp)
    h_bs_diff = h_bs - np.ones(h_bs.shape) * wilting_point
    alpha = np.multiply(k_srs, h_bs_diff) / (-krs * wilting_point)  # [1]
    # print(alpha)
    # print("alpha", np.nanmin(alpha), np.nanmax(alpha))

    # Omega: root system averaged stress factor
    # suf_ = hm.get_suf(rs_age + sim_time)
    suf_ = hm.suf
    suf = peri.aggregate(suf_[0,:])
    # print(suf)
    # print("suf", np.min(suf), np.max(suf), np.sum(suf))
    alphaSUF = np.multiply(alpha, suf)
    omega = np.nansum(alphaSUF)  # note that nan are treated as 0
    print("omega", omega)
    # print("alphaSUF", np.nanmin(alphaSUF), np.nanmax(alphaSUF))

    # Omega_c: critical stress factor
    tp = trans * sinusoidal(t) / area  # potential tranpiration [cm3 day-1] -> [cm day-1]
    # print("tp", tp)
    omega_c = tp / (-wilting_point * krs)
    print("max uptake", (-wilting_point * krs), tp)
    print("omega_c", omega_c)
    print("omega / omega_c", omega / omega_c)

    # Sink, stressed
    q_s = alphaSUF * tp / omega_c
    # print("q_s", np.nansum(q_s), np.nanmin(q_s), np.nanmax(q_s))

    # Sink, unstressed
    denumerator = np.multiply(h_bs_diff, np.nansum(np.divide(alphaSUF, h_bs_diff)))
    print("denumerator", np.nansum(denumerator), np.nanmin(denumerator), np.nanmax(denumerator))
    print("- term: ", np.nansum(np.divide(alphaSUF, denumerator) * (omega / omega_c - 1) * tp))
    q_us = alphaSUF * tp / omega_c - np.divide(alphaSUF, denumerator) * (omega / omega_c - 1) * tp
    # print("q_us", np.nansum(q_us), np.nanmin(q_us), np.nanmax(q_us))

    print("pot", tp, "q_us", np.nansum(q_us), "q_s", np.nansum(q_s))

    if omega < omega_c:
        print("stressed")
        q = q_s
    else:
        print("unstressed")
        q = q_us

    start_time_soil = timeit.default_timer()

    fluxes = make_source(q, area)
    s.setSource(fluxes)
    s.solve(dt)  # |\label{l7xa:soil_model_end}|

    final_time = timeit.default_timer()

    x_.append(t)
    y_.append(-np.nansum(q) * area)  # |\label{l7xa:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l7xa:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], potential {:g}, actual {:g}; [{:g}, {:g}] cm soil at {:g} days"
            .format(tp * area, np.nansum(q) * area, np.min(h_bs), np.max(h_bs), s.simTime))

    print("wall times:", (start_time_ao - start_time_soil) / (start_time_ao - final_time), (start_time_soil - final_time) / (start_time_ao - final_time))

    if i % 10 == 0:  # |\label{l7xa:write}|
        vp.write_soil("results/example72_{:06d}".format(i // 10), s, min_b, max_b, cell_number)
        vp.write_plant("results/example72_{:06d}".format(i // 10), hm.ms.plant())  # |\label{l7xa:write_end}|
        # vp.plot_roots_and_soil(hm.ms.mappedSegments(), "matric potential", hx, s, True, np.array(min_b), np.array(max_b), cell_number) # BETTER output

    t += dt  # [day]

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l7xa:timing}|

""" VTK visualisation """  # |\label{l7xa:plots}|
# vp.plot_roots_and_soil(hm.ms.mappedSegments(), "matric potential", hx, s, True, np.array(min_b), np.array(max_b), cell_number)

""" Transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiratio
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
plt.show()
