""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.Perirhizal import PerirhizalPython  # |\label{l73:perirhizal}|
import functional.van_genuchten as vg
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import figure_style
import timeit


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


""" Parameters """  # |\label{l73c:param}|
min_b = [-35., -10., -50.]  # [cm]
max_b = [35., 10., 0.]  # [cm]
cell_number = [17, 5, 50]  # ~[4*4*1] cm3

path = "../../modelparameter/structural/rootsystem/"
name = "Zeamays_synMRI_modified"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Zeamays_synMRI.xml  <<<<-------
trans = 250  # cm3 /day (sinusoidal) = mL/day
wilting_point = -15000  # cm
rs_age = 21  # root system initial age [day]

loam = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
initial = -400  # cm

sim_time = 7.5  # [day]
dt = 360. / (24 * 3600)  # [days]  # |\label{l73c:param_end}|

""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsSP())  # |\label{l73c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.setParameter("Soil.SourceSlope", "1000")  # |\label{l73c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l73c:soil_end}|

""" Initialize xylem model """
plant = pb.MappedPlant()  # |\label{l73c:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  # |\label{l73c:domain}|
plant.setGeometry(sdf)  # |\label{l73c:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l73c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
# params.plot_conductivities(True)
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = wilting_point  # |\label{l73c:hydraulic_end}|

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z])  # |\label{l73c:coupling}|
plant.setSoilGrid(picker)
plant.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False, False)  # |\label{l73c:rectgrid}|
plant.initialize(True)
plant.simulate(rs_age, True)
hm.test()  # |\label{l73c:test}|

""" Perirhizal initialization """
peri = PerirhizalPython(hm.ms)  # |\label{l73c:peri}|
# peri.set_soil(vg.Parameters(loam))  # |\label{l73c:perisoil}|
peri.open_lookup("results/hydrus_loam")  # |\label{l73c:peritable}|

outer_r = peri.get_outer_radii("length")  # |\label{l73c:outer}|
inner_r = peri.ms.radii
rho_ = np.divide(outer_r, np.array(inner_r))  # |\label{l73c:rho}|

""" Numerical solution (a) """
start_time = timeit.default_timer()
x_, y_, z_ = [], [], []

hs = s.getSolutionHead_()  # [cm] matric potential # |\label{l73c:initial_hs}|
hs_ = hm.ms.getHs(hs)  # [cm] matric potential per segment  |\label{l73c:geths}|
hsr = hs_.copy()  # initial values for fix point iteration # |\label{l73c:initial_hsr}|

N = round(sim_time / dt)
t = 0.

for i in range(0, N):  # |\label{l73c:loop}|

    hx = hm.solve(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  # |\label{l73c:initial_hx}|
    hx_old = hx.copy()

    kr_ = hm.params.getKr(rs_age + t)  # |\label{l73c:update_kr}|
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

    err = 1.e6
    c = 0
    while err > 100. and c < 100:  # |\label{l73c:fixpoint}|

        """ interpolation """
        hsr = peri.soil_root_interface_potentials(hx[1:], hs_, inner_kr_, rho_)  # |\label{l73c:interpolation}|

        """ xylem matric potential """
        hx = hm.solve_again(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  # |\label{l73c:hydraulic_hsr}|
        err = np.linalg.norm(hx - hx_old)
        hx_old = hx.copy()

        c += 1  # |\label{l73c:fixpoint_end}|

    water = s.getWaterVolume()  # |\label{l73c:domain_water}|
    fluxes = hm.radial_fluxes(rs_age + t, hx, hsr, cells = False)  # |\label{l73c:fluxes}|
    s.setSource(hm.sumSegFluxes(fluxes))  # TODO will be moved to MappedSegments # |\label{l73c:soil_fluxes}|
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt  # |\label{l73c:domain_water_end}|

    hs = s.getSolutionHead()  # per cell |\label{l73c:new_hs}|
    hs_ = hm.ms.getHs(hs)  # per segment |\label{l73c:new_hs2}|

    x_.append(t)  # |\label{l73c:results}|
    y_.append(hm.get_transpiration(rs_age + t, hx.copy(), hsr.copy()))  # cm3/day
    z_.append(soil_water)  # cm3/day |\label{l73c:results_end}|

    n = round(float(i) / float(N) * 100.)  # |\label{l73c:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g} iterations, soil hs [{:g}, {:g}], interface [{:g}, {:g}] cm, root [{:g}, {:g}] cm, {:g} days"
          .format(c, np.min(hs), np.max(hs), np.min(hsr), np.max(hsr), np.min(hx), np.max(hx), s.simTime))

    if i % 10 == 0:  # |\label{l73c:write}|
        vp.write_soil("results/example73_{:06d}".format(i // 10), s, min_b, max_b, cell_number)
        vp.write_plant("results/example73_{:06d}".format(i // 10), hm.ms.plant())  # |\label{l73c:write_end}|

    t += dt

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l73c:timing}|

""" VTK visualisation """  # |\label{l73c:plots}|
vp.plot_roots_and_soil(hm.ms.mappedSegments(), "pressure head", hx, s, True, np.array(min_b), np.array(max_b), cell_number, name)

""" Transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, trans * sinusoidal(x_), 'k')  # potential transpiration
ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
np.save("results/" + name + "_fp", np.vstack((x_, -np.array(y_), np.array(z_))))  # |\label{l72c:npsave}|
plt.show()

