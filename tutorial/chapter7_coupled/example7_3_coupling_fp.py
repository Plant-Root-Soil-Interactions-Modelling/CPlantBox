"""coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled"""

import timeit

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
from plantbox.functional.Perirhizal import PerirhizalPython  # |\label{l73:perirhizal}|
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
from plantbox.visualisation import figure_style
import plantbox.visualisation.vtk_plot as vp
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding)


def sinusoidal(t):
    """sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0


# Parameters |\label{l73c:param}|
box_min = [-35.0, -10.0, -50.0]  # cm
box_max = [35.0, 10.0, 0.0]  # cm
cell_number = [17, 5, 50]  # ~4*4*1 cm3

path = "../../modelparameter/structural/rootsystem/"
filename = "Zeamays_synMRI_modified"
t_pot = 250  # cm3 /day = mL day-1
wilting_point = -15000  # cm
plant_age = 21  # root system h_s_initial age (day)

loam = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
h_s_initial = -400  # cm

sim_time = 7.5  # days
dt = 360.0 / (24 * 3600)  # days  |\label{l73c:param_end}|

# Initialize macroscopic soil model
s = RichardsWrapper(RichardsSP())  # |\label{l73c:soil}|
s.initialize()
s.createGrid(box_min, box_max, cell_number, periodic=True)  # cm
s.setHomogeneousIC(h_s_initial, True)  # total potential (cm)
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.setParameter("Soil.SourceSlope", "1000")  # |\label{l73c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l73c:soil_end}|

# Initialize xylem model
plant = pb.MappedPlant()  # |\label{l73c:soil_plant}|
plant.readParameters(path + filename + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, box_max[2] - box_min[2])  # |\label{l73c:domain}|
plant.setGeometry(sdf)  # |\label{l73c:soil_plant_end}|

# root hydraulic properties
params = PlantHydraulicParameters()  # |\label{l73c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
# params.plot_conductivities(True)
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = wilting_point  # |\label{l73c:hydraulic_end}|


# Coupling (map indices) #
def picker(x, y, z):
    """soil grid cell index for positon (_x, _y, z)"""
    return s.pick([x, y, z])  # |\label{l73c:coupling}|


plant.setSoilGrid(picker)
plant.setRectangularGrid(pb.Vector3d(box_min), pb.Vector3d(box_max), pb.Vector3d(cell_number), False, False)  # |\label{l73c:rectgrid}|
plant.initialize(True)
plant.simulate(plant_age, True)
hm.test()  # |\label{l73c:test}|

# Perirhizal initialization
peri = PerirhizalPython(hm.ms)  # |\label{l73c:peri}|
# peri.set_soil(vg.Parameters(loam))  # |\label{l73c:perisoil}|
peri.open_lookup("results/hydrus_loam")  # |\label{l73c:peritable}|

outer_r = peri.get_outer_radii("length")  # |\label{l73c:outer}|
inner_r = peri.ms.radii
rho_ = np.divide(outer_r, np.array(inner_r))  # |\label{l73c:rho}|

# Numerical solution
start_time = timeit.default_timer()
sim_times_, t_act_, q_soil_ = [], [], []

h_s = s.getSolutionHead_()  # matric potential (cm) |\label{l73c:initial_hs}|
h_s_ = hm.ms.getHs(h_s)  # matric potential per segment (cm) |\label{l73c:geths}|
h_sr = h_s_.copy()  # h_s_initial values for fix point iteration # |\label{l73c:initial_hsr}|

n_steps = round(sim_time / dt)
t = 0.0

for i in range(0, n_steps):  # |\label{l73c:loop}|
    h_x = hm.solve(plant_age + t, -t_pot * sinusoidal(t), h_sr, cells=False)  # |\label{l73c:initial_hx}|
    h_x_old = h_x.copy()

    kr_ = hm.params.getKr(plant_age + t)  # |\label{l73c:update_kr}|
    inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

    err = 1.0e6
    c = 0
    while err > 100.0 and c < 100:  # |\label{l73c:fixpoint}|
        # interpolation
        h_sr = peri.soil_root_interface_potentials(h_x[1:], h_s_, inner_kr_, rho_)  # |\label{l73c:interpolation}|

        # xylem matric potential
        h_x = hm.solve_again(plant_age + t, -t_pot * sinusoidal(t), h_sr, cells=False)  # |\label{l73c:hydraulic_hsr}|
        err = np.linalg.norm(h_x - h_x_old)
        h_x_old = h_x.copy()

        c += 1  # |\label{l73c:fixpoint_end}|

    water = s.getWaterVolume()  # |\label{l73c:domain_water}|
    fluxes = hm.radial_fluxes(plant_age + t, h_x, h_sr, cells=False)  # |\label{l73c:fluxes}|
    s.setSource(hm.sumSegFluxes(fluxes))  # TODO will be moved to MappedSegments # |\label{l73c:soil_fluxes}|
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt  # |\label{l73c:domain_water_end}|

    h_s = s.getSolutionHead()  # per cell |\label{l73c:new_hs}|
    h_s_ = hm.ms.getHs(h_s)  # per segment |\label{l73c:new_hs2}|

    sim_times_.append(t)  # |\label{l73c:results}|
    t_act_.append(hm.get_transpiration(plant_age + t, h_x.copy(), h_sr.copy()))  # cm3/day
    q_soil_.append(soil_water)  # cm3/day |\label{l73c:results_end}|

    n = round(i / n_steps * 100)  # |\label{l73c:progress}|
    print(f"[{'*' * n}{' ' * (100 - n)}], {c:g} iterations, soil h_s [{np.min(h_s):g}, {np.max(h_s):g}], interface [{np.min(h_sr):g}, {np.max(h_sr):g}] cm, root [{np.min(h_x):g}, {np.max(h_x):g}] cm, {s.simTime:g} days")

    if i % 10 == 0:  # |\label{l73c:write}|
        vp.write_soil(f"results/example73_{i // 10:06d}", s, box_min, box_max, cell_number)
        vp.write_plant(f"results/example73_{i // 10:06d}", hm.ms.plant())  # |\label{l73c:write_end}|

    t += dt

print("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l73c:timing}|

# VTK visualisation  |\label{l73c:plots}|
vp.plot_roots_and_soil(hm.ms.mappedSegments(), "pressure head", h_x, s, True, np.array(box_min), np.array(box_max), cell_number, filename)

# Transpiration over time
fig, ax1 = figure_style.subplots12(1, 1)
ax1.plot(sim_times_, t_pot * sinusoidal(sim_times_), "k")  # potential transpiration
ax1.plot(sim_times_, -np.array(t_act_), "g")  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(sim_times_, np.cumsum(-np.array(t_act_) * dt), "c--")  # cumulative transpiration (neumann)
ax1.set_xlabel("Time (day)")
ax1.set_ylabel("Transpiration (mL day$^{-1}$) per plant")
ax1.legend(["Potential", "Actual", "Cumulative"], loc="upper left")
np.save("results/" + filename + "_fp", np.vstack((sim_times_, -np.array(t_act_), np.array(q_soil_))))  # |\label{l72c:npsave}|
plt.show()
