"""coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled"""

from pathlib import Path
import timeit

import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np

import plantbox as pb
from plantbox.functional.Perirhizal import PerirhizalPython  # |\label{l73:perirhizal}|
from plantbox.functional.PlantHydraulicModel import HydraulicModel_Doussan
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
import plantbox.visualisation.vtk_plot as vp
from rosi.richards import RichardsWrapper  # Python part
from rosi.rosi_richards import RichardsSP  # C++ part (Dumux binding)

from plantbox.visualisation import figure_style

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
max_rank = comm.Get_size()


def sinusoidal(t):
    """sinusoidal function (used for transpiration) (integral over one day is 1)"""
    return np.sin(2.0 * np.pi * np.array(t) - 0.5 * np.pi) + 1.0


# Parameters|\label{l73c:param}|
min_b = [-35.0, -10.0, -50.0]  # cm
max_b = [35.0, 10.0, 0.0]  # cm
cell_number = [17, 5, 50]  # ~4*4*1 cm3

path = "../../modelparameter/structural/rootsystem/"
filename = "Zeamays_synMRI_modified"  # "Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Zeamays_synMRI.xml
trans = 250  # cm3 day-1 (sinusoidal) = mL day-1
wilting_point = -15000  # cm
rs_age = 21  # root system initial age (day)

loam = [0.078, 0.43, 0.036, 1.56, 24.96]  # hydrus loam
initial = -400  # cm

sim_time = 7.5  # (day)
dt = 360.0 / (24 * 3600)  # days |\label{l73c:param_end}|

# Initialize macroscopic soil model
s = RichardsWrapper(RichardsSP())  # |\label{l73c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # cm
s.setHomogeneousIC(initial, True)  # total potential (cm)
s.setTopBC("noFlux")
s.setBotBC("noFlux")
s.setVGParameters([loam])
s.setParameter("Soil.SourceSlope", "1000")  # |\label{l73c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l73c:soil_end}|

# Initialize xylem model
plant = pb.MappedPlant(1)  # |\label{l73c:soil_plant}|
plant.readParameters(path + filename + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 2.)  # |\label{l73c:domain}|
plant.setGeometry(sdf)  # |\label{l73c:soil_plant_end}|

# root hydraulic properties
params = PlantHydraulicParameters()  # |\label{l73c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
# params.plot_conductivities(True)
hm = HydraulicModel_Doussan(plant, params)
hm.wilting_point = wilting_point  # |\label{l73c:hydraulic_end}|


# Coupling (map indices)
def picker(x, y, z):
    """soil grid cell index for positon (x, y, z)"""
    return s.pick([x, y, z])  # |\label{l73c:coupling}|


plant.setSoilGrid(picker)
plant.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False, False)  # |\label{l73c:rectgrid}|
plant.initialize(True)
plant.simulate(rs_age, True)
hm.test()  # |\label{l73c:test}|

# Perirhizal initialization
hs = s.getSolutionHead()
if rank == 0:
    peri = PerirhizalPython(hm.ms)  # |\label{l73c:peri}|
    # peri.set_soil(vg.Parameters(loam))  # |\label{l73c:perisoil}|
    home = Path.home()
    peri.open_lookup("results/hydrus_loam")  # |\label{l73c:peritable}|

    outer_r = peri.get_outer_radii("length")  # |\label{l73c:outer}|
    inner_r = peri.ms.radii
    rho_ = np.divide(outer_r, np.array(inner_r))  # |\label{l73c:rho}|

    # Numerical solution (a)
    start_time = timeit.default_timer()
    x_, y_, z_ = [], [], []
    hs_ = hm.ms.getHs(hs)  # matric potential per segment (cm) |\label{l73c:geths}|
    hsr = hs_.copy()  # initial values for fix point iteration # |\label{l73c:initial_hsr}|

source_water = None
hx = None
n_steps = round(sim_time / dt)
t = 0.0
for i in range(0, n_steps):  # |\label{l73c:loop}|
    if rank == 0:
        hx = hm.solve(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  # |\label{l73c:initial_hx}|
        hx_old = hx.copy()

        kr_ = hm.params.getKr(rs_age + t)  # |\label{l73c:update_kr}|
        inner_kr_ = np.multiply(inner_r, kr_)  # multiply for table look up; here const

        err = 1.0e6
        c = 0
        while err > 100.0 and c < 100:  # |\label{l73c:fixpoint}|

            # interpolation #
            hsr = peri.soil_root_interface_potentials(hx[1:], hs_, inner_kr_, rho_)  # |\label{l73c:interpolation}|

            # xylem matric potential
            hx = hm.solve_again(rs_age + t, -trans * sinusoidal(t), hsr, cells = False)  # |\label{l73c:hydraulic_hsr}|
            err = np.linalg.norm(hx - hx_old)
            hx_old = hx.copy()

            c += 1  # |\label{l73c:fixpoint_end}|

        fluxes = hm.radial_fluxes(rs_age + t, hx, hsr, cells = False)  # |\label{l73c:fluxes}|
        source_water = hm.sumSegFluxes(fluxes)

    water = s.getWaterVolume()  # |\label{l73c:domain_water}|
    source_water = comm.bcast(source_water, root = 0)
    s.setSource(source_water)  # TODO will be moved to MappedSegments # |\label{l73c:soil_fluxes}|
    s.solve(dt)
    soil_water = (s.getWaterVolume() - water) / dt  # |\label{l73c:domain_water_end}|

    hs = s.getSolutionHead()
    wc = s.getWaterContent()
    if rank == 0:
        x_.append(t)  # |\label{l73c:results}|
        y_.append(hm.get_transpiration(rs_age + t, hx.copy(), hsr.copy()))  # cm3 day-1
        z_.append(soil_water)  # cm3 day-1 |\label{l73c:results_end}|
        n = round(i / n_steps * 100)  # |\label{l73c:progress}|
        print(f"[{'*' * n}{' ' * (100 - n)}], {c:g} iterations, soil hs [{np.min(hs):g}, {np.max(hs):g}], interface [{np.min(hsr):g}, {np.max(hsr):g}] cm, root [{np.min(hx):g}, {np.max(hx):g}] cm, {s.simTime:g} days")

        if i % 10 == 0:  # |\label{l73c:write}|
            vp.write_soil_mpi(f"results/example73_{i // 10:06d}", hs, wc, min_b, max_b, cell_number)
            vp.write_plant(f"results/example73_{i // 10:06d}", hm.ms.plant())  # |\label{l73c:write_end}|

    t += dt

if rank == 0:
    print("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l73c:timing}|

# VTK visualisation |\label{l73c:plots}|
hx = comm.bcast(hx, root = 0)
vp.plot_roots_and_soil(hm.ms.mappedSegments(), "pressure head", hx, s, True, np.array(min_b), np.array(max_b), cell_number, filename)

if rank == 0:
    # Transpiration over time
    fig, ax1 = figure_style.subplots11()
    ax1.plot(x_, trans * sinusoidal(x_), "k")  # potential transpiration
    ax1.plot(x_, -np.array(y_), "g")  # actual transpiration (neumann)
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt), "c--")  # cumulative transpiration (neumann)
    ax1.set_xlabel("Time (day)")
    ax1.set_ylabel("Transpiration (mL day$^{-1}$) per plant")
    ax1.legend(["Potential", "Actual", "Cumulative"], loc = "upper left")
    np.save("results/" + filename + "_fp", np.vstack((x_, -np.array(y_), np.array(z_))))  # |\label{l72c:npsave}|
    plt.show()
