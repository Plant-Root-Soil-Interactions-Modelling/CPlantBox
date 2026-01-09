"""Plant including rhizosphere models"""

from datetime import datetime
import timeit

import matplotlib as mpl
from matplotlib.dates import DateFormatter, HourLocator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import plantbox as pb
from plantbox.functional.Photosynthesis import PhotosynthesisPython
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
import plantbox.functional.van_genuchten as vg  # van Genuchten model for soil hydraulic properties
from plantbox.visualisation import figure_style
import plantbox.visualisation.vtk_plot as vp
from rosi import helpful
from rosi.rhizo_models_plant import RhizoMappedSegments
from rosi.richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
from rosi.rosi_richardsnc import RichardsNCSP  # C++ part (Dumux binding)
from rosi.rosi_richardsnc_cyl import RichardsNCCylFoam  # C++ part (Dumux binding), macroscopic soil model


def getWeatherData(t):
    """retrieves the weather data for the current sim_time"""
    diffDt = abs(pd.to_timedelta(weather_data["time"]) - pd.to_timedelta(t % 1, unit="d"))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    return weather_data.iloc[line_data]


# Main parameters
path = "../../modelparameter/structural/plant/"
filename = "Triticum_aestivum_test_2021"
plant_age = 14.3  # root system h_s_initial age (day)
sim_time = 14.8
dt = 20 / 60 / 24  # days
n_steps = int((sim_time - plant_age) / dt)

# Weather data
path_weather = "../../modelparameter/functional/climate/"
weather_data = pd.read_csv(path_weather + "Selhausen_weather_data.txt", delimiter="\t")

# Bulk soil
box_min = [-4.0, -4.0, -24.0]
box_max = [4.0, 4.0, 0.0]
cell_number = [4, 4, 12]  # spatial resolution
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96]
vg_loam = vg.Parameters(hydrus_loam)
h_s_initial = -600  # cm
nitrate_initial_values = np.array([5.0e-3]) / 0.43 / 1000  #  kg m-3 -> g L-1


# Initialize macroscopic soil model
def setSoilParams(s):
    """sets DuMux soil parameters"""
    s.results_dir = "./results/"
    s.setParameter("Component.MolarMass", "6.2e-2")  # kg mol-1
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # m2 s-1
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Flux.UpwindWeight", "1")
    s.setVGParameters([hydrus_loam])
    s.wilting_point = -10000  # cm
    s.eps_regularization = 1e-10


s = RichardsWrapper(RichardsNCSP())
s.initialize()
s.createGrid(box_min, box_max, cell_number, periodic=True)  # cm
s.setHomogeneousIC(h_s_initial, True)  # total potential (cm)
s.setICZ_solute(nitrate_initial_values)  # step-wise function, ascending order
s.setTopBC("noFlux")
s.setBotBC("freeDrainage")
s.setTopBC_solute("noFlux")
s.setBotBC_solute("outflow")  # todo: add mass balance control including outflow of W and NO3
helpful.setDefault(s)  # other input parameters for the solver
setSoilParams(s)
s.initializeProblem()
s.setCriticalPressure(s.wilting_point)
s.setRegularisation(s.eps_regularization, s.eps_regularization)  # needs to be low when using sand parameters

# Initialize plant model
plant = pb.MappedPlant(1)
plant.readParameters(path + filename + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, box_max[2] - box_min[2] - 2.0)
plant.setGeometry(sdf)

# plant hydraulic properties
params = PlantHydraulicParameters()
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")
hm = PhotosynthesisPython(plant, params)
hm.wilting_point = s.wilting_point
path = "../../modelparameter/functional/plant_photosynthesis/"
hm.read_photosynthesis_parameters(filename=path + "photosynthesis_parameters2025")


# Coupling (map indices)
def picker(x, y, z):
    """soil grid cell index for positon (x, y, z)"""
    return s.pick([x, y, z])


plant.setSoilGrid(picker, noChanges=True)
plant.initialize()
plant.simulate(plant_age, False)

# Perirhizal zone models |\label{l74:perirhizal_models_start}|
rs = RhizoMappedSegments(soilModel=s, ms=plant, hm=hm, RichardsNCCylFoam=RichardsNCCylFoam)


def setSoilParamsCyl(s):
    """sets soil parameters and Michaelis Menten parameters"""
    setSoilParams(s)
    rs_uptake_vmax = 2.7e-6  # g cm-2 day-1, Roose and Kirk (2009)
    rs_uptake_km = 3.1e-6  # g cm-3, Roose and Kirk (2009)
    s.setInnerBC_solute(8)  # Michaelis Menten uptake
    s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(rs_uptake_vmax))  # active uptake parameters
    s.setParameter("RootSystem.Uptake.Km", s.dumux_str(rs_uptake_km))


rs.setSoilParam = setSoilParamsCyl  # |\label{l74:perirhizal_models_end}|

# Simulation loop
start_time = timeit.default_timer()
sim_times_, t_act_ = [], []
net_flux_water = np.zeros(np.prod(cell_number))
net_flux_solute = np.zeros(np.prod(cell_number))

for i in range(n_steps):  # |\label{l74:loop_start}|
    # Weather variables
    weatherData_i = getWeatherData(plant_age)

    # Plant growth |\label{l74:simulate_plant_start}|
    plant_age += dt
    plant.simulate(dt, False)
    rs.update()  # |\label{l74:simulate_plant_end}|

    # Plant transpiration  |\label{l74:plant_transpi_start}|
    h_rsi = rs.get_inner_heads(weatherData_i)  # inner values of the perirhizal models (cm)
    hm.pCO2 = weatherData_i["co2"]
    es = hm.get_es(weatherData_i["Tair"])
    ea = es * weatherData_i["RH"]

    hm.solve(
        sim_time=plant_age,
        rsx=h_rsi,
        cells=False,
        ea=ea,
        es=es,
        PAR=weatherData_i["PAR"] * (24 * 3600) / 1e4,  # mol photons m-2 s-1 -> mol photons cm-2 d-1
        TairC=weatherData_i["Tair"],
    )

    proposed_inner_fluxes_water = hm.radial_fluxes()  # cm3 day-1
    h_x = hm.get_water_potential()  # |\label{l74:plant_transpi_end}|

    # Perirhizal zone models  |\label{l74:perirhizal_start}|
    proposed_outer_fluxes_water = rs.splitSoilVals(soilVals=net_flux_water, compId=0, dt=dt)
    proposed_outer_fluxes_solute = rs.splitSoilVals(soilVals=net_flux_solute, compId=1, dt=dt)
    rs.solve(
        dt,
        proposed_inner_fluxes_water,  # inner BC water
        proposed_outer_fluxes_water,  # outer BC water
        proposed_outer_fluxes_solute,
    )  # outer BC solute 1 # |\label{l74:perirhizal_end}|

    # Bulk soil  |\label{l74:soil_model_start}|
    soil_source_water = rs.sumSegFluxes(rs.getRealisedInnerFluxes(0))  # cm3 day-1 per soil cell
    soil_source_solute = rs.sumSegFluxes(rs.getRealisedInnerFluxes(1))  # g day-1 per soil cell

    s.setSource(soil_source_water, 0)  # cm3 day-1, in richards.py
    s.setSource(soil_source_solute, 1)  # g day-1, in richards.py

    s.solve(dt, saveInnerFluxes_=True)  # |\label{l74:soil_model_end}|

    # Post processing
    rs.check1d3dDiff()  # |\label{l74:1d3d_diff_start}|

    # intracell exchange
    net_fluxes = -np.array([s.getFlowsPerCell(nc) for nc in range(s.numComp)])  # < 0 means leave the cell, > 0 means enter the cell
    net_flux_water = net_fluxes[0] - rs.alldiff1d3dCNW[0] / dt  # per soil cell (cm3 day-1)
    net_flux_solute = net_fluxes[1] - rs.alldiff1d3dCNW[1] / dt  # per soil cell (g day-1)  |\label{l74:1d3d_diff_end}|

    sim_times_.append(datetime.strptime(weatherData_i["time"], "%H:%M:%S"))
    t_act_.append(float(np.sum(hm.get_transpiration())))  # |\label{l74:transpi}|

    n = round(float(i) / float(n_steps - 1) * 100.0)
    h_soil = s.getSolutionHead()  # matric potential within the soil (cm)
    h_rsi_soil = np.delete(rs.get_inner_heads(weatherData_i), rs.airSegs)  # remove air segments
    print(f"[{'*' * n}{' ' * (100 - n)}], [{np.min(h_soil):g}, {np.max(h_soil):g}] cm bulk soil, [{np.min(h_rsi_soil):g}, {np.max(h_rsi_soil):g}] cm root-soil interface, [{np.min(h_x):g}, {np.max(h_x):g}] cm plant xylem at {weatherData_i['time']}")  # |\label{l74:info}|

print("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")

# VTK visualisation
vp.plot_plant_and_soil(hm.ms, "Xylem potential (cm)", h_x, s, False, np.array(box_min), np.array(box_max), cell_number, filename, sol_ind=1)

# Transpiration over time
fig, ax1 = figure_style.subplots12(1, 1)
ax1.plot(sim_times_, np.array(t_act_), "g")  # actual transpiration
ax2 = ax1.twinx()
ax2.plot(sim_times_, np.cumsum(np.array(t_act_) * dt), "c")  # cumulative transpiration
ax1.set_xlabel("Time (hh:mm)")
ax1.set_ylabel("Actual transpiration rate (mL day$^{-1}$)", color="g")
ax1.tick_params(axis="y", colors="g")
ax2.yaxis.label.set_color("c")
ax2.set_ylabel("Cumulative transpiration (mL)", color="c")
ax2.tick_params(axis="y", colors="c")
ax1.xaxis.set_major_locator(HourLocator(range(0, 25, 1)))
ax1.xaxis.set_major_formatter(DateFormatter("%H:%M"))
ax1.fmt_xdata = DateFormatter("%H:%M")
fig.autofmt_xdate()
plt.show()
