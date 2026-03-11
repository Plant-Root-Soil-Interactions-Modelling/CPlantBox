"""Example of the photosynthesis module, using real data from the Selhausen lysimeter setup"""

from datetime import datetime  # |\label{l43:imports}|
from matplotlib.dates import DateFormatter, HourLocator
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

import plantbox as pb
from plantbox.functional.Photosynthesis import PhotosynthesisPython  # |\label{l43:importsPhotosynthesis}|
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style


def getWeatherData(t):
    """get the weather data for time t"""
    diffDt = abs(pd.to_timedelta(weather_data["time"]) - pd.to_timedelta(t % 1, unit="d"))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    return weather_data.iloc[line_data]


# Parameters and variables
plant_age = 14.5  # plant age (day) |\label{l43:Parameters}|
sim_time = 10.0 / 60.0 / 24.0
depth = 60  # soil depth (cm)
Hs = -1000  # top soil matric potential (cm)

# Weather data
path = "../../modelparameter/functional/climate/"
weather_data = pd.read_csv(path + "Selhausen_weather_data.txt", delimiter="\t")  # |\label{l43:Tereno}|


# Soil
def picker(_x, _y, z):
    """soil grid cell index for positon (_x, _y, z)"""
    return max(int(np.floor(-z)), -1)  # aboveground nodes get index -1


soil_domain = pb.SDF_PlantContainer(np.inf, np.inf, depth, True)  # to avoid root growing aboveground
p_s = np.linspace(Hs, Hs - depth, depth)  # water potential per soil layer |\label{l43:SoilEnd}|

# Plant
plant = pb.MappedPlant(2)  # |\label{l43:plant}|
path = "../../modelparameter/structural/plant/"
filename = "Triticum_aestivum_test_2021"
plant.readParameters(path + filename + ".xml")

plant.setGeometry(soil_domain)  # creates soil space to stop roots from growing out of the soil
plant.setSoilGrid(picker)

plant.initialize(False)
plant.simulate(plant_age, False)  # |\label{l43:plantEnd}|

# Plant hydraulic properties
params = PlantHydraulicParameters()  # |\label{l43:hydraulicparams}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")  # |\label{l43:hydraulic_end}|
hm = PhotosynthesisPython(plant, params)  # |\label{l43:PhotosynthesisPython}|

path = "../../modelparameter/functional/plant_photosynthesis/"
hm.read_photosynthesis_parameters(filename=path + "photosynthesis_parameters")  # |\label{l43:read}|
# hm.write_photosynthesis_parameters(filename=path+"photosynthesis_parametersNew")   # |\label{l43:write}|

# Weather variables
weatherData_i = getWeatherData(plant_age)  # |\label{l43:weather}|

# Plant growth
plant_age += sim_time
plant.simulate(sim_time, False)  # |\label{l43:plant}|

# Plant transpiration and photosynthesis
hm.pCO2 = weatherData_i["co2"]
es = hm.get_es(weatherData_i["Tair"])
ea = es * weatherData_i["RH"]

hm.solve(
    sim_time=plant_age,
    rsx=p_s,
    cells=True,
    ea=ea,
    es=es,
    PAR=weatherData_i["PAR"] * (24 * 3600) / 1e4,  # (mol m-2 s-1) -> (mol cm-2 d-1)
    TairC=weatherData_i["Tair"],
    verbose=0,
)  # |\label{l43:solve}|

# Post processing
h_x = hm.get_water_potential()  # |\label{l43:results}|
fluxes = hm.radial_fluxes() # cm3/day
surfs = np.multiply(np.array(plant.segLength()), 2 * np.array(plant.radii) * np.pi)  # plant segment side surface [cm2]
fluxes = np.divide(fluxes, surfs)  # we convert to [cm3/(cm2 day)]

""" plot results """

# Additional vtk plot
ana = pb.SegmentAnalyser(hm.ms.mappedSegments())  # |\label{l43:sa}|
ana.addData("h_x", h_x)  # xylem potentials (cm)
ana.addData("radial_flux", fluxes)  
vp.plot_plant(ana, "radial_flux")  # |\label{l43:sa_end}|

# Output for Paraview
ana.write(
    "results/example4_3_planthydraulics.vtp",  # |\label{l43:paraview}|
    types=["radius", "subType", "age", "h_x", "radial_flux"],
)
