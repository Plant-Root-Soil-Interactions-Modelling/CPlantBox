""" Example of the photosynthesis module, using real data from the Selhausen lysimeter setup """
import sys; sys.path.append("../.."); sys.path.append("../../src/")  # |\label{l43:imports}|
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np 
import plantbox as pb
import visualisation.vtk_plot as vp
from datetime import datetime 
from matplotlib.dates import DateFormatter, HourLocator
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Photosynthesis import PhotosynthesisPython  # |\label{l43:importsPhotosynthesis}|

def getWeatherData(sim_time):
    diffDt = abs(pd.to_timedelta(weatherData['time']) - pd.to_timedelta(sim_time % 1,unit='d'))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    return weatherData.iloc[line_data]  # get the weather data for the current time step

""" Parameters and variables """
plant_age = 14  # plant age [day]       # |\label{l43:Parameters}|
sim_time = 1. # [day]
dt = 10./60./24.
N = int(sim_time/dt)
depth = 60  # soil depth [cm]
Hs = -1000  #  top soil matric potential [cm]

""" Weather data """
path = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(path + 'Selhausen_weather_data.txt', delimiter = "\t")  # |\label{l43:Tereno}|

""" soil """
soilSpace = pb.SDF_PlantContainer(np.inf, np.inf, depth, True)  # to avoid root growing aboveground
picker = lambda x, y, z: max(int(np.floor(-z)), -1)  # abovegroud nodes get index -1
p_s = np.linspace(Hs, Hs - depth, depth)  # water potential per soil layer  # |\label{l43:SoilEnd}|

""" plant """
plant = pb.MappedPlant()  # |\label{l43:plant}|
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021"
plant.readParameters(path + name + ".xml")

plant.setGeometry(soilSpace)  # creates soil space to stop roots from growing out of the soil
plant.setSoilGrid(picker)

plant.initialize(False)
plant.simulate(plant_age, False)  # |\label{l43:plantEnd}|


""" plant hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l43:hydraulicparams}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")  # |\label{l6h:hydraulic_end}|
hm = PhotosynthesisPython(plant, params)  # |\label{l43:PhotosynthesisPython}|

path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename = path + "photosynthesis_parameters")  # |\label{l43:read}|
# hm.write_photosynthesis_parameters(filename=path+"photosynthesis_parametersNew")   # |\label{l43:write}|

results = {'transpiration':[], 'gco2':[], 'An':[], 'Vc':[], 'Vj':[] }
for i in range(N):  # |\label{l43:loop}|
    """ Weather variables """
    weatherData_i = getWeatherData(plant_age)   # |\label{l52:weather}|

    """ Plant growth """
    plant_age += dt
    plant.simulate(dt, False)  # |\label{l52:plant}|   

    """ Plant transpiration and photosynthesis """
    hm.pCO2 = weatherData_i['co2']
    es = hm.get_es(weatherData_i['Tair'])
    ea = es * weatherData_i['RH']

    hm.solve(sim_time = plant_age, rsx = p_s, cells = True,
             ea = ea, es = es, 
             PAR = weatherData['PAR'][i] * (24 * 3600) / 1e4, # [mol/m2/s] -> [mol/cm2/d]
             TairC = weatherData_i['Tair'],
             verbose = 0)  # |\label{l43:solve}|
             
    """ Post processing """             
    hx = hm.get_water_potential()  # |\label{l43:results}|
    results['transpiration'].append(np.sum(hm.get_transpiration()) / 18 * 1e3)  # [cm3/day] * [mol/cm3] * [mmol/mol]
    results['An'].append(np.sum(hm.get_net_assimilation()) * 1e3)
    results['Vc'].append(np.sum(hm.get_Vc()) * 1e3)
    results['Vj'].append(np.sum(hm.get_Vj()) * 1e3)  # |\label{l43:resultsEnd}|
    
    print('at', weatherData['time'][i], 'mean water potential (cm)', np.round(np.mean(hx)),
            '\n\tin [mmol d-1], net assimilation:', np.round(np.sum(hm.get_net_assimilation()) * 1e3, 2),
            'transpiration:', np.round(np.sum(hm.get_transpiration()) / 18 * 1e3, 2))


time = [datetime.strptime(tt, '%H:%M:%S') for tt in weatherData['time']]
with plt.rc_context({
    'axes.labelsize': 15,     # axis labels
    'xtick.labelsize': 15,    # x-tick labels
    'ytick.labelsize': 15,    # y-tick labels
    'legend.fontsize': 14     # legend text
    }):
    fig, axs = plt.subplots(2, 2)  # |\label{l43:plot}|
    axs[0, 1].plot(time, weatherData['PAR'] * 1e3 * (24 * 3600) / 1e4, 'k', label =  'PAR [mmol cm-2 d-1]')
    axs[0, 1].plot(time, weatherData['Tair'] / 6.2, 'tab:red', label = 'T [°C]')
    axs[0, 1].set(ylabel = 'PAR\n[mmol cm-2 d-1]')
    secax_y = axs[0, 1].secondary_yaxis('right',
                                 functions=(lambda v: v*6.2, lambda v: v/6.2))
    secax_y.set_ylabel("T [°C]", color="tab:red")
    secax_y.tick_params(axis="y", colors="tab:red")
    secax_y.spines["right"].set_color("tab:red")

    axs[0, 0].plot(time, weatherData['RH'], 'tab:blue')
    axs[0, 0].set(ylabel = 'Relative humidity\n[-]')

    axs[1, 1].plot(time, results['An'], 'g', lw = 3 ,label = 'Net', zorder=1)
    axs[1, 1].plot(time, results['Vj'], 'k', label = 'Electron transport-limited')
    axs[1, 1].plot(time, results['Vc'], 'tab:red', label = 'Carboxilation-limited')
    axs[1, 1].set(xlabel = 'Time', ylabel = 'Total assimilation\n[mmol CO2 d-1]')
    axs[1, 1].legend(loc="upper left",frameon= False)

    axs[1, 0].plot(time, results['transpiration'], 'tab:blue')
    axs[1, 0].set(xlabel = 'Time', ylabel = 'Total transpiration\n[mmol H2O d-1]')
    axs[1, 0].xaxis.set_major_locator(MaxNLocator(5))

    for ax in axs.flatten():
        ax.xaxis.set_major_locator(HourLocator(range(0, 25, 6)))
        ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
        ax.fmt_xdata = DateFormatter('%H:%M')
        fig.autofmt_xdate()
        
    plt.show()
