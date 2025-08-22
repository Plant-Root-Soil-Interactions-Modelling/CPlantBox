import os
import sys
sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../modelparameter/functional")
import plantbox as pb
from functional.xylem_flux import XylemFluxPython  # Python hybrid solver
import numpy as np
import visualisation.vtk_plot as vp # for quick 3d vizualisations
import matplotlib.pyplot as plt
from functional.phloem_flux import PhloemFluxPython  
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from plant_photosynthesis.wheat_FcVB_Giraud2023adapted import *
from plant_hydraulics.wheat_Giraud2023adapted import *
from plant_sucrose.wheat_phloem_Giraud2023adapted import *
import numpy as np
import pandas as pd
from datetime import datetime 
from matplotlib.dates import DateFormatter, HourLocator

def getWeatherData(sim_time):
    diffDt = abs(pd.to_timedelta(weatherData['time']) - pd.to_timedelta(sim_time % 1,unit='d'))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    return weatherData.iloc[line_data]  # get the weather data for the current time step

""" Parameters and variables """
plant_age = 7.3 # [day] init simtime
sim_time = 0.5 # [day]
dt = 1./24.
N = int(sim_time/dt)
depth = 60
p_mean = -600 # mean soil water potential [cm]

""" Weather data """
path = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(path + 'Selhausen_weather_data.txt', delimiter = "\t")  

""" plant """
plant = pb.MappedPlant(seednum = 2) 
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021" #"Triticum_aestivum_adapted_2023"
plant.readParameters(path + name + ".xml")

sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )
plant.setGeometry(sdf) 
verbose = False
plant.initialize(verbose )
plant.simulate(plant_age, verbose)

""" soil """
p_bot = p_mean + depth/2
p_top = p_mean - depth/2
sx = np.linspace(p_top, p_bot, depth) #soil water potential per voxel
picker = lambda x,y,z : max(int(np.floor(-z)),-1) 
plant.setSoilGrid(picker)  

""" Plant functional properties """
params = PlantHydraulicParameters()  # |\label{l52:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")  # |\label{l52:hydraulic_end}|
hm = PhloemFluxPython(plant, params, psiXylInit = min(sx),ciInit = weatherData['co2'][0]*0.5) # |\label{l52:phloempy}|
hm.wilting_point = -10000
path = '../../modelparameter/functional/'
hm.read_photosynthesis_parameters(filename = path + "plant_photosynthesis/photosynthesis_parameters2025")  # |\label{l52:read}|
hm.read_phloem_parameters(filename =path + "plant_sucrose/phloem_parameters2025")  # |\label{l52:read2}| 
# list_data = hm.get_phloem_data_list() # option of data that can be obtained from the phloem model
# hm.write_phloem_parameters(filename= 'phloem_parameters')  # |\label{l52:read_end}| 

time = []
cumulAssimilation = 0.
cumulTranspiration = 0.
Q_Rm_is, Q_Gr_is, Q_Exud_is, Q_Water_is = [], [], [], []

""" Simulation loop """
for i in range(N):
    """ Weather variables """
    weatherData_i = getWeatherData(plant_age)   # |\label{l52:weather}|

    """ Plant growth """
    plant_age += dt
    plant.simulate(dt, False)  # |\label{l52:plant}|   
    
    """ Plant transpiration and photosynthesis """
    hm.pCO2 = weatherData_i['co2']
    es = hm.get_es(weatherData_i['Tair'])
    ea = es * weatherData_i['RH']

    hm.solve(sim_time = plant_age, rsx = sx, cells = True,  
             ea = ea, es = es, 
             PAR = weatherData_i['PAR'] * (24 * 3600) / 1e4, 
             TairC = weatherData_i['Tair'],
             verbose = 0)  # |\label{l52:solve}|    
    
    """ Plant inner carbon balance """
    hm.solve_phloem_flow(plant_age, dt,  weatherData_i['Tair'])  # |\label{l52:balance}|       
                 
    """ Post processing """                   
    cumulAssimilation  += np.sum(hm.get_net_assimilation())  * dt #  [mol CO2]
    cumulTranspiration += np.sum(hm.get_transpiration()) * dt # [cm3]
    
    C_ST = hm.get_phloem_data( data = "sieve tube concentration")
    # cumulative
    Q_Rm    = hm.get_phloem_data( data = "maintenance respiration", doSum = True) 
    Q_Exud  = hm.get_phloem_data( data = "exudation", doSum = True)
    Q_Gr    = hm.get_phloem_data( data = "growth", doSum = True)
    Q_out   = Q_Rm + Q_Exud + Q_Gr
    # last time step
    Q_Rm_i   = hm.get_phloem_data( data = "maintenance respiration", last = True, doSum = True) 
    Q_Exud_i = hm.get_phloem_data( data = "exudation", last = True, doSum = True)
    Q_Gr_i   = hm.get_phloem_data( data = "growth", last = True, doSum = True)
    Q_out_i  = Q_Rm_i + Q_Exud_i + Q_Gr_i
    
    n = round(float(i) / float(N - 1) * 100.)
    print("\n[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "]")
    print("\t\tat ", int(np.floor(plant_age)),"d", int((plant_age%1)*24),"h, PAR:",  round(weatherData_i['PAR']*1e6),"mumol m-2 s-1")
    print("cumulative: transpiration {:5.2e} [cm3]\tnet assimilation {:5.2e} [mol]".format(cumulTranspiration, cumulAssimilation))
    print("sucrose concentration in sieve tube (mol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), max(C_ST)))     
    print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(Q_Rm_i/Q_out_i*100, 
         Q_Gr_i/Q_out_i*100,Q_Exud_i/Q_out_i*100))
    print("total aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(Q_Rm/Q_out*100, 
         Q_Gr/Q_out*100, Q_Exud/Q_out*100))   
    
    time.append(datetime.strptime(weatherData_i['time'], '%H:%M:%S'))
    Q_Rm_is.append(Q_Rm_i/dt)
    Q_Exud_is.append(Q_Exud_i/dt)
    Q_Gr_is.append(Q_Gr_i/dt)
    Q_Water_is.append(np.sum(hm.get_transpiration()))
    
    
""" Plot results """
fig, axs = plt.subplots(2,2)
axs[0,0].plot(time, Q_Rm_is)
axs[0,0].set(xlabel="Time [hh:mm]", ylabel='total respiration rate (mol/day)')
axs[1,0].plot(time, Q_Gr_is, 'tab:red')
axs[1,0].set(xlabel="Time [hh:mm]", ylabel='total growth rate (mol/day)')
axs[0,1].plot(time, Q_Exud_is , 'tab:brown')
axs[0,1].set(xlabel="Time [hh:mm]", ylabel='total exudation\nrate (mol/day)')
axs[1,1].plot(time, Q_Water_is , 'tab:brown')
axs[1,1].set(xlabel="Time [hh:mm]", ylabel='total transpiration\nrate (cm3/day)')
for ax in axs.flatten():
    ax.xaxis.set_major_locator(HourLocator(range(0, 25, 1)))
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.fmt_xdata = DateFormatter('%H:%M')
    fig.autofmt_xdate()
fig.tight_layout()
plt.show()
