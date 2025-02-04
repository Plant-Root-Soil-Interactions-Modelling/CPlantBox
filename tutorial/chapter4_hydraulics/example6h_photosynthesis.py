""" Example of the photosynthesis module, using real data from the Selhausen lysimeter setup """
import sys; sys.path.append("../.."); sys.path.append("../../src/") # |\label{6h:imports}|
import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  
from functional.Photosynthesis import PhotosynthesisPython  # |\label{6h:importsPhotosynthesis}|

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import numpy as np # |\label{6h:importsEnd}|

""" Parameters and variables """
plant_age = 14  # plant age [day]       # |\label{6h:Parameters}|
kx = 4.32e-2  # axial conductivity [cm3/day]
kr = 1.728e-4  # radial conductivity [1/day]
gmax = 0.004  #  cm3/day radial conductivity of leaves = stomatal conductivity [1/day]
Hs = -1000  #  top soil matric potential [cm]
depth = 60 # soil depth [cm]

sim_init = '00:00:00'
sim_end = '23:50:00'
path = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(path + 'Selhausen_weather_data.txt', delimiter = "\t")      # |\label{6h:Tereno}|
line_init = weatherData.index[weatherData['time'] == sim_init].tolist()[0]
line_end = weatherData.index[weatherData['time'] == sim_end].tolist()[0]       # |\label{6h:ParametersEnd}|


""" soil """
min_ = np.array([-5, -5, -60])  # |\label{6h:Soil}|
max_ = np.array([9, 4, 0])
# res_ = np.array([5, 5, 5])
soilSpace = pb.SDF_PlantContainer(max_[0] - min_[0], max_[1] - min_[1], max_[2] - min_[2], True) #to avoid root growing aboveground
soil_index = lambda x,y,z : max(int(np.floor(-z)),-1) # abovegroud nodes get index -1
p_s = np.linspace(Hs, Hs - depth, depth) # water potential per soil layer  # |\label{6h:SoilEnd}|



""" plant """
plant = pb.MappedPlant() # |\label{6h:plant}|
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021" 
plant.readParameters(path + name + ".xml")

plant.setGeometry(soilSpace) # creates soil space to stop roots from growing out of the soil
plant.setSoilGrid(soil_index)

verbose = False
plant.initialize(verbose)
plant.simulate(plant_age,verbose) # |\label{6h:plantEnd}|

""" plant hydraulic properties """
params = PlantHydraulicParameters()  # |\label{6h:hydraulicparams}|
params.set_kr_const(kr,organType = pb.root)
params.set_kr_const(0.,organType = pb.stem)
params.set_kr_const(gmax,organType = pb.leaf)
params.set_kx_const(kx,organType = -1)
hm = PhotosynthesisPython(plant, params)  # |\label{6h:PhotosynthesisPython}|

path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename=path+"photosynthesis_parameters")  # |\label{6h:read}|
hm.pCO2 = weatherData['co2'][0]
# hm.write_photosynthesis_parameters(filename=path+"photosynthesis_parametersNew")   # |\label{6h:write}|

results = {'transpiration':[],'An':[],'Vc':[],'Vj':[] }
for i in range(line_init, line_end): # |\label{6h:loop}|
    
    dt = (pd.to_timedelta(weatherData['time'][i]) - pd.to_timedelta(weatherData['time'][max(0,i-1)])).total_seconds()/(60*60*24)
    plant_age += dt
    plant.simulate(dt,verbose)
    
    es = hm.get_es(weatherData['Tair'][i])
    ea = es * weatherData['RH'][i]
    
    
    hm.solve(sim_time = plant_age, rsx = p_s, cells = True,
             ea =  ea, es = es, PAR = weatherData['PAR'][i]* (24 * 3600) / 1e4, TairC=  weatherData['Tair'][i])  # |\label{6h:solve}|
    
    hx = hm.get_water_potential() # |\label{6h:results}|
    results['transpiration'].append(hm.get_tot_transpiration()/18 * 1e3) # [cm3/day] * [mol/cm3] * [mmol/mol]
    results['An'].append(np.sum(hm.get_net_assimilation()) * 1e3)
    results['Vc'].append(np.sum(hm.get_Vc()) * 1e3)
    results['Vj'].append(np.sum(hm.get_Vj()) * 1e3) # |\label{6h:resultsEnd}|
    
    print('at',weatherData['time'][i],'mean water potential (cm)',np.round(np.mean(hx)),
            '\n\tin [mmol d-1], net assimilation:',np.round(np.sum(hm.get_net_assimilation()) * 1e3,2),
            'transpiration:',np.round(hm.get_tot_transpiration()/18 * 1e3,2))
    

timePlot = weatherData[line_init:line_end]['time']

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(timePlot, results['An'])
axs[0, 0].set(xlabel = '', ylabel = 'total net actual\nassimilation rate[mmol CO2 d-1]')
axs[0, 0].xaxis.set_major_locator(MaxNLocator(5))
axs[1, 0].plot(timePlot, results['Vc'], 'tab:red')
axs[1, 0].set(xlabel = 'time', ylabel = 'total gross carboxilation-limited\nassimilation rate [mmol CO2 d-1]')
axs[1, 0].xaxis.set_major_locator(MaxNLocator(5))
axs[0, 1].plot(timePlot, results['Vj'], 'tab:brown')
axs[0, 1].set(xlabel = '', ylabel = 'total gross electron transport-limited\nassimilation rate [mmol CO2 d-1]')
axs[0, 1].xaxis.set_major_locator(MaxNLocator(5))
axs[1, 1].plot(timePlot, results['transpiration'], 'tab:brown')
axs[1, 1].set(xlabel = 'time', ylabel = 'total transpiration (mmol H2O d-1)') 
axs[1, 1].xaxis.set_major_locator(MaxNLocator(5))
plt.show()
