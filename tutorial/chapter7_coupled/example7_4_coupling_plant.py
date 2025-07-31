""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l74:imports}|
from rosi_richardsnc import RichardsNCSP   # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam   # C++ part (Dumux binding), macroscopic soil model
from richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
from rhizo_models_plant import RhizoMappedSegments
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Photosynthesis import PhotosynthesisPython # |\label{l74:imports_end}|
import functional.van_genuchten as vg  # van Genuchten model for soil hydraulic properties
import helpful  # helper functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeit
from datetime import datetime 
from matplotlib.dates import DateFormatter, HourLocator


""" Main parameters """  # |\label{l74:param}|
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021" 
plant_age = 14.3  # root system initial age [day]
sim_end = 14.8
dt = 20/60/24 # d
N = int((sim_end - plant_age)/dt)

""" Weather data """
pathWeather = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(pathWeather + 'Selhausen_weather_data.txt', delimiter = "\t")  # |\label{6h:Tereno}|

""" Bulk soil """
min_b = [-4., -4., -24.]
max_b = [4., 4., 0.]
cell_number = [4 , 4, 12]  # [1] spatial resolution
loam = [0.08, 0.43, 0.04, 1.6, 50]
vg_loam = vg.Parameters(loam)  
initial = -600  # cm
nitrate_initial_values = np.array([5.e-3]) / 0.43 / 1000 #  [kg/m3] -> [g/L]


""" Initialize macroscopic soil model """
def setSoilParams(s):
    s.results_dir = './results/'  
    s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Flux.UpwindWeight", "1")
    s.setVGParameters([loam])
    s.wilting_point = -10000 # cm

s = RichardsWrapper(RichardsNCSP())  # |\label{l74:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setICZ_solute(nitrate_initial_values)  # step-wise function, ascending order
s.setTopBC("noFlux")
s.setBotBC("freeDrainage")
s.setTopBC_solute("noFlux")
s.setBotBC_solute("outflow") # todo: add mass balance contorl including outflow of W and NO3
helpful.setDefault(s) # other input parameters for the solver
setSoilParams(s)
s.initializeProblem()
s.setCriticalPressure(s.wilting_point)  # |\label{l74:soil_end}|

""" Initialize plant model """
plant = pb.MappedPlant(1)  # |\label{l74:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  # |\label{l74:domain}|
plant.setGeometry(sdf)  # |\label{l74:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l74:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")  # |\label{l74:hydraulic_end}|
hm = PhotosynthesisPython(plant, params)
hm.wilting_point = s.wilting_point  # |\label{l74:hydraulic_end}|
path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename = path + "photosynthesis_parameters2025")  # |\label{6h:read}|

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z]) # |\label{l74:coupling}|
plant.setSoilGrid(picker, noChanges = True) 
plant.initialize()
plant.simulate(plant_age, False)

""" Perirhizal zone models """
rs = RhizoMappedSegments(soilModel = s, ms = plant, hm = hm, RichardsNCCylFoam = RichardsNCCylFoam) 

def setSoilParamsCyl(s):
    setSoilParams(s)
    RS_Uptake_Vmax = 2.7e-6 # [g cm-2 day-1], Roose and Kirk (2009) 
    RS_Uptake_km =  3.1e-6 # [g cm-3], Roose and Kirk (2009)                 
    s.setInnerBC_solute(8) # Michaelis Menten uptake
    s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(RS_Uptake_Vmax))   # active uptake parameters
    s.setParameter("RootSystem.Uptake.Km", s.dumux_str(RS_Uptake_km))

rs.setSoilParam = setSoilParamsCyl

""" Simulation loop """
start_time = timeit.default_timer()
x_, y_ = [], []
net_flux_water = np.zeros(np.prod(cell_number))
net_flux_solute = np.zeros(np.prod(cell_number))

for i in range(N):  # |\label{l74:loop}|
    """ Weather variables """
    diffDt = abs(pd.to_timedelta(weatherData['time']) - pd.to_timedelta(plant_age % 1,unit='d'))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    weatherData_ = weatherData.iloc[line_data]  # get the weather data for the current time step

    """ Plant growth """
    plant_age += dt
    plant.simulate(dt, False)  # |\label{l74:plant}|    
    rs.update()

    """ Plant transpiration """
    h_rsi = rs.get_inner_heads(weatherData_)  # inner values of the perirhizal models [cm]
    soil_k = np.divide(vg.hydraulic_conductivity(h_rsi, vg_loam), rs.radii)  # [cm3/day/cm] hydraulic conductivity at the root-soil interface, perirhizal model
    hm.pCO2 = weatherData_['co2']
    es = hm.get_es(weatherData_['Tair'])
    ea = es * weatherData_['RH']

    hm.solve(sim_time = plant_age, rsx = h_rsi, cells = False,
            ea = ea, es = es, 
            PAR = weatherData_['PAR'] * (24 * 3600) / 1e4, # [mol photons m-2 s-1] -> [mol photons cm-2 d-1] 
            TairC = weatherData_['Tair'])  # |\label{l74:solve}|
    
    proposed_inner_fluxes_water = hm.radial_fluxes() # [cm3/day] 
    h_xylem = hm.get_water_potential()
        
    """ Perirhizal zone models """
    proposed_outer_fluxes_water  = rs.splitSoilVals(soilVals = net_flux_water,  compId = 0, dt = dt)
    proposed_outer_fluxes_solute = rs.splitSoilVals(soilVals = net_flux_solute, compId = 1, dt = dt)         
    rs.solve(dt, proposed_inner_fluxes_water,  # inner BC water
                 proposed_outer_fluxes_water,  # outer BC water
                 proposed_outer_fluxes_solute) # outer BC solute 1

    """ Bulk soil """
    soil_source_water  = rs.sumSegFluxes(rs.getRealisedInnerFluxes(0)) # [cm3/day]  per soil cell # TODO: move sumSegFluxes to mapped segments
    soil_source_solute = rs.sumSegFluxes(rs.getRealisedInnerFluxes(1)) # [g/day]  per soil cell
    
    s.setSource(soil_source_water , 0) # [cm3/day], in richards.py
    s.setSource(soil_source_solute, 1) # [g/day], in richards.py

    s.solve(dt)  # |\label{l74:soil_model_end}|

    
    """ Post processing """
    rs.check1d3dDiff()
   
    # intracell exchange
    net_fluxes = - np.array([s.getFluxesPerCell(nc)  for nc in range(s.numComp)]) # < 0 means leave the cell, > 0 means enter the cell
    net_flux_water  = net_fluxes[0] - rs.alldiff1d3dCNW[0] /dt # [cm3/day] per soil cell
    net_flux_solute = net_fluxes[1] - rs.alldiff1d3dCNW[1] /dt # [g/day] per soil cell

    x_.append(datetime.strptime(weatherData_['time'], '%H:%M:%S'))
    y_.append(float(hm.get_transpiration()))  # |\label{l74:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l74:progress}|
    h_soil = s.getSolutionHead()  # pressure head in the soil [cm]
    h_rsi_soil = np.delete(rs.get_inner_heads(weatherData_),rs.airSegs)  # remove air segments
    
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm bulk soil, [{:g}, {:g}] cm root-soil interface, [{:g}, {:g}] cm plant xylem at {}"
        .format(np.min(h_soil), np.max(h_soil), np.min(h_rsi_soil), np.max(h_rsi_soil), np.min(h_xylem), np.max(h_xylem), weatherData_['time']))



print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l74:timing}|

""" VTK visualisation """  # |\label{l74:plots}|
vp.plot_plant_and_soil(hm.ms, "xylem pressure head (cm)", h_xylem, s,
                        False, np.array(min_b), np.array(max_b), cell_number, name,
                        sol_ind = 1)

""" Transpiration over time """
fig, ax1 = plt.subplots()
line1, = ax1.plot(x_, np.array(y_), 'g', label='Actual $[mL~d^{-1}]$')  # actual transpiration
ax2 = ax1.twinx()
line2, = ax2.plot(x_, np.cumsum(np.array(y_) * dt), 'c--', label='Cumulative $[mL]$')  # cumulative transpiration
ax1.set_xlabel("Time [hh:mm]")
ax1.set_ylabel("Transpiration per plant")
ax1.tick_params(axis='y', colors='g')
ax2.yaxis.label.set_color('c')
ax2.tick_params(axis='y', colors='c')
ax1.legend(handles=[line1, line2], loc = 'upper left')
ax1.xaxis.set_major_locator(HourLocator(range(0, 25, 1)))
ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
ax1.fmt_xdata = DateFormatter('%H:%M')
fig.autofmt_xdate()
plt.show()
