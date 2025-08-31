""" Plant including rhizosphere models """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  
from rosi_richardsnc import RichardsNCSP   # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam   # C++ part (Dumux binding), macroscopic soil model
from richards_flat import RichardsFlatWrapper as RichardsWrapper  # Python part, macroscopic soil model
from rhizo_models_plant import RhizoMappedSegments
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Photosynthesis import PhotosynthesisPython 
import functional.van_genuchten as vg  # van Genuchten model for soil hydraulic properties
import helpful  # helper functions
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeit
from datetime import datetime 
from matplotlib.dates import DateFormatter, HourLocator
import matplotlib as mpl
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()

SMALL_SIZE = 20
MEDIUM_SIZE = 20
BIGGER_SIZE = 20
plt.rc('font', size = SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize = SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize = MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize = SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize = SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize = BIGGER_SIZE)  # fontsize of the figure title
mpl.rcParams['mathtext.default'] = 'regular'

def getWeatherData(sim_time):
    diffDt = abs(pd.to_timedelta(weatherData['time']) - pd.to_timedelta(sim_time % 1,unit='d'))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    return weatherData.iloc[line_data]  # get the weather data for the current time step

""" Main parameters """  
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021" 
plant_age = 14.3  # root system initial age [day]
sim_end = 14.8
dt = 20/60/24 # d
N = int((sim_end - plant_age)/dt)

""" Weather data """
pathWeather = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(pathWeather + 'Selhausen_weather_data.txt', delimiter = "\t")  

""" Bulk soil """
min_b = [-4., -4., -24.]
max_b = [4., 4., 0.]
cell_number = [4, 4, 12]  # [1] spatial resolution
hydrus_loam = [0.078, 0.43, 0.036, 1.56, 24.96] 
vg_loam = vg.Parameters(hydrus_loam)  
initial = -600  # cm
nitrate_initial_values = np.array([5.e-3]) / 0.43 / 1000 #  [kg/m3] -> [g/L]


""" Initialize macroscopic soil model """
def setSoilParams(s):
    s.results_dir = './results/'  
    s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
    s.setParameter("Newton.EnableChop", "True")
    s.setParameter("Flux.UpwindWeight", "1")
    s.setVGParameters([hydrus_loam])
    s.wilting_point = -10000 # cm
    s.eps_regularization = 1e-10

s = RichardsWrapper(RichardsNCSP()) 
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
s.setCriticalPressure(s.wilting_point)  
s.setRegularisation(s.eps_regularization, s.eps_regularization) # needs to be low when using sand parameters. 

""" Initialize plant model """
plant = pb.MappedPlant(1) 
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  
plant.setGeometry(sdf) 

""" plant hydraulic properties """
params = PlantHydraulicParameters()  
params.read_parameters("../../modelparameter/functional/plant_hydraulics/wheat_Giraud2023adapted")
hm = PhotosynthesisPython(plant, params)
hm.wilting_point = s.wilting_point  
path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename = path + "photosynthesis_parameters2025") 

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z]) 
plant.setSoilGrid(picker, noChanges = True) 
plant.initialize()
plant.simulate(plant_age, False)

""" Perirhizal zone models """ # |\label{l74:perirhizal_models_start}| 
rs = RhizoMappedSegments(soilModel = s, ms = plant, hm = hm, RichardsNCCylFoam = RichardsNCCylFoam) 

def setSoilParamsCyl(s):
    setSoilParams(s)
    RS_Uptake_Vmax = 2.7e-6 # [g cm-2 day-1], Roose and Kirk (2009) 
    RS_Uptake_km =  3.1e-6 # [g cm-3], Roose and Kirk (2009)                 
    s.setInnerBC_solute(8) # Michaelis Menten uptake
    s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(RS_Uptake_Vmax))   # active uptake parameters
    s.setParameter("RootSystem.Uptake.Km", s.dumux_str(RS_Uptake_km))
            

rs.setSoilParam = setSoilParamsCyl # |\label{l74:perirhizal_models_end}| 

""" Simulation loop """
start_time = timeit.default_timer()
x_, y_ = [], []
net_flux_water = np.zeros(np.prod(cell_number))
net_flux_solute = np.zeros(np.prod(cell_number))
proposed_outer_fluxes_water =  None
proposed_outer_fluxes_solute = None
proposed_inner_fluxes_water = None 
h_xylem = None

for i in range(N):  # |\label{l74:loop_start}|
    """ Weather variables """
    weatherData_i = getWeatherData(plant_age) 

    """ Plant growth """  # |\label{l74:simulate_plant_start}|  
    plant_age += dt
    plant.simulate(dt, False)   
    rs.update()  # |\label{l74:simulate_plant_end}|  

    """ Plant transpiration """ # |\label{l74:plant_transpi_start}|   
    h_rsi = rs.get_inner_heads(weatherData_i)  # inner values of the perirhizal models [cm]  
    hm.pCO2 = weatherData_i['co2']
    es = hm.get_es(weatherData_i['Tair'])
    ea = es * weatherData_i['RH']
    
    if rank == 0:
        hm.solve(sim_time = plant_age, rsx = h_rsi, cells = False,
                ea = ea, es = es, 
                PAR = weatherData_i['PAR'] * (24 * 3600) / 1e4, # [mol photons m-2 s-1] -> [mol photons cm-2 d-1] 
                TairC = weatherData_i['Tair'])  
        
        proposed_inner_fluxes_water = hm.radial_fluxes() # [cm3/day] 
        h_xylem = hm.get_water_potential() # |\label{l74:plant_transpi_end}| 
        
            
    """ Perirhizal zone models """ # |\label{l74:perirhizal_start}| 
    proposed_outer_fluxes_water  = rs.splitSoilVals(soilVals = net_flux_water,  compId = 0, dt = dt) 
    proposed_outer_fluxes_solute = rs.splitSoilVals(soilVals = net_flux_solute, compId = 1, dt = dt)   
    
    rs.solve(dt, proposed_inner_fluxes_water,  # inner BC water
                 proposed_outer_fluxes_water,  # outer BC water
                 proposed_outer_fluxes_solute) # outer BC solute 1 # |\label{l74:perirhizal_end}| 

    """ Bulk soil """ # |\label{l74:soil_model_start}|
    realisedInnerFlows_water = rs.getRealisedInnerFluxes(0)
    realisedInnerFlows_solute = rs.getRealisedInnerFluxes(1)
    
    soil_source_water  = comm.bcast(rs.sumSeg(realisedInnerFlows_water), root = 0) # [cm3/day]  per soil cell
    soil_source_solute = comm.bcast(rs.sumSeg(realisedInnerFlows_solute), root = 0) # [g/day]  per soil cell

    s.setSource(soil_source_water.copy() , 0) # [cm3/day], in richards.py
    s.setSource(soil_source_solute.copy(), 1) # [g/day], in richards.py

    s.solve(dt, saveInnerFluxes_ = True)  # |\label{l74:soil_model_end}|
    
    """ Post processing """
    rs.check1d3dDiff() # |\label{l74:1d3d_diff_start}|
   
    # inter-cell exchange
    net_fluxes = -np.array([s.getFlowsPerCell(nc)  for nc in range(s.numComp)]) # < 0 means leave the cell, > 0 means enter the cell
    
    if rank == 0:
        net_flux_water  = net_fluxes[0] - rs.alldiff1d3dCNW[0] /dt # [cm3/day] per soil cell
        net_flux_solute = net_fluxes[1] - rs.alldiff1d3dCNW[1] /dt # [g/day] per soil cell # |\label{l74:1d3d_diff_end}|

    h_soil = s.getSolutionHead()  # pressure head in the soil [cm]
    h_rsi_soil = rs.get_inner_heads(weatherData_i)
    
    if rank == 0:
        h_rsi_soil = np.delete(h_rsi_soil,rs.airSegs)  # remove air segments
        x_.append(datetime.strptime(weatherData_i['time'], '%H:%M:%S'))
        y_.append(float(np.sum(hm.get_transpiration())))  # |\label{l74:transpi}|

        n = round(float(i) / float(N - 1) * 100.)  
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm bulk soil, [{:g}, {:g}] cm root-soil interface, [{:g}, {:g}] cm plant xylem at {}"
            .format(np.min(h_soil), np.max(h_soil), np.min(h_rsi_soil), np.max(h_rsi_soil), np.min(h_xylem), np.max(h_xylem), weatherData_i['time']))  # |\label{l74:info}|


if rank == 0:
    print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s") 

""" VTK visualisation """  
h_xylem = comm.bcast(h_xylem, root = 0) 
vp.plot_plant_and_soil(hm.ms, "xylem pressure head (cm)", h_xylem, s,
                        False, np.array(min_b), np.array(max_b), cell_number, name,
                        sol_ind = 1)

if rank == 0:
    """ Transpiration over time """
    fig, ax1 = plt.subplots()
    ax1.plot(x_, np.array(y_), 'g')  # actual transpiration
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(np.array(y_) * dt), 'c')  # cumulative transpiration
    ax1.set_xlabel("Time [hh:mm]")
    ax1.set_ylabel("Actual transpiration rate $[mL~d^{-1}]$", color='g')
    ax1.tick_params(axis='y', colors='g')
    ax2.yaxis.label.set_color('c')
    ax2.set_ylabel("Cumulative transpiration $[mL]$", color='c')
    ax2.tick_params(axis='y', colors='c')
    ax1.xaxis.set_major_locator(HourLocator(range(0, 25, 1)))
    ax1.xaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax1.fmt_xdata = DateFormatter('%H:%M')
    fig.autofmt_xdate()
    plt.show()
