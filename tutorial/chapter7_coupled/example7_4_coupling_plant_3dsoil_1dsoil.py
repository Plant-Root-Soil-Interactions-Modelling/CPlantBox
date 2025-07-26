""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l41:imports}|
from rosi_richardsnc import RichardsNCSP   # C++ part (Dumux binding)
from rosi_richardsnc_cyl import RichardsNCCylFoam   # C++ part (Dumux binding), macroscopic soil model
from richards import RichardsWrapper  # Python part, macroscopic soil model
from rhizo_models_plant import RhizoMappedSegments
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Photosynthesis import PhotosynthesisPython # |\label{l41:imports_end}|
import functional.van_genuchten as vg  # van Genuchten model for soil hydraulic properties
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeit
from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()



""" Parameters """  # |\label{l71c:param}|
results_dir = './results/'
path = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(path + 'Selhausen_weather_data.txt', delimiter = "\t")  # |\label{6h:Tereno}|
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Anagallis_femina_Leitner_2010 TODO <<<<-------

plant_age = 14.  # root system initial age [day]
sim_end = 15.
dt = 1/60/24 # d
N = (sim_end - plant_age)/dt

""" soil """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 25]  # [1] spatial resolution
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -600  # cm
nitrate_z = [0., -15., -15., -25.]  # top soil layer of 30 cm
nitrate_initial_values = np.array([5.e-3, 5.e-3, 1.e-3, 1.e-3]) / 0.43 / 1000  #  [kg/m3] -> [g/L]


""" Initialize macroscopic soil model """
def setSoilParams(s):
    s.results_dir = results_dir  
    s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
    s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
    s.setVGParameters([loam])
    # s.setParameter("Newton.EnableChop", "True")
    # s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
    # s.setParameter("Soil.SourceSlope", "1000")  # |\label{l71c:regularisation}|
    s.wilting_point = -10000 # cm

s = RichardsWrapper(RichardsNCSP())  # |\label{l71c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1])  # step-wise function, ascending order
s.setTopBC("noFlux")
s.setBotBC("freeDrainage")
s.setTopBC_solute("noFlux")
s.setBotBC_solute("outflow") # todo: add mass balance contorl including outflow of W and NO3
setSoilParams(s)
s.initializeProblem()
s.setCriticalPressure(s.wilting_point)  # |\label{l71c:soil_end}|
cell_volumes = s.getCellVolumes()  # cm3
cell_volumes = comm.bcast(cell_volumes, root = 0)

""" Initialize xylem model """
plant = pb.MappedPlant(0)  # |\label{l71c:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  # |\label{l71c:domain}|
plant.setGeometry(sdf)  # |\label{l71c:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l71c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
hm = PhotosynthesisPython(plant, params)
hm.wilting_point = s.wilting_point  # |\label{l71c:hydraulic_end}|
path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename = path + "photosynthesis_parameters")  # |\label{6h:read}|
hm.pCO2 = weatherData['co2'][0]

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z]) # |\label{l71c:coupling}|
plant.setSoilGrid(picker, noChanges = True) 
plant.initialize()
plant.simulate(plant_age, False)
hm.test()  # |\label{l71c:test}|

""" rhizosphere models """
rs = RhizoMappedSegments(soilModel = s, ms = plant, hm = hm) 
split_type = 0  # type 0 == volume, type 1 == surface, type 2 == length

def setSoilParamsCyl(s):
    setSoilParams(s)
    s.setOuterBC("fluxCyl", 0.)  #  [cm/day] Neumann boundary condition
    s.setInnerBC("fluxCyl", 0.) 
    RS_Uptake_Vmax = 2.7e-6 # [g cm-2 day-1], Roose and Kirk (2009) 
    RS_Uptake_km =  3.1e-6 # [g cm-3], Roose and Kirk (2009)                 
    s.setInnerBC_solute(8) # Michaelis Menten uptake
    s.setParameter("RootSystem.Uptake.Vmax", s.dumux_str(RS_Uptake_Vmax))   # active uptake parameters
    s.setParameter("RootSystem.Uptake.Km", s.dumux_str(RS_Uptake_km))

rs.setSoilParam = setSoilParamsCyl

""" Numerical solution """
start_time = timeit.default_timer()
x_, y_ = [], []
net_flux_water = np.zeros(cell_volumes.shape)
net_flux_solute = np.zeros(cell_volumes.shape)

for i in range(N):  # |\label{6h:loop}|
    diffDt = abs(pd.to_timedelta(weatherData['time']) - pd.to_timedelta(plant_age % 1,unit='d'))
    line_data = np.where(diffDt == min(diffDt))[0][0]
    plant_age += dt
    plant.simulate(dt, False)  # |\label{l71c:plant}|
    rs.update()

    hx = rs.get_inner_heads()  # matric potential at the root soil interface, i.e. inner values of the cylindric models [cm]
    soil_k = np.divide(vg.hydraulic_conductivity(hx, loam), rs.radii)

    es = hm.get_es(weatherData['Tair'][line_data])
    ea = es * weatherData['RH'][line_data]

    if rank == 0:
        hm.solve(sim_time = plant_age, rsx = hx, cells = False,
                ea = ea, es = es, PAR = weatherData['PAR'][line_data] * (24 * 3600) / 1e4, 
                TairC = weatherData['Tair'][line_data])  # |\label{6h:solve}|
        proposed_inner_fluxes_water = np.array(plantModel.outputFlux)# [cm3/day] 
    else:
        proposed_inner_fluxes_water = None

    proposed_outer_fluxes_water = rs.splitSoilFluxes(net_flux_water / dt, split_type)
    proposed_outer_fluxes_solute = rs.splitSoilFluxes(net_flux_solute / dt, split_type)
    proposed_inner_fluxes_solute = np.zeros(proposed_outer_fluxes_solute)

    rs.solve(dt, proposed_inner_fluxes_water , # inner BC water
                                proposed_outer_fluxes_water, # outer BC water
                                [proposed_inner_fluxes_solute], # inner BC solute 1
                                [proposed_outer_fluxes_solute], # outer BC solute 1
                                )
    
    realised_inner_fluxes_water = rs.getXcyl(data2share=rs.seg_fluxes_limited, reOrder = True) 
    realised_inner_fluxes_solute = rs.getXcyl(data2share=rs.seg_fluxes_limited_CN_Out[0], reOrder = True) 
            
    soil_source_water = rs.sumSegFluxes(realised_inner_fluxes_water)  # [cm3/day]  per soil cell
    soil_source_solute = rs.sumSegFluxes(realised_inner_fluxes_solute)  # [g/day]  per soil cell

    s.setSource(soil_source_water.copy(), 0) # TODO UNITS ???? [cm3/day], in richards.py
    s.setSource(soil_source_solute.copy(), 1) # TODO UNITS ???? [cm3/day], in richards.py

    s.solve(dt)  # |\label{l71c:soil_model_end}|
    
    
    # intracell exchange
    net_fluxes = -np.array([s.getFluxesPerCell(nc) * dt for nc in range(s.numComp)]) # < 0 means leave the cell, > 0 means enter the cell
    net_flux_water = net_fluxes[0]
    net_flux_solute = net_fluxes[1:]
    #bulkSoil_sources = np.array([s.getSource(nc) * s.getCellVolumes() * dt for nc in range(s.numComp)]) # < 0 means net sink, > 0 means net source  
    #sources_wat_from3d =  bulkSoil_sources[0]# cm3
    #sources_CN_from3d =  bulkSoil_sources[1:] # mol

    x_.append(plant_age)
    y_.append(float(hm.get_transpiration()))  # |\label{l71c:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l71c:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm plant at {:g} days {:g}"
            .format(np.min(hs), np.max(hs), np.min(hx), np.max(hx), s.simTime, hx[0]))

    raise Exception

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l71c:timing}|

""" VTK visualisation """  # |\label{l71c:plots}|
vp.plot_roots_and_soil(hm.ms, "pressure head", hx, s,
                        True, np.array(min_b), np.array(max_b), cell_number, name)

""" Transpiration over time """
fig, ax1 = plt.subplots()
ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration (neumann)
ax2 = ax1.twinx()
ax2.plot(x_, np.cumsum(-np.array(y_) * dt), 'c--')  # cumulative transpiration (neumann)
ax1.set_xlabel("Time [d]")
ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
ax1.legend([ 'Actual', 'Cumulative'], loc = 'upper left')
np.savetxt(name, np.vstack((x_, -np.array(y_))), delimiter = ';')
plt.show()


