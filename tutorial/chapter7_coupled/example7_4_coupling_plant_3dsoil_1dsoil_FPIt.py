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
from richards_no_mpi import RichardsNoMPIWrapper  # Python part of cylindrcial
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.Photosynthesis import PhotosynthesisPython # |\label{l41:imports_end}|
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeit



""" Parameters """  # |\label{l71c:param}|
sim_init = '00:00:00'
sim_end = '23:50:00'
path = "../../modelparameter/functional/climate/"
weatherData = pd.read_csv(path + 'Selhausen_weather_data.txt', delimiter = "\t")  # |\label{6h:Tereno}|
line_init = weatherData.index[weatherData['time'] == sim_init].tolist()[0]
line_end = weatherData.index[weatherData['time'] == sim_end].tolist()[0]  # |\label{6h:ParametersEnd}|
N = line_end - line_init
path = "../../modelparameter/structural/plant/"
name = "Triticum_aestivum_test_2021"  #"Anagallis_femina_Leitner_2010"  # Zea_mays_1_Leitner_2010, Anagallis_femina_Leitner_2010 TODO <<<<-------
wilting_point = -10000  # cm
plant_age = 14  # root system initial age [day]

""" soil """
min_b = [-4., -4., -25.]
max_b = [4., 4., 0.]
cell_number = [8, 8, 25]  # [1] spatial resolution
loam = [0.08, 0.43, 0.04, 1.6, 50]
initial = -600  # cm
nitrate_z = [0., -15., -15., -25.]  # top soil layer of 30 cm
nitrate_initial_values = np.array([5.e-3, 5.e-3, 1.e-3, 1.e-3]) / 0.43 / 1000  #  [kg/m3] -> [g/L]


""" Initialize macroscopic soil model """
s = RichardsWrapper(RichardsNCSP())  # |\label{l71c:soil}|
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
s.setHomogeneousIC(initial, True)  # [cm] total potential
s.setICZ_solute(nitrate_initial_values[::-1], nitrate_z[::-1])  # step-wise function, ascending order
s.setTopBC("noFlux")
s.setBotBC("freeDrainage")
s.setTopBC_solute("noFlux")
s.setBotBC_solute("outflow") # todo: add mass balance contorl including outflow of W and NO3
s.setParameter("Component.MolarMass", "6.2e-2")  # [kg/mol]
s.setParameter("Component.LiquidDiffusionCoefficient", "1.9e-9")  # [m2/s]
s.setVGParameters([loam])
# s.setParameter("Newton.EnableChop", "True")
# s.setParameter("Newton.EnableAbsoluteResidualCriterion", "True")
s.setParameter("Soil.SourceSlope", "1000")  # |\label{l71c:regularisation}|
s.initializeProblem()
s.setCriticalPressure(wilting_point)  # |\label{l71c:soil_end}|

""" Initialize xylem model """
plant = pb.MappedPlant()  # |\label{l71c:soil_plant}|
plant.readParameters(path + name + ".xml")
sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1)  # |\label{l71c:domain}|
plant.setGeometry(sdf)  # |\label{l71c:soil_plant_end}|

""" root hydraulic properties """
params = PlantHydraulicParameters()  # |\label{l71c:hydraulic}|
params.read_parameters("../../modelparameter/functional/plant_hydraulics/couvreur2012")
hm = PhotosynthesisPython(plant, params)
hm.wilting_point = wilting_point  # |\label{l71c:hydraulic_end}|
path = '../../modelparameter/functional/plant_photosynthesis/'
hm.read_photosynthesis_parameters(filename = path + "photosynthesis_parameters")  # |\label{6h:read}|
hm.pCO2 = weatherData['co2'][0]

""" Coupling (map indices) """
picker = lambda x, y, z: s.pick([x, y, z]) # |\label{l71c:coupling}|
plant.setSoilGrid(picker)
plant.initialize()
plant.simulate(plant_age, False)
hm.test()  # |\label{l71c:test}|

""" Numerical solution """
start_time = timeit.default_timer()
t = 0.
x_, y_ = [], []

for i in range(line_init+1, line_end):  # |\label{6h:loop}|

    dt = (pd.to_timedelta(weatherData['time'][i]) - pd.to_timedelta(weatherData['time'][max(0, i - 1)])).total_seconds() / (60 * 60 * 24)
    plant_age += dt
    plant.simulate(dt, False)  # |\label{l71c:plant}|
    hs = s.getSolutionHead()  # |\label{l71c:hs}|

    
    es = hm.get_es(weatherData['Tair'][i])
    ea = es * weatherData['RH'][i]

    hm.solve(sim_time = plant_age, rsx = hs, cells = True,
             ea = ea, es = es, PAR = weatherData['PAR'][i] * (24 * 3600) / 1e4, TairC = weatherData['Tair'][i])  # |\label{6h:solve}|

    hx = hm.get_water_potential()  # |\label{l71c:hx}|

    fluxes = hm.soil_fluxes(plant_age, hx, hs)  # |\label{l71c:soil_model}|
    s.setSource(fluxes)
    s.solve(dt)  # |\label{l71c:soil_model_end}|


    x_.append(t)
    y_.append(float(hm.get_transpiration()))  # |\label{l71c:results}|

    n = round(float(i) / float(N) * 100.)  # |\label{l71c:progress}|
    print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], [{:g}, {:g}] cm soil [{:g}, {:g}] cm plant at {:g} days {:g}"
            .format(np.min(hs), np.max(hs), np.min(hx), np.max(hx), s.simTime, hx[0]))

    t += dt  # [day]
    raise Exception

print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l71c:timing}|

""" VTK visualisation """  # |\label{l71c:plots}|
vp.plot_roots_and_soil(hm.ms, "pressure head", hx, s, True, np.array(min_b), np.array(max_b), cell_number, name)

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

