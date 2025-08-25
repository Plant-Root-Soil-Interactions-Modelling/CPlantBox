""" coupling with DuMux as solver for the soil part, dumux-rosi must be installed & compiled """
import sys; sys.path.append("../.."); sys.path.append("../../src/")
sys.path.append("../../../dumux-rosi/build-cmake/cpp/python_binding/")  # dumux python binding
sys.path.append("../../../dumux-rosi/python/modules/")  # python wrappers

import plantbox as pb
import visualisation.vtk_plot as vp
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.Perirhizal import PerirhizalPython 
import functional.van_genuchten as vg
from rosi_richards import RichardsSP  # C++ part (Dumux binding)
from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import timeit
from scipy import interpolate 
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors
from pyevtk.hl import gridToVTK 
import os


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.
    
def theta2H(vg,theta):#(cm3/cm3) to cmH2O
    thetar =vg[0]# 
    thetas = vg[1]#
    alpha = vg[2]
    n = vg[3]#
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)
    
def soil_vg_(name:str):
    """ 
    Van Genuchten parameter for soil from Hydrus1D, 
    4D look up tables are created with thes script create_sra_table_v2
    """
    soil = {}
    soil["hydrus_loam"] = [0.078, 0.43, 0.036, 1.56, 24.96]
    soil["hydrus_clay"] = [0.068, 0.38, 0.008, 1.09, 4.8]
    soil["hydrus_sand"] = [0.045, 0.43, 0.145, 2.68, 712.8]
    soil["hydrus_sandyloam"] = [0.065, 0.41, 0.075, 1.89, 106.1]
    table_name = "table_{:s}".format(name)  # name for 4D look up table ########################
    return soil[name], table_name
    
def maize_(res:float):
    """ parameters for maize simulation """
    min_b = [-75/2, -15/2, -150.] #[cm]
    max_b = [75/2, 15/2, 0.]  #[cm]
    cell_number = [int(75/res), int(15/res), int(150/res)]
    return min_b, max_b, cell_number
    
def wheat_(res:float):
    """ parameters for wehat simulation """
    min_b = [-15/2, -3/2, -150.] #[cm]
    max_b = [15/2, 3/2, 0.]  #[cm]
    cell_number = [int(15/res), int(3/res), int(150/res)] 
    return min_b, max_b, cell_number
    
def set_scenario(plant_, res, soil_, initial_, trans_, rs_age, inf_:bool, evap_:bool):
    """ 
    Sets up a Scenario     
    
    plant            plant name: 'maize' or ?
    res              grid side length in (cm)
    wc_ini           initial water content [-] 
    soil             name of the soil 4D look up table, see soil_vg_() for the names and corresponding VG parameters
    trans            transpiration in (cm/d)  
    rs_age           root system age when RWU computation starts 
    """
    assert plant_ == "maize" or plant_ == "wheat", "plant should be 'maize', or 'wheat'"
    assert res == "high" or res == "low", "resolution should be high or low"
    assert soil_ in ["hydrus_loam", "hydrus_clay", "hydrus_sand", "hydrus_sandyloam"], "soil should be 'hydrus_loam', 'hydrus_clay', 'hydrus_sand' or 'hydrus_sandyloam' "
    assert (inf_ == False and evap_ == True) or (inf_ == True and evap_ == False) or (inf_ == False and evap_ == False), 'Infiltration must be false if evaporation is true or the other way round'

    # Hidden parameters
    wilting_point = -15000  # cm
    target = 150 #cm, length of target domain for gpr max in x, y, z
    wc_root = 0.8 #mean root water content (-)
    inf = 0.1 #cm/d - excess water that does no infiltrate is treated as runoff, i.e. is not accounted for 
    evap = -0.1 #cm/d

    # Numeric parameters
    random_seed = 1  # random seed

    if res == 'high': 
        res_ = 1
    elif res == 'low': 
        res_ = 3

    soil, table_name = soil_vg_(soil_)
    if plant_ == "maize":
        min_b, max_b, cell_number = maize_(res_)
        param_name = "Zeamays_synMRI_modified.xml"
        rh_params = 'couvreur2012'
    elif plant_ == "wheat":
        min_b, max_b, cell_number = wheat_(res_)
        param_name = "wheat_Morandage.xml"
        rh_params = 'wheat_Giraud2023adapted'

        
    cellvol = ((max_b[0]-min_b[0])/cell_number[0])*((max_b[1]-min_b[1])/cell_number[1])*((max_b[1]-min_b[1])/cell_number[1])
    X = np.linspace(min_b[0], max_b[0], cell_number[0])
    Y = np.linspace(min_b[1], max_b[1], cell_number[1])
    Z = np.linspace(min_b[2], 0, cell_number[2])
    
    
    trans = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1]) * trans_  # cm3 / day
    if initial_<0: #SWP
        initial = initial_
    elif 0<initial_<1: #wc
        initial = theta2H(soil,initial_)  # cm
    print('Initial soil water potential = ', str(initial), ' cm')
    
    """ Initialize macroscopic soil model """
    s = RichardsWrapper(RichardsSP()) 
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic = True)  # [cm]
    s.setHomogeneousIC(initial, False)  # [cm] constant total potential
    if inf_ == True: 
        s.setTopBC("atmospheric", 0.5, [[0, 1e10], [inf, inf]]) #infiltration 
    elif evap_ == True: 
        s.setTopBC("atmospheric", 0.5, [[0, 1e10], [evap, evap]]) #evaporation
    else: 
        s.setTopBC("noFlux")
    s.setBotBC("freeDrainage") #before: no flux
    s.setVGParameters([soil])
    s.setParameter("Soil.SourceSlope", "1000") 
    s.initializeProblem()
    s.setCriticalPressure(wilting_point)  

    """ Initialize xylem model """
    plant = pb.MappedPlant()  
    plant.readParameters("../../modelparameter/structural/rootsystem/" + param_name)
    if plant_ == "maize":
        print("maize dx modified")
        params = plant.getOrganRandomParameter(pb.root)
        for i in range(0, len(params)): 
            params[i].dx = 1
    elif plant_ == "wheat":
        print("wheat dx modified")
        params = plant.getOrganRandomParameter(pb.root)
        for i in range(0, len(params)): 
            params[i].dx = 1
    
    sdf = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.1) 
    plant.setSeed(random_seed)
    plant.setGeometry(sdf)  

    """ root hydraulic properties """
    params = PlantHydraulicParameters()  
    params.read_parameters("../../modelparameter/functional/plant_hydraulics/"+rh_params)
    # params.plot_conductivities(True)
    hm = HydraulicModel_Doussan(plant, params)
    hm.wilting_point = wilting_point  

    """ Coupling (map indices) """
    picker = lambda x, y, z: s.pick([x, y, z])
    plant.setSoilGrid(picker)
    plant.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b), pb.Vector3d(cell_number), False, False)  
    plant.initialize(True)
    plant.simulate(rs_age, True)
    hm.test()

    """ Perirhizal initialization """
    peri = PerirhizalPython(hm.ms)  
    # peri.set_soil(vg.Parameters(loam)) 
    peri.open_lookup("lookup/"+soil_) 

    outer_r = peri.get_outer_radii("length")  
    inner_r = peri.ms.radii
    rho_ = np.divide(outer_r, np.array(inner_r))  
  
    
    return inner_r, rho_, wilting_point, soil, s, peri, hm, plant, target, cell_number, cellvol, X, Y, Z, wc_root, res_
    
def write_npz(name, t, wc, hs, frac, rootvol, wc_stitch, hs_stitch, frac_stitch, rootvol_stitch): 
    
    cell_number = [np.shape(wc)[0],np.shape(wc)[1],np.shape(wc)[2]]
    cell_number_stitch = [np.shape(wc_stitch)[0],np.shape(wc_stitch)[1],np.shape(wc_stitch)[2]]
    
    directory = 'npz_'+name
    if not os.path.exists('results/'+directory):
        os.makedirs('results/'+directory)
    
    np.savez('results/'+directory+'/'+name+'_day'+str(int(t)), wc, hs, frac, rootvol, cell_number)
    np.savez('results/'+directory+'/'+name+'_stitched_day'+str(int(t)), wc_stitch, hs_stitch, frac_stitch, rootvol_stitch, cell_number_stitch)
    
def write_vtr(name, t, target, X, Y, Z, wc_root, wc, hs, rootvol, wc_stitch, hs_stitch, rootvol_stitch, plant, res): 

    X_stitch = np.linspace(-target/2, target/2, int(target/res))
    Y_stitch = np.linspace(-target/2, target/2, int(target/res))
    Z_stitch = np.linspace(-target, 0, int(target/res))
    
    wc_root_soil = wc+wc_root*rootvol
    wc_root_soil_stitch = wc_stitch+wc_root*rootvol_stitch
    
    cell_number = [np.shape(wc)[0],np.shape(wc)[1],np.shape(wc)[2]]
    cell_number_stitch = [np.shape(wc_stitch)[0],np.shape(wc_stitch)[1],np.shape(wc_stitch)[2]]
    
    directory = 'vtr_'+name
    if not os.path.exists('results/'+directory):
        os.makedirs('results/'+directory)
    
    plant.write('results/'+directory+'/'+name+'.vtp')
    gridToVTK('results/'+directory+'/'+name+'_day'+str(int(t)), X, Y, Z, pointData = {"Water content":wc, "Soil water potential":hs, "Root volume":rootvol,"water content root and soil":wc_root_soil })
    gridToVTK('results/'+directory+'/'+name+'_stitch_day_'+str(int(t)), X_stitch, Y_stitch, Z_stitch, pointData = {"Water content":wc_stitch,"Soil water potential":hs_stitch,"Root volume":rootvol_stitch, "water content root and soil":wc_root_soil_stitch})
