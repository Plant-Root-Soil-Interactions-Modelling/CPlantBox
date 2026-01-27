"""
Run simulation
"""
import sys; 
from pathlib import Path

#edit CPlantBox dir
CPlantBox_dir =  "../.." # "/data2model_0807/CPlantBox"

from mpi4py import MPI; comm = MPI.COMM_WORLD; rank = comm.Get_rank(); max_rank = comm.Get_size()
sys.path.append(CPlantBox_dir)
sys.path.append( CPlantBox_dir+"/src")
sys.path.append( CPlantBox_dir+"/gui/viewer")


import plantbox as pb
import visualisation.vtk_plot as vp
import viewer_conductivities
from functional.PlantHydraulicParameters import PlantHydraulicParameters  # |\label{l42:imports}|
from functional.PlantHydraulicModel import HydraulicModel_Meunier_large, HydraulicModel_Meunier   # |\label{l42:imports_end}|


import numpy as np
from structural.Plant import PlantPython
from structural.MappedOrganism import MappedPlantPython
from functional.Perirhizal import PerirhizalPython as Perirhizal
import matplotlib.pyplot as plt
import copy
import os
import pickle
from scipy.stats import gaussian_kde
import time
import time
from collections import defaultdict
from functional.xylem_flux import sinusoidal2

def tic():
    return time.perf_counter()

def toc(t0):
    return time.perf_counter() - t0

##### add soil
sys.path.append(CPlantBox_dir+"/../dumux-rosi/python/modules");
sys.path.append(CPlantBox_dir+"/../dumux-rosi/build-cmake/cpp/python_binding/");
from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PlantHydraulicModel import HydraulicModel_Meunier
from functional.Perirhizal import PerirhizalPython as Perirhizal
import functional.van_genuchten as vg


from richards import RichardsWrapper  # Python part
import numpy as np
import matplotlib.pyplot as plt
import timeit

#def sinusoidal(t):
#    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.


def make_source(q_plant, area, ncell_soil):
    q = np.flip(np.pad(q_plant, (0, ncell_soil - len(q_plant))))
    s = {}
    for i in range(0, len(q)):
        #assert not np.isnan(q[i])
        if not np.isnan(q[i]):#: # if there are not roots in that layer
            if q[i] != 0:
                s[i] = -q[i] * area

    return s

""" Parameters """  
yr_to_BEDD = 1225 # BEDD/yr
day_to_BEDD = yr_to_BEDD/365. # BEDD/day
periodic = True
trans = 2000 # mL/day (sinusoidal)
wilting_point = -15000  # cm
loam = [0.1406, 0.4148, 0.04052, 1.32416, 38.43, -2.067]  
sp = vg.Parameters(loam)  # needed for Perirhizal class

simtime_soil = 7 * 4 * 6  # from April to October


length_x = 100
length_y = 200
area = length_x * length_y  # [cm2], TODO: check
soilSpace = pb.SDF_PlantContainer(1e6, 1e6,1e6, True)  # to avoid root growing aboveground

dz = 1.
volumes = dz * area
picker = lambda x, y, z: -int(z) 
####

def do_soil_simulation(plant, hm, hm_meunier, depth, suf, #plant_ordered
                       name_str, scenario, simSoil):
    is_E0_50 = False#(scenario == 0) and (simSoil == 50)
    if is_E0_50:
        from rosi_richards import RichardsSPnum as RichardsSP # C++ part (Dumux binding) RichardsSPnum as
    else:
        from rosi_richards import RichardsSP # C++ part (Dumux binding) RichardsSPnum as 
    depth_soil = depth #750 #
    print('depth_soil',depth_soil)
    assert depth_soil >= depth
    min_b = [-length_x/2., -length_y/2., -depth_soil]  # [cm]
    max_b = [length_x/2., length_y/2., 0.]  # [cm]
    cell_number = [1, 1, abs(int(depth_soil))]  # [cm3]

    """ Initialize macroscopic soil model """
    s = RichardsWrapper(RichardsSP())  
    s.initialize()
    s.createGrid(min_b, max_b, cell_number, periodic)  # [cm]
    
    if scenario == 0:
        initial =  -330  
    elif scenario == 1:
        initial = -5000         # cm # cf https://www.sciencedirect.com/science/article/pii/S0304380015002343#fig0010
    else:
        raise Exception
        
    s.setHomogeneousIC(initial, True) 
    s.setTopBC("noFlux")
    s.setVGParameters([loam])
    
    if False:
        s.setParameter("Soil.SourceSlope", "1000") 
        s.setBotBC("noFlux") #"freeDrainage")
        dt_soil = 3600. / (24 * 3600) * 10.
    else:
        if is_E0_50:
            s.MaxRelativeShift = 1e-6
            s.setParameter("Newton.MaxRelativeShift", str(s.MaxRelativeShift))
        s.MaxTimeStepDivisions = 20
        s.setParameter("Newton.MaxTimeStepDivisions",
                         str( s.MaxTimeStepDivisions) )  
        s.MaxSteps = 50
        s.setParameter("Newton.MaxSteps",
                         str( s.MaxSteps) ) 
        s.setParameter("Soil.SourceSlope", "1000") 
        s.setBotBC("freeDrainage")
        dt_soil = 3600. / (24 * 3600) * 1.
        dt_soil_night = dt_soil * 12.
        dt_soil_day = dt_soil
        s.ddt = 1.e-3
    
    N_soil = round(simtime_soil / dt_soil)
    s.maxDt = 3600. / (3600.*24.)
    s.initializeProblem(s.maxDt)
    
    s.setCriticalPressure(wilting_point)  
    #s.setRegularisation(1e-3,1e-3)
    depth_array = np.flip([cc[2] for cc in s.getCellCenters()])[:len(suf)] #plant_ordered
    #print(depth_array)
    # assert abs(float(s.getParameter(s.param_group + "VanGenuchten.L")) - loam[-1]) < 1e-5 # ok
    
    '''
    fig, ax1 = plt.subplots()
    ax1.plot(suf, depth_array, label = "SUF [-]")
    ax1.set_xlabel("Fraction [-]")
    ax1.set_ylabel("Depth [cm]")
    plt.show()
    '''
    
    """ Plant """
    # Alpha: root system averaged stress factor
    # krs, _ = hm.get_krs(rs_age + sim_time)  # [cm2/day] (could be precomputed for static case)
    krs = hm.krs #get_krs_forsoil()
    krs = krs / area
    
    """ Numerical solution """
    start_time = timeit.default_timer()
    top_ind = s.pick([0., 0., -0.5])
    dt_rain = 0.
    t_day = 0.
    x_, y_ = [], []
    h_bss, qs = [], []
    top_new = []
    torain = 0.
    verbose = False
    first_error = True
    for i in range(0, N_soil): 
        start_time_ao = timeit.default_timer()

        h_bs = np.flip(s.getSolutionHead())[:len(suf)]
        if verbose:
            print('suf',suf)
            print('h_bs',h_bs)
        h_bs += depth_array #np.array(plant.matric2total(h_bs))
        if verbose:
            print(h_bs)
        h_sr = np.ones(h_bs.shape) * wilting_point
        if verbose:
            print('h_sr',h_sr,'wilting_point',wilting_point)
    
        k_srs = hm.get_soil_rootsystem_conductance(plant.getSimTime(), area, suf, h_bs, wilting_point, sp)
        if verbose:
            print('k_srs',k_srs, hm.krs)#get_krs_forsoil())
        h_bs_diff = h_bs - np.ones(h_bs.shape) * wilting_point
        if verbose:
            print('h_bs_diff',h_bs_diff)
        alpha = np.multiply(k_srs, h_bs_diff) / (-krs * wilting_point)  # [1]
        if verbose:
            print('alpha',alpha)

        # Omega: root system averaged stress factor
        alphaSUF = np.multiply(alpha, suf)
        if verbose:
            print('alphaSUF',alphaSUF)
        omega = np.nansum(alphaSUF)  # note that nan are treated as 0
        if verbose:
            print('omega',omega)
            
        # Omega_c: critical stress factor
        tp = trans * sinusoidal2(t_day, dt_soil) / area  # potential tranpiration [cm3 day-1] -> [cm day-1]
        #tp = trans * 1. / area  # potential tranpiration [cm3 day-1] -> [cm day-1]; print("att: remove dummy transpiration")
        
        if (abs(t_day% 1 * 24. - 18.) < 1e-6 ):
            print("do night", tp)
            dt_soil = dt_soil_night
        else:
            dt_soil = dt_soil_day
            
            
        if verbose:
            print("tp", tp, sinusoidal2(t_day, dt_soil))
        omega_c = tp / (-wilting_point * krs)
        if verbose:
            print('omega_c',omega_c)

        # Sink, stressed
        if(tp > 0):
            q_s = alphaSUF * tp / omega_c
        else:
            q_s = np.zeros(len(alphaSUF))
        #assert not np.isnan(q_s).any(), "q_us array contains NaN values!"
        if verbose:
            print('qs',qs)
        # Sink, unstressed
        denumerator = np.multiply(h_bs_diff, np.nansum(np.divide(alphaSUF, h_bs_diff)))
        if(tp > 0):
            q_us = alphaSUF * tp / omega_c - np.divide(alphaSUF, denumerator) * (omega / omega_c - 1) * tp
        else:
            q_us = np.zeros(len(alphaSUF))
        #assert not np.isnan(q_us).any(), "q_us array contains NaN values!"
        if verbose:
            print('q_us',q_us)

        if omega < omega_c:
            print("stressed", s.ddt,end = ", ")
            q = q_s
        else:
            print("unstressed", s.ddt,end = ", ")
            q = q_us

        start_time_soil = timeit.default_timer()

        fluxes = make_source(q, 
                             area, depth_soil)
        
        s.setSource(fluxes)
        
        if (scenario == 1) and (dt_rain > 7.) and (abs(t_day% 1 * 24. - 6.) < 1e-6 ): # 12.5 mm over 24h
            dt_rain = 0
            torain = 24.
            
        if torain > 0:
            s.setTopBC("constantFlux", 12.5/10.)  #  [cm/day] atmospheric is with surface run-off   |\label{l61i:top_bc}|
        else:
            s.setTopBC("noFlux")            
        torain -= (24.*dt_soil)
        s.solve(dt_soil)  # [day]
        
        topflux = s.getNeumann(top_ind)
        top_new.append(topflux)
        final_time = timeit.default_timer()

        x_.append(t_day)
        y_.append(-np.nansum(q) * area) 
        qs.append(q)
        h_bss.append(h_bs)

        n = round(float(i) / float(N_soil) * 100.)  
        print(topflux , torain, dt_rain)
        print("[" + ''.join(["*"]) * n + ''.join([" "]) * (100 - n) + "], {:g}, potential {:g}, actual {:g}; [{:g}, {:g}] cm soil, {:g} flux at {:g} days"
                .format(sinusoidal2(t_day, dt_soil), tp * area, np.nansum(q) * area, np.min(h_bs), np.max(h_bs), topflux, s.simTime))

        t_day  += dt_soil  # [day]
        dt_rain  += dt_soil

        #print ("Coupled benchmark solved in ", timeit.default_timer() - start_time, " s")  # |\label{l7xa:timing}|
        if verbose:
            print("qq",q)
    h_bs -= depth_array # total 2 matric
    rx = hm.solve_again(dt_soil,  -trans * sinusoidal2(t_day, dt_soil), hm.get_hs(h_bs), cells = False)
    #rx_meunier = hm_meunier.solve(dt_soil,  trans * sinusoidal2(t_day, dt_soil), hm.get_hs(h_bs), cells = False)
    '''
    """ Transpiration over time """
    fig, ax1 = plt.subplots()
    ax1.plot(x_, trans * sinusoidal2(x_, dt_soil), 'k')  # potential transpiration
    ax1.plot(x_, -np.array(y_), 'g')  # actual transpiration
    ax2 = ax1.twinx()
    ax2.plot(x_, np.cumsum(-np.array(y_) * dt_soil), 'c--')  # cumulative transpiratio
    ax1.set_xlabel("Time [d]")
    ax1.set_ylabel("Transpiration $[mL d^{-1}]$ per plant")
    ax1.legend(['Potential', 'Actual', 'Cumulative'], loc = 'upper left')
    plt.show()
    
    
    '''
    vp.plot_roots_and_soil(hm.ms.mappedSegments(), "pressure head",rx, s, periodic, min_b, max_b, 
                           cell_number, filename = name_str, interactiveImage = False, dorender = False)  # VTK vizualisation
    
    plant_soil = {
            'time':x_,
            'transpiration':y_,
            'h_bs':h_bs,
            'q':q,
            'depth_array':depth_array,
        'top_new':top_new
            }
    return plant_soil
    

def get3Dshape(plant,title_ = 'wine',data={}, saveOnly = True, show = "subType"):    
    preparedraw = time.time() 
    ana = pb.SegmentAnalyser(plant.mappedSegments()) 
    segOs = plant.getSegmentOrigins(-1, all = False)
    '''
    Lignification status, Survival, fine roots
    '''
    lignification = plant.getSubStatus() #[sO.lignificationStatus() for sO in segOs]
    aliveSegs = [sO.isAlive() for sO in segOs] #todo: removes
    is_fine_root = [sO.getParameter('is_fine_root') for sO in segOs]
    ana.addData('aliveSegs', aliveSegs)
    ana.addData('lignification', lignification)
    ana.addData('is_fine_root', is_fine_root)
    p_names = ['subType','aliveSegs','lignification','is_fine_root',"creationTime","id"] 
    for dd in data.keys():
        ana.addData(dd, data[dd])
        p_names.append(dd)
    #ana.filter('alive', 1)
    initdraw = time.time()
    vp.plot_roots(ana, show,p_names, win_title = title_, render = not saveOnly)
    #print("drawing: %s, %s" % (int(initdraw-preparedraw), int(time.time() - initdraw)), end=", ")
    

long_root_types = np.array([1,2,3,4,5])
fine_root_types = np.array([6,7,8,9])
subtypes = max(fine_root_types) #max(long_root_types)
orders = {'main' : [2], 'sub' : [3], 'subsub' : [4,5]}


kr =[viewer_conductivities.convert_radial(4.0e-7), 
      viewer_conductivities.convert_radial(7.0e-8),
      viewer_conductivities.convert_radial(4.0e-8)] # suberization status m s-1 MPa-1 => [1 / day]

# for K in m4 s-1 MPa-1 => [cm3 / day]
Kax_a =  {'B' : viewer_conductivities.convert_axial(0.04749/100.), 
          'D' : viewer_conductivities.convert_axial(0.1539/100.), 
          'E' : viewer_conductivities.convert_axial(0.02622/100.)}
Kax_b =  {'B' : 2.06437, 'D' : 2.2410, 'E' : 1.98847}

    
def run_benchmark(xx, genotype = 'B', rep_input = -1, simSoil = 50, extraName = "default", scenario = 0): #llambdao, kko,
    
    fast_imfp, fast_mfp = vg.create_mfp_lookup(sp, wilting_point = wilting_point, n = 1501)  # needed for Perirhizal class
    start_time = time.time()
    output = []
    file_path =  CPlantBox_dir + '/experimental/wine_fichtl/rsml/RSML_year1/'
    file_names = [name for name in os.listdir(file_path) if name[0] == genotype ] # repetitions of same genotype
    doBiCGSTAB = False

    data_file = file_names[rep_input]
    reps = 1# len(file_names)
    doVTP = (rep_input == 0)
    for rep in range(reps):
        N = 50
        if doProfile:
            N = 15
        outputs_12 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0
                }
        outputs_349 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0
                }
        outputs_50 = {
                'num':[0. for i in range(subtypes)],
                'length':[0. for i in range(1,subtypes)],
                'ratio':0,
                'kde_main':[0. for i in range(50)],
                'kde_sub':[0. for i in range(50)],
                'kde_subsub':[0. for i in range(50)]
                }
                
        outpouts_mean = {
            'year'+str(i + 1): copy.deepcopy(outputs_12) if i < 2 else copy.deepcopy(outputs_349) if i < 49 else copy.deepcopy(outputs_50)
            for i in range(N)
        }
        
        
        
        plant = MappedPlantPython(1) #PlantPython() # pb.MappedPlant() #

        # Open plant and root parameter from a file
        plant.readParameters( CPlantBox_dir + '/experimental/wine_fichtl/results/xmlFiles/' + genotype + "-wineV2.xml")
        plant.setSoilGrid(picker)
        plant.setGeometry(soilSpace) 


        no_thin = 500
        thin_root_P = [[0.],[0.05],[1.]]
        thin_root_Page = [0., 1225 * 5, 1225 * 10]
        ps = plant.getOrganRandomParameter(1)[0]
        ps.Lmax_unsuberized = 4.
        ps.Lmax_suberized = 10.
        #ps.seedPos = pb.Vector3d(0.,0.,-10.) # useless
        for ii, pp in enumerate(plant.getOrganRandomParameter(2)):            
            if ii == 0:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD * 50.#params['ldelays0'] #200#*5
                pp.successorP = [[[1.]]]# [[params['successorP0']]]
                pp.successorNo = [10]#[params['successorNo0']] 
                pp.successorST = [[[2]]]
                pp.successorOT = [[[2]]]
                
            elif ii == 1:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays = yr_to_BEDD *0.5 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[[1]]]#[[params['successorP0']]]
                pp.successorNo = [1] #[params['successorNo0']] 
                pp.successorST = [[[2]]]
                pp.successorOT = [[[2]]]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
                
            elif ii == 2:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50] 
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.5],[0.5],[0.1]], thin_root_P] #, [0.1, 0]]#[[params['successorP0']]]
                pp.successorST = [ [[ii + 1],[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #[[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorOT = [ [[2],[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[0., 1225 * 2, 1225 * 100],thin_root_Page]
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 26.8 
                pp.rlt_winter_min =  3.169
                pp.tropismN = 0.5  
                
                ratioChange = pp.ln / pp.dx
                pp.ln = pp.dx
                pp.lns = 0
                pp.successorNo = [pp.successorNo[0],int(pp.successorNo[1]/ratioChange)] 
                pp.successorP = [ [[0.5/ratioChange],[0.5/ratioChange],[0.1/ratioChange]], thin_root_P]
                
            elif ii == 3:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50]  #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.5],[0.5],[0.5]],thin_root_P] # [[0.5, 0.], [0.05, 0]]#[[params['successorP0']]]
                pp.successorST = [ [[ii + 1],[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #
                pp.successorOT = [ [[2],[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[0., 1225 * 2, 1225 * 10],thin_root_Page]
                # pp.r *= 2
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                #pp.k_survive = 1.367 
                #pp.lambda_survive =  4.995 
                pp.rlt_winter_max = 13
                pp.rlt_winter_min =  0
                
                ratioChange = pp.ln / pp.dx
                pp.ln = pp.dx
                pp.lns = 0
                pp.successorNo = [pp.successorNo[0],int(pp.successorNo[1]/ratioChange)] 
                pp.successorP = [ [[0.5/ratioChange],[0.5/ratioChange],[0.5/ratioChange]], thin_root_P]
                
            elif ii == 4:      
                pp.ldelay_v  = [0,0]
                pp.ldelays_v = [yr_to_BEDD , yr_to_BEDD*50] #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorNo = [1,no_thin] #[params['successorNo0']] 
                pp.successorP = [ [[0.05]],thin_root_P] #[[0.05, 0.], [0.05, 0.]]#[[params['successorP0']]]
                pp.successorST =  [ [[ii + 1],[ii + 1]], [[ii + 4],[ii + 4],[ii + 4]]] #[[ii + 1, ii + 4], [ii + 1, ii + 4]]
                pp.successorOT = [ [[2],[2]], [[2],[2],[2]]] 
                pp.successorP_age = [[],thin_root_Page]
                #pp.k_survive = 1.261
                #pp.lambda_survive = 1.367
                #pp.rlt_winter_max = 13.9 
                #pp.rlt_winter_min =  0.
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 13/2
                pp.rlt_winter_min =  0
                
                ratioChange = pp.ln / pp.dx
                pp.ln = pp.dx
                pp.lns = 0
                pp.successorNo = [pp.successorNo[0],int(pp.successorNo[1]/ratioChange)] 
                pp.successorP = [ [[0.05/ratioChange]], thin_root_P]
                
            elif ii == 5:      
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays =yr_to_BEDD * 50 #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [thin_root_P]#[[params['successorP0']]]
                #pp.successorNo = [no_thin] #[params['successorNo0']] 
                pp.successorST = [[[ii + 4],[ii + 4],[ii + 4]]]
                pp.successorOT = [[[2],[2],[2]]]
                pp.successorP_age = [thin_root_Page]
                #pp.k_survive = 1.261
                #pp.lambda_survive = 1.367
                #pp.rlt_winter_max = 13.9 
                #pp.rlt_winter_min =  0.
                pp.k_survive = 1.261
                pp.lambda_survive = 9.426
                pp.rlt_winter_max = 13/2/2
                pp.rlt_winter_min =  0
                
                ratioChange = pp.ln / pp.dx
                pp.ln = pp.dx
                pp.lns = 0
                pp.successorNo = [int(no_thin/ratioChange)]
                
            if (ii <= 5):   
                pp.a = 0.093/2
                pp.a_gr = 0.083/2/yr_to_BEDD
                
                
            else: # fine roots
                pp.r = 1. # don t know the value, i just want something high
                
                # was not used for r300
                pp.ldelay  = 0*yr_to_BEDD
                pp.ldelays =yr_to_BEDD  #yr_to_BEDD * params['ldelays0'] #200#*5
                pp.successorP = [[[1.]]]#[[params['successorP0']]]
                pp.successorNo = [1]#no_thin] #[params['successorNo0']] 
            
        saveall = ((simSoil == 50) and (scenario == 0))
        if saveall:     
            plant.writeParameters("./results/xmlFiles/"+genotype + "-wineV3.xml")
            
        '''
        start simulation
        '''
        # select randomly one of the rsml files of year 1
        plant.initialize_static(file_path + data_file, [0, 1])  # 0 is shoot, 1 are static roots

        # the static laterals 2 and 3 are replaced with the growing lateral 2
        ld, ld1 = plant.set_identical_laterals([0, 1], [1, 2, 3], 2)
        # plt.hist(np.array(ld1)/yr_to_BEDD, density = False, bins=30)
        # plt.title("Creation time of the main roots")
        # plt.show()
        # plt.hist(np.array(ld)/yr_to_BEDD, density = False, bins=30)
        # plt.title("Creation time of the main roots")
        # plt.show()

        plant.initialize_static_laterals()
        # plant.disableExtraNode() # TODO? troubleshoot the statis start for a mappedplant
        #plant.betaN = 5000
          
        all_real_lengths = []
        all_real_subtypes = []
        
        all_lengths = []
        all_ages = []
        all_subtypes = []
        all_alive = []
        SUFs, RLDs = [], []
        # all_alive_long = []
        # all_alive_short = []
        dt = yr_to_BEDD  # ~1 yr
        plant.do_simulate(0, False) # just to store data in the mapped segment
        
        if doVTP and ((simSoil == 50) and (scenario == 0)):
            get3Dshape(plant,title_ = "./results/part1/vtp/"+extraName+'/'+genotype+"0", data = {}, saveOnly = True)
            
        param = PlantHydraulicParameters()
        #param.set_kr_const(1.728e-4)          
        #param.set_kx_const(4.32e-2 )
        param.set_kr_suberize_dependent(kr)          
        param.set_kx_radius_dependent([Kax_a[genotype],Kax_b[genotype]])
        hm = HydraulicModel_Doussan(plant, param) #HydraulicModel_Meunier(plant, param) # _large
        hm_meunier = HydraulicModel_Meunier(plant, param)
        hm.doBiCGSTAB = doBiCGSTAB
        #hm_meunier.doBiCGSTAB = doBiCGSTAB
        #peri = Perirhizal(plant)
        assert len(plant.get_nodes()) == (len(plant.get_segments()) +1)
        reptimeBU = time.time()
        timing = defaultdict(float)
        total_loop_time = 0.0

        for i in range(N):
            iter_start = tic()
            print('age', i, end=", ", flush=True)

            # ---------------- survival test ----------------
            t0 = tic()
            plant.survivalTest()
            timing['survival'] += toc(t0)

            # ---------------- growth ----------------
            t0 = tic()
            plant.do_simulate(dt, False)
            timing['growth'] += toc(t0)
            
            # ---------------- RLD ----------------
            t0 = tic()
            ana = pb.SegmentAnalyser(plant.mappedSegments())
            minZ = max(np.fromiter(hm.ms.cell2seg.keys(), dtype=int)) + 1
            if saveall: 
                timing['rld'] += toc(t0)
                t0 = tic()
                RLDs.append(
                    ana.distributionFast("lengths", plant) #("length", 0., minZ, int(-minZ), False)
                )
                timing['rld2'] += toc(t0)

            # ---------------- SUF ----------------
            if (( i + 1 ) == simSoil) or saveall: 
                # t0 = tic()
                try:
                    hm.update(i * dt, dz, abs(int(minZ)), ana, volumes, fast_imfp, fast_mfp)
                    suf_ = hm.get_suf() #get_suf(i * dt)
                    #suf_meunier = hm_meunier.get_suf(i * dt)
                    assert suf_.min() >= 0 
                except:
                    print("issue with suf computation")
                    get3Dshape(
                        plant,
                        title_="./results/part1/vtp/" + extraName + '/' + genotype + str(i+1) +"_"+ str(rep_input) +"_error",
                        saveOnly=True,data ={'SUF': suf_ #,'psiXyl':hm.psiXyl
                                            },show = 'SUF',# 
                    )
                    #print('np.unique(hm.psiXyl)',np.unique(hm.psiXyl))
                    #hm.test()
                    plant.writeRSML("./results/outputSim/" + extraName + '/' + genotype + str(i+1) +"_"+ str(rep_input) +"_error.rsml")
                    raise Exception


                timing['suf'] += toc(t0)
                t0 = tic()
                ana.addData('SUF', suf_)
                #ana.addData('SUF_meunier', suf_meunier)
                suf1d = np.array(ana.distributionFast("SUF", plant))
                #suf1d_meunier = np.array(ana.distributionFast("SUF_meunier", plant))
                SUFs.append(suf1d)
                timing['suf_distrib'] += toc(t0)
            
            '''
            fig, ax1 = plt.subplots()
            depth_array = np.array([-ii for ii in range(len(suf1d))])
            #ax1.plot(suf1d_meunier, depth_array,linewidth=6, label = "SUF_menier [-]",color='g')
            ax1.plot(suf1d, depth_array, label = "SUF_dussan [-]",color='r')
            plt.show()
            #print(suf1d_meunier)
            #print(suf1d)
            #raise Exception
            '''
            # ---------------- processing ----------------
            if saveall: 
                t0 = tic()
                orgs_all = plant.getOrgans(2, False)
                orgs_ = []

                for oo in orgs_all:
                    distance = 0.
                    ooc = oo
                    stc = ooc.param().subType
                    while stc > 2:
                        distance += ooc.getParent().getLength(ooc.parentNI)
                        ooc = ooc.getParent()
                        stc = ooc.param().subType
                    if distance < 175:
                        orgs_.append(oo)

                all_ages = np.array([org.getAge()/yr_to_BEDD for org in orgs_ if org.isAlive()]) #years
                all_alive.append(np.array([org.isAlive() for org in orgs_ if org.isAlive()]))
                all_subtypes.append(np.array([org.param().subType for org in orgs_ if org.isAlive()]))
                all_lengths.append(np.array([org.getLength() for org in orgs_ if org.isAlive()]))

                all_real_subtypes.append(np.array([org.param().subType for org in orgs_all if org.isAlive()]))
                all_real_lengths.append(np.array([org.getLength() for org in orgs_all if org.isAlive()]))

                is_fine_roots = np.isin(all_real_subtypes[i], fine_root_types)
                len_fine = sum(all_real_lengths[i][is_fine_roots]) if np.any(is_fine_roots) else 0.0
                len_long = sum(all_real_lengths[i][~is_fine_roots])
                value = len_fine / (len_long + len_fine)

                timing['processing'] += toc(t0)

                # ---------------- VTP ----------------
                if doVTP:
                    t0 = tic()
                    get3Dshape(
                        plant,
                        title_="./results/part1/vtp/" + extraName + '/' + genotype + str(i+1),
                        saveOnly=True, show = 'SUF',
                        data ={'SUF': suf_ #, 'psiXyl':hm.psiXyl
                              }
                    )
                    timing['vtp'] += toc(t0)

            # ---------------- SOIL ----------------
            t0 = tic()
            if (( i + 1 ) == simSoil): 
                name_str = "part1/vtp/" + extraName + '/' + genotype + '_soil_' + str(scenario) + "_"  + str(i+1) #'outputSim/'+extraName+'/'+ genotype + str(rep) +'_soil' + str(i + 1) 
                soil_simulation = do_soil_simulation(plant, hm, hm_meunier,minZ, suf1d, #np.flip(), #flip suf, as the plant picker is ordered opposit to the soil picker
                                                     name_str, scenario, simSoil)
                print('/results/outputSim/'+extraName+'/'+ genotype + str(rep) + '_soil_' + str(scenario) + "_"  + str(i + 1) +'.pkl')
                with open(CPlantBox_dir + '/experimental/wine_fichtl/results/outputSim/'+extraName+'/'+ genotype + str(rep) + '_soil_' + str(scenario) + "_"  + str(i + 1) +'.pkl','wb') as f:
                     pickle.dump(soil_simulation,f, protocol=pickle.HIGHEST_PROTOCOL)
                        
                if ((simSoil != 50) or (scenario != 0)): 
                    print('Goodbye')
                    return 0
            
            timing['soil'] += toc(t0)
            
            # ---------------- iteration total ----------------
            iter_time = toc(iter_start)
            timing['total'] += iter_time
            total_loop_time += iter_time

            print(f"iter_time: {iter_time:.2f}s", [(tt, np.round(timing[tt],1)) for tt in timing])
                            

        if saveall:     
            # print('postprocessing')
            #print("--- %s seconds for plant development---" % (time.time() - start_time),rep)
            outpouts_mean['SUF'] = SUFs
            outpouts_mean['RLDs'] = RLDs
            for year in range(N):
                '''
                Lengths
                '''
                for st in range(1,subtypes): # from 2 to 5
                    value = sum(all_lengths[year][all_subtypes[year] == (st + 1)])
                    outpouts_mean['year'+str(year+1)]['length'][st-1] = value  
                '''
                Num
                '''
                for st in range(subtypes): # from 2 to 5            
                    value = sum(all_alive[year][all_subtypes[year] == (st + 1)])
                    outpouts_mean['year'+str(year+1)]['num'][st] = value
                #outpouts_mean['year'+str(year+1)]['num'][subtypes] = sum(all_alive[year][all_subtypes[year] > subtypes])
                '''
                Ratio
                ''' 
                len_fine = 0
                is_fine_roots = np.isin(all_real_subtypes[year],fine_root_types)
                if sum(is_fine_roots) > 0:
                    len_fine = sum(all_real_lengths[year][is_fine_roots])
                len_long = sum(all_real_lengths[year][np.invert(is_fine_roots)]) # ignore the length of 1 as seem different between xml and xlsx
                value = len_fine/(len_long + len_fine)
                outpouts_mean['year'+str(year+1)]['ratio'] = value     


            '''
            Age distribution at last time step
            '''
            idy = N - 1


            for key, value in orders.items():
                # List to store KDE y-values for each FileName

                x = xx[key]

                root_age = all_ages[np.isin(all_subtypes[idy],value)] #years
                try:
                    if len(root_age) == 0:
                        y = np.zeros(len(x))
                    elif len(root_age) == 1:
                        y = np.zeros(len(x))
                        indx_y = np.where(abs(x - root_age[0]) == min(abs(x - root_age[0])))[0][0]
                        y[indx_y] = 1.
                    else:
                        kde = gaussian_kde(root_age)
                        y = kde(x)
                except:
                    print('issue root_age',root_age)
                    raise Exception

                outpouts_mean['year'+str(N)]['kde_' +key] = y

            def flatten_values(obj):
                if isinstance(obj, dict):
                    for v in obj.values():
                        yield from flatten_values(v)
                elif isinstance(obj, (list, tuple, np.ndarray)):
                    for v in obj:
                        yield from flatten_values(v)
                else:
                    yield obj
            output.append(outpouts_mean)#list(flatten_values(outpouts_mean)))

            print('output',output)
            print("--- %s seconds for plant development---" % ((time.time() - start_time)/(rep+1)))
            # print("--- %s seconds for total---" % (time.time() - start_time),'at',(rep+1),'/',(reps))

    #if getdict:
    #    return outpouts_mean
        
    # print('outpouts_mean',outpouts_mean)
    #f = [i for s in output for i in s]
    return output #np.array(f)
        


    
if __name__ == '__main__': 
    #     jid=$(sbatch --nodelist=node10 --parsable wineSoil_24cpu256.sh B 0 "$i" "$1")
    genotype = sys.argv[1]
    rep = int(sys.argv[2])
    
    if len(sys.argv) > 3:
        simSoil = int(sys.argv[3])
    else:
        simSoil = 50
    scaling = np.array([ 1,  2,  3,  5,  8, 13, 20, 32, 50]) # 9
    
    scenario = 0
    if len(sys.argv) > 4:
        scenario = int(sys.argv[4])

    extraName = "defaultsoil"
    if len(sys.argv) > 5:
        extraName = sys.argv[5]


    print('sys.argv',sys.argv, genotype, rep, simSoil,scenario, extraName, flush = True)

    if ((simSoil == 50) and (scenario == 0)):   
        directory ='./results/outputSim/'+extraName
        os.makedirs(directory, exist_ok=True)
        directory ='./results/part1/vtp/'+extraName
        os.makedirs(directory, exist_ok=True)

        with open(CPlantBox_dir + '/experimental/wine_fichtl/results/objectiveData/measurements'+ genotype +'InitXX.pkl','rb') as f:
            xx = pickle.load(f)
    else:
        xx = None 

    doProfile = False
    if doProfile:
        import cProfile
        import pstats, io
        pr = cProfile.Profile()
        pr.enable()

    output = run_benchmark(xx, genotype, rep, simSoil, extraName,scenario)
    if doProfile:
        pr.disable()
        filename = './results/profile'+str(rank)+'.prof' 
        s = io.StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('tottime')
        ps.dump_stats(filename)

    saveall = ((simSoil == 50) and (scenario == 0))
    if ((simSoil == 50) and (scenario == 0)):     
        with open(CPlantBox_dir + '/experimental/wine_fichtl/results/outputSim/'+extraName+'/'+ genotype + str(rep) +'.pkl','wb') as f:
             pickle.dump(output,f, protocol=pickle.HIGHEST_PROTOCOL)

    #comm.Barrier()
    #if rank == 0:
    #    import gatherAndPlotOutputs
    #    gatherAndPlotOutputs.plot_all()

