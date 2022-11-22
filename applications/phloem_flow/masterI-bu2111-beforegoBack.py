""" water movement within the root (static soil) """


import sys; 
#directoryN = "/"+sys.argv[0].split('.')[0]+"/"
sys.path.append("../.."); sys.path.append("../../src/python_modules")
CPBdir = "../.."
#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs
import gc
from phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np


#import matplotlib.pyplot as plt
import time
import vtk_plot as vp
isCluster = (os.environ['HOME'] == '/home/m.giraud')



#       qr, qs, alpha, n, ks (in [cm/d])
#vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2
def stairs(t):
    hours = t%1
    coef = int(hours<=(12/24))
    return coef

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh


def weather(simDuration, Qmax):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    
    Qmin = 0;# Qmax = 1000e-6 #458*2.1
    Tmin = 10; Tmax = 22
    specificHumidity = 0.0097
    Pair = 1010.00 #hPa
    thetaInit = 30/100

    coefhours = stairs(simDuration)#sinusoidal(simDuration)
    TairC_ = Tmin + (Tmax - Tmin) * coefhours
    Q_ = Qmin + (Qmax - Qmin) * coefhours
    cs = 350e-6 #co2 paartial pressure at leaf surface (mol mol-1)
    #RH = 0.5 # relative humidity
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    RH = 0.3#qair2rh(specificHumidity, es, Pair)
    
    pmean =-100#theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,
                    'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
    #print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        


def setKrKx_xylem(TairC, RH): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    mu = mu * hPa2cm #hPa d to cmh2o d 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #radius of xylem type^4 * number per bundle
    rad_x_l_1   = (0.0015 **4) * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_x_s_1   = (0.0017 **4) * 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_x_r0_1  = (0.0015 **4) * 4    
    rad_x_r12_1 = (0.00041**4) * 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_x_r3_1  = (0.00068**4) * 1      

    # axial conductivity [cm^3/day]        
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8)  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8) 
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) # 4.32e-1

    #radial conductivity [1/day],
    kr_l  = 3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 6.37e-5 * hPa2cm 
    kr_r1 = 7.9e-5  * hPa2cm 
    kr_r2 = 7.9e-5  * hPa2cm  
    kr_r3 = 6.8e-5  * hPa2cm 
    l_kr = 100 #cm
                                     
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ,kr_s,kr_s ],[kr_l,kr_l,kr_l,kr_l]], kr_length_=l_kr)
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s,kz_s,kz_s ],[kz_l,kz_l,kz_l,kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm

    r.psi_air = p_a #*MPa2hPa #used only with xylem
    
def setKrKx_phloem(): #inC

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #numPerBundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1
    a_ST    = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    a_ST = np.array([np.array(xi) for xi in a_ST], dtype=object)#*2
    #radius of phloem type^4 * number per bundle
    rad_s_l   = numL  * (a_ST[2][0] ** 4)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    rad_s_s   = numS  * (a_ST[1][0] ** 4) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    rad_s_r0  = numr0 * (a_ST[0][0] ** 4) #* 4    
    rad_s_r12 = numr1 * (a_ST[0][1] ** 4) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    rad_s_r3  = numr3 * (a_ST[0][2] ** 4) #* 1      

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 1 #*1000#Thompson 2003a
    kz_l   = VascBundle_leaf * rad_s_l   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * rad_s_s   * np.pi /8 * beta
    kz_r0  = VascBundle_root * rad_s_r0  * np.pi /8 * beta
    kz_r12 = VascBundle_root * rad_s_r12 * np.pi /8 * beta
    kz_r3  = VascBundle_root * rad_s_r3  * np.pi /8 * beta
    
    #print([[kz_r0,kz_r12,kz_r12,kz_r3],[kz_s,kz_s ],[kz_l]])
    #raise Exception
    #radial conductivity [1/day],
    kr_l  = 0.#3.83e-4 * hPa2cm# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 = 0#5e-2
    kr_r1 = 0#5e-2
    kr_r2 = 0#5e-2
    kr_r3 = 0#5e-2
    l_kr = 0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s, kr_s,kr_s],[kr_l,kr_l,kr_l,kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_l,kz_l, kz_l,kz_l],[kz_l,kz_l,kz_l,kz_l]])
    
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi#(0.00039 **2) #* 4    
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi# (0.00068**2) #* 1    
    
    Perimeter_s_l   = numL*VascBundle_leaf *(a_ST[2][0])* 2 * np.pi# (0.00025 **2)# * 2; rad_x_l_2   = (0.0005 **4) * 2   
    Perimeter_s_s   = numS *VascBundle_stem * (a_ST[1][0])* 2 * np.pi#(0.00019 **2) #* 3; rad_x_s_2   = (0.0008 **4) * 1     
    Perimeter_s_r0  = numr0 *VascBundle_root * (a_ST[0][0])* 2 * np.pi#(0.00039 **2) #* 4    
    Perimeter_s_r12 = numr1*VascBundle_root * (a_ST[0][1])* 2 * np.pi#(0.00035**2) #* 4; rad_x_r12_2 = (0.00087**4) * 1
    Perimeter_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2])* 2 * np.pi# (0.00068**2) #* 1  
    #print(a_ST[2][0],a_ST[1][0],a_ST[0][0],a_ST[0][1],a_ST[0][2])
    #r.a_ST = a_ST #to check for water equilibrium assumption
    #tot surface/np.pi of sieve tube  (np.pi added after)
    #r.a_ST_eqs = [[rad_s_r0,rad_s_r12,rad_s_r12,rad_s_r0],[rad_s_s,rad_s_s],[rad_s_l]]
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s,Across_s_s,Across_s_s],[Across_s_l,Across_s_l,Across_s_l,Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])

from pathlib import Path
def runSim(directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, UseRatiothresholdAux,
           Qmax_, thresholdSuc,
           GrRatio ,  maxLBud , budGR,L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D = 0.,kss=0.2,kaa=1.):
    
    useCWGr = True
    dt_lastWrote = time.time() - dt_write * 2
    dtSIM_lastWrote = -dtSIM_write*2
    temp_time  = time.time()
    outputsDict_array = {}
    outputsDict_float = {}
    #print("num ",thread,temp_time , start_time, "has started. This sampling run took %5.4f seconds." % (temp_time - start_time))
    doDecapitation = (nodeD > 0)
    directoryN = directoryN_
    strPRate_ = str(np.round(PRate_*10))#[:-2]
    strThA = str(np.round(thresholdAux*10))#[:-2]
    strThRatioA = str(np.round(RatiothresholdAux*10))#[:-2]
    
    strQ = str(np.round(Qmax_*1e6))#[:-2]
    strThS = str(np.round(thresholdSuc*10))#[:-2]
    
    strDecap = str(nodeD)#[0]
    
    dir4allResults ="_"+str(thread) #"_"+ strPRate_+ "_"+strThA+ "_"+strThRatioA + "_"+ strQ + "_"+strThS+"_"+strDecap+"_"+str(thread)
    dir4allResults=dir4allResults.replace(".", "o")
    def write_file_array(name, data):
        if doPrint :
            name2 = 'results'+ directoryN+ name+ dir4allResults+ '.txt'
            with open(name2, 'a') as log:
                log.write(','.join([num for num in map(str, data)])  +'\n')
        if doDict:
            if name not in outputsDict_array:
                outputsDict_array[name] = list()
            outputsDict_array[name].append(data)
        

    def write_file_float(name, data):
        if doPrint:
            name2 = 'results'+ directoryN+ name+ dir4allResults+ '.txt'
            with open(name2, 'a') as log:
                log.write(repr( data)  +'\n')
        if doDict:
            if name not in outputsDict_float:
                outputsDict_float[name] = np.array([])
            outputsDict_float[name] = np.append(outputsDict_float[name],data )

    def delete_file(name):
        if doPrint:
            name2 = 'results'+ directoryN+ name+dir4allResults+ '.txt'
            my_file = Path(name2)
            if my_file.is_file():
                os.remove(name2)
            else:
                print("file running missing",thread, name2)
    
    namef = 'results'+ directoryN+ "running"+ dir4allResults+ '.txt'
    my_file = Path(namef)
    if my_file.is_file():
        print("input", np.array([directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, UseRatiothresholdAux,
           Qmax_, thresholdSuc,
           GrRatio ,  maxLBud , budGR,L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D ]))
        print(namef)        
        print("runningfile already exists")
        raise Exception
    
    doPrintbu = doPrint
    doDictbu = doDict
    doDict = False
    doPrint = True
    write_file_float("running", thread)  
    write_file_array("input", np.array(["directoryN_","doVTP", "verbosebase",
           "PRate_", "thresholdAux", "RatiothresholdAux", "UseRatiothresholdAux",
           "Qmax_", "thresholdSuc",
           "GrRatio" ,  "maxLBud" , "budGR","L_dead_threshold" ,
           "nodeD", "thread",  
           "testTime", "dtBefore", "dtAfter", "start_time", "dt_write","dtSIM_write",
           "doPrint" ,"doDict", "auxin_D","kss", "kaa" ])) 
    write_file_array("input", np.array([directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, UseRatiothresholdAux,
           Qmax_, thresholdSuc,
           GrRatio ,  maxLBud , budGR,L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D,kss,kaa ])) 
    # if RatiothresholdAux == 0 and UseRatiothresholdAux and activeAtThreshold_auxin:
    #     print(thread, RatiothresholdAux,UseRatiothresholdAux ,activeAtThreshold_auxin)
    #     print("issue ratio threshold")
    #     write_file_float("stop", thread)
    #     raise Exception
    doPrint = doPrintbu
    doDict = doDictbu
    weatherInit = weather(0, Qmax_)
    
    simDuration = 0#simInit # [day] init simtime
    depth = 60
    dt = dtBefore #1h
    verbose = False

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/plant/"
    name = "UQ_simple_stem_bud"#"smallPlant_mgiraud"#"morning_glory_UQ"#

    pl.readParameters(path + name + ".xml")


    #raise Exception
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    #pl.activeAtThreshold_auxin = activeAtThreshold_auxin
    #pl.activeAtThreshold = activeAtThreshold_suc
    pl.initialize(verbose = False)#, stochastic = False)
    
    while (len(pl.getOrgans(3, False)) ==0 ) or (pl.getOrgans(3, False)[0].getNumberOfLinkingNodes() <3):
        pl.simulate(dt, False)#, "outputpm15.txt")
        simDuration += dt
        
    simMax = 100#simDuration + testTime #if !doDecapitation

    """ Coupling to soil """
    min_b = [-3./2, -12./2, -61.]#distance between wheat plants
    max_b = [3./2, 12./2, 0.]
    cell_number = [6, 24, 61]#1cm3? 
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = -100#weatherInit["p_mean"]#mean matric potential [cm] pressure head

    # sx = s.getSolutionHead()  # inital condition, solverbase.py
    p_mean = weatherInit["p_mean"]#-187
    p_bot = p_mean + depth/2
    p_top = p_mean - depth/2
    sx = np.linspace(p_top, p_bot, depth)
    picker = lambda x,y,z : max(int(np.floor(-z)),-1) #abovegroud nodes get index -1

    #picker = lambda x, y, z: s.pick([x, y, z])    
    pl.setSoilGrid(picker)  # maps segment


    """ Parameters phloem and photosynthesis """
    global r
    r = PhloemFluxPython(pl,psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#

    setKrKx_phloem()
    r.g0 = 8e-2
    r.VcmaxrefChl1 =2
    r.VcmaxrefChl2 =7
    r.a1 = 2#0.8/0.2
    r.a3 = 2.2
    r.Rd_ref = 0#2e-6
    r.alpha = 0.27
    r.theta = 0.51
    r.k_meso = 1e-3#1e-4
    r.setKrm2([[0]])#2e-4
    r.setKrm1([[1.3e-3]])#3e-03#([[2.5e-2]])
    #r.setRhoSucrose([[0.51],[0.65],[0.56]])
    rho_org = [[1.34],[1.44],[1.56]]#g C/gDW?
    density = 0.17 #g DW/cm3?
    #density and rho in mmol suc or carbon?
    rho_org = np.array([np.array(xi) for xi in rho_org])*density/2 #/2 => glucose to sucrose
    r.setRhoSucrose(rho_org)
    grRate = [[4.,2.,1.0,4.],[2.,2.,3.],[3.,3.,3.,3.]]
    grRate = np.array([np.array(xi) for xi in grRate],dtype=object)
    grRate[1:] *= GrRatio #[1:]

    r.setRmax_st(grRate)
    r.KMfu = 0.2
    #r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    #r.k_gr = 1#0
    r.sameVolume_meso_st = True
    r.sameVolume_meso_seg = False
    r.withInitVal =True
    r.initValST = 0.#0.15
    r.initValMeso = 0.#0.2
    r.beta_loading = 0.0#15
    r.Vmaxloading = 0.3 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 0.2
    r.Gr_Y = 0.75
    r.CSTimin = 0.05#
    #r.surfMeso=0.0025
    r.cs = weatherInit["cs"]

    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-14
    r.rtol = 1e-10
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 60.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4

    #turn off effect of water limitation
    r.psiMin = -10000000000*(1/0.9806806)
    r.psiMax = -1000000000*(1/0.9806806)

    """ for post processing """
    structSumInit = 0
    orgs_all = r.plant.getOrgans(-1, True)

    for org in orgs_all:
        if org.organType() < 2:
            raise Exception("ot < 3")
        structSumInit += org.orgVolume(-1,False) * r.rhoSucrose_f(
            int(org.getParameter("subType")-1),
                           org.organType())

    AnSum = 0
    Nt = len(r.plant.nodes) 
    Q_Rmbu      = np.array([0.])

    Q_Grbu      = np.array([0.])
    Q_Gr4bu      = np.array([0.])
    Q_Gr3bu      = np.array([0.])
    Q_Gr2bu      = np.array([0.])
    Q_Grbuth      = np.array([0.])
    Q_Gr4buth      = np.array([0.])
    Q_Gr3buth      = np.array([0.])
    Q_Gr2buth      = np.array([0.])

    Q_Exudbu    = np.array([0.])
    Q_Rmmaxbu   = np.array([0.])
    Q_Grmaxbu   = np.array([0.])
    Q_Exudmaxbu = np.array([0.])
    Q_STbu      = np.array([0.])
    Q_Parbu      = np.array([0.])
    Q_mesobu    = np.array([0.])
    volSegbu =  np.array([0.])
    NOrg = r.plant.getNumberOfOrgans()
    delta_ls_bu = np.full(NOrg, 0.)
    delta_ls_max = 0


    Ntbu = 1
    Q_in  = 0
    Q_out = 0


    orgs_all = r.plant.getOrgans(-1, True)
    volOrgbu = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
    volOrg = np.array([org.orgVolume(-1,False) for org in orgs_all]) 
    ot_orgs = np.array([org.organType() for org in orgs_all])
    st_orgs = np.array([org.getParameter("subType") for org in orgs_all])

    volOrgini = np.array([org.orgVolume(-1,False) for org in orgs_all])

    volOrgi_th = 0.
    lenOrgbu = np.array([org.getLength(False) for org in orgs_all]) 
    lenOrg = np.array([org.getLength(False) for org in orgs_all]) 
    lenOrgi_th = 0.
    Orgidsbu = np.array([org.getId() for org in orgs_all])
    Orgids = np.array([org.getId() for org in orgs_all]) #true:realized
    #raise Exception
    ö=0

                                 
    volOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)]),
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)]), 
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])]) 
      
    volOrgini2 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])
    volOrgini3 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])
    volOrgini4 = np.array([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])
    sucOrgini_type =  np.array([sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(2, True)])*r.rhoSucrose_f(0,2),
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(3, True)])*r.rhoSucrose_f(0,3), 
                                sum([org.orgVolume(-1,False) for org in r.plant.getOrgans(4, True)])*r.rhoSucrose_f(0,4)]) 
    typeOrg_unit =  np.array([org.organType() for org in r.plant.getOrgans(-1, True)])
    idOrg_unit =  np.array([org.getId() for org in r.plant.getOrgans(-1, True)])
    sucOrg_unit =  np.array([org.orgVolume(-1,False)*r.rhoSucrose_f(int(org.getParameter("subType")-1),org.organType()) for org in r.plant.getOrgans(-1, True)])
    sucOrgini_unit = sucOrg_unit
    #print(volOrgini2_type, volOrgini_type)
    sucOrg_type = sucOrgini_type
    #sucOrg2_type =sucOrgini2_type
    #volOrg2_type = volOrgini2_type

    volOrg_type = volOrgini_type
    orgs_roots = r.plant.getOrgans(2, True)
    orgs_ln = np.array([])
    orgs_st = np.array([])
    for org_ in orgs_roots:
        if org_.getParameter("subType") == 2 :
            orgs_ln= np.append(orgs_ln,len(org_.param().ln) ) 
            orgs_st= np.append(orgs_st,org_.getParameter("subType")  ) 


    orgs_all = r.plant.getOrgans(-1, True)
    volOrgbu = np.array([np.full( (len(org.getNodeIds())-1),org.getId())  for org in orgs_all], dtype=object) 


    #1h for 1d when dxMin = 0.3

    AnSum = 0
    Q_ST_init = np.array([])
    Q_meso_init  = np.array([])
    Q_AuxinInit  = np.array([])
    Q_Gr4bu =Q_Gr3bu=Q_Gr2bu=[0]
    deltasucorgbu = np.array([])
    AnSum = 0
    Nt = len(r.plant.nodes) 
    Q_Rmbu      = np.array([0.])

    Q_Grbu      = np.array([0.])
    Q_Gr4bu      = np.array([0.])
    Q_Gr3bu      = np.array([0.])
    Q_Gr2bu      = np.array([0.])
    Q_Grbuth      = np.array([0.])
    Q_Gr4buth      = np.array([0.])
    Q_Gr3buth      = np.array([0.])
    Q_Gr2buth      = np.array([0.])

    Q_Exudbu    = np.array([0.])
    Q_Rmmaxbu   = np.array([0.])
    Q_Grmaxbu   = np.array([0.])
    Q_Exudmaxbu = np.array([0.])
    Q_STbu      = np.array([0.])
    Q_Parbu     = np.array([0.])
    Q_mesobu    = np.array([0.])


    r.doTroubleshooting = False

    r.auxin_threshold = thresholdAux
    r.auxin_D = auxin_D #e-6#muM /d  #e-6#3e-7
    r.auxin_P = PRate_ #muM /d  #*1e-6# 6e-7
    r.auxin_alpha = 0.002727303 #0.002704272#0.00245669#0.02 #5#e-6# 2.4e-3
    r.initValAuxin = 0 #0.2
    r.deleteAtRootTip = True
    r.plant.useCWGr = useCWGr 
    InAuxin = 0
    #r.canStartActivating = False
    r.CSTthreshold = thresholdSuc
    r.StopLoss = True
    if thread == 0 :
        print(r.L_dead_threshold, r.plant.maxLBud)
    r.L_dead_threshold = L_dead_threshold
    r.plant.maxLBud = maxLBud
    r.plant.budGR = budGR
    orgs_stems = r.plant.getOrgans(3, True)
    r.computeBerth = lambda ss_, aa_: ((aa_+kaa)/(ss_*10+kss))*(1-(0.15/(ss_*10+kss)))*2 
    
    def killChildren(orgToKill):
        toFill = orgToKill.getNodeIds()
        orgToKill.alive = False
        orgToKill.budStage = -1
        if orgToKill.getNumberOfChildren() > 0:
            toFill_ = [killChildren(orgToKill.getChild(ni)) for ni in range(orgToKill.getNumberOfChildren())]
            toFill_ = [item for sublist in toFill_ for item in sublist]
            toFill = toFill + toFill_
        return toFill
    
    
    
    burninDuration = simDuration
    r.burnInTime = True#activeAtThreshold_auxin

    changedSimMax = False
    oo = 0
    while  simDuration < simMax:#
        #if(ö>=120):
        
        temp_time = time.time()
        #print("num ",thread,temp_time , start_time, "has started. This sampling run took %5.4f seconds." % (temp_time - start_time),'simDuration:',simDuration )
        numLNodes = r.plant.getOrgans(3, False)[0].getNumberOfLinkingNodes()
        Mstem = r.plant.getOrgans(3, False)[0]
        kids4distbase = np.array([Mstem.getChild(nkdb) for nkdb in range(Mstem.getNumberOfChildren())])
        tempstst = np.array([nkdb.getParameter("subType") for nkdb in kids4distbase ])
        kids4distbase = kids4distbase[tempstst == 2][-2:]
        distbase4decap = np.array([nkdb.getParent().getLength(nkdb.parentNI) for nkdb in kids4distbase])
        
        if  (numLNodes >= nodeD) and (not r.burnInTime)and (nodeD !=0) and ((distbase4decap[1] - distbase4decap[0])>=1): 
            org_ = r.plant.getOrgans(3, False)[0] #get first (==main) stem
            org_.active = False
            org_.budStage = -1
            toKil = np.array([org_.getNodeId(org_.getNumberOfNodes()-1)])#org_.getNodeId(org_.getNumberOfNodes()-2),
            kid_pni = np.array([oo.getNodeId(0) for oo in r.plant.getOrgans(-1, True)])
            if len(kid_pni) > 0:
                selectKids =np.where( kid_pni ==toKil[0] )[0]#np.concatenate((,np.where(kid_pni==toKil[1])[0]))
                dying_kids = np.array(r.plant.getOrgans(-1, True))[selectKids]
                k_ =[killChildren(ni) for ni in dying_kids]
                toFill_ = [item for sublist in k_ for item in sublist]
                toFill_= np.concatenate((toFill_,toKil))
                toKil  = np.array(np.unique(toFill_),dtype = np.int32)
            r.plant.node_Decapitate = toKil
            
            doDecapitation = False
            if not changedSimMax:
                simMax = simDuration + testTime #end 7 days after decapitation
                changedSimMax = True
                dt = dtAfter #1MIN

        if doDecapitation and (numLNodes > nodeD): 
            raise Exception("too many linking nodes")
        
        if (not doDecapitation) and (not changedSimMax) and (nodeD ==0) and (numLNodes ==8): 
            simMax = simDuration + testTime
            changedSimMax = True
            dt = dtAfter #1MIN
            
        if r.burnInTime:
            weatherX = weather(burninDuration, Qmax_)
        else:
            weatherX = weather(simDuration, Qmax_)
            

        r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"
            
            
        setKrKx_xylem(weatherX["TairC"], weatherX["RH"])
        
        #r.maxLoop = 3
        #r.doTroubleshooting = False
        r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
            verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
        #print(r.psiXyl)
        #raise Exception
        if False:#verbosebase:
            segIdx = r.get_segments_index(4)
            idleafBlade = np.where(np.array(r.plant.leafBladeSurface)[segIdx] > 0)
            print("An",np.mean(np.array(r.An)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#An.append
            print("Rd",r.Rd*1e6, "mumol CO2 m-2 s-1")#An.append
            print("Vc",np.mean(np.array(r.Vc)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#Vc.append
            #print("Vcmax",np.mean(r.Vcmax)*1e6, "mumol CO2 m-2 s-1")#Vc.append
            print("Vj",np.mean(np.array(r.Vj)[idleafBlade])*1e6, "mumol CO2 m-2 s-1")#Vj.append
            #print("Vjmax",np.mean(r.Vjmax)*1e6, "mumol CO2 m-2 s-1")#Vj.append
            print("gco2",np.mean(np.array(r.gco2)[idleafBlade]), "mol CO2 m-2 s-1")#gco2.append
            print("gh2o",np.mean(np.array(r.gco2)[idleafBlade])*1.6, "mol H2O m-2 s-1")#gco2.append
            print("cics",np.mean(np.array(r.ci)[idleafBlade])/r.cs,"mol mol-1")#cics.append
            print("ci",np.mean(np.array(r.ci)[idleafBlade]), "mol mol-1")#fw.append
            print("deltagco2",np.mean(np.array(r.deltagco2)[idleafBlade]), "mol mol-1")#fw.append
            print("fw",np.mean(np.array(r.fw)[idleafBlade]), "-")#fw
            
        ots = np.concatenate((np.array([0]), r.get_organ_types()))#per node
        leavesSegs = np.where(ots[1:] ==4)
        fluxes = np.array(r.outputFlux)
        fluxes_leaves = fluxes[leavesSegs]
        errLeuning = sum(r.outputFlux)
            
        """ dumux """    
        fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
        #r.doTroubleshooting = (ö == 22)
        #r.Ag4Phloem = np.full(len(r.Ag4Phloem),0)
        AnSum += np.sum(r.Ag4Phloem)*dt
        startphloem= simDuration
        endphloem = startphloem + dt
        stepphloem = 1
        
        # save_stdout = sys.stdout
        # sys.stdout = open('trash', 'w')
        verbose_phloem = False
        
        filename = "results"+ directoryN +"inPM_"+str(ö)+dir4allResults+ '.txt'
        #print("startpm")
        
        r.useStemTip = True
        doStartPM = True
        if doStartPM:
            r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
        else:
            r.Q_out = np.full(Nt*10,0.)
            r.C_ST = np.full(Nt,0.)
            r.vol_ST = np.full(Nt,1.)
            r.Fl = np.full(Nt,0.)
            r.vol_Meso = np.full(Nt,1.)
            r.SucSTLost= np.full(Nt,0.)
            r.SucMesoLost= np.full(Nt,0.)
            r.manualAddST= np.full(Nt,0.)
            r.manualAddMeso= np.full(Nt,0.)
            r.AuxinSource= [int(0) for i in range(Nt)]
            r.AuxinLost= np.full(Nt,0.)
            r.manualAddAux= np.full(Nt,0.)
            r.C_Auxin= np.full(Nt,0.)
            r.C_AuxinOut= np.full(Nt,0.)
            r.JAuxin_ST2= np.full(Nt,0.)
            r.Delta_JA_ST= np.full(Nt,0.)
        
        Q_ST    = np.array(r.Q_out[0:Nt])
        Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
        Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
        Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
        Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])
        
        
        C_ST    = np.array(r.C_ST)
        Q_Par   = np.full(len(C_ST),0.)#r.Q_out[(Nt*8):(Nt*9)])
        Fl      = np.array(r.Fl)
        volST   = np.array(r.vol_ST)
        volMeso   = np.array(r.vol_Meso)
        C_Par   = Q_Par/volST
        C_meso  = Q_meso/volMeso
        Q_in   += sum(np.array(r.AgPhl)*dt)
        Q_out   = Q_Rm + Q_Exud + Q_Gr
        Q_STLost = np.array(r.SucSTLost)
        Q_MesoLost = np.array(r.SucMesoLost)
        manualAddST =  np.array(r.manualAddST)
        manualAddMeso =  np.array(r.manualAddMeso)
        
        error   = sum(Q_ST +Q_meso + Q_out + Q_STLost + Q_MesoLost)- Q_in - sum(Q_ST_init)  - sum(Q_meso_init) - sum(manualAddST) -sum(manualAddMeso)
        
        Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
        Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])
        Q_Exudmax     = np.array(r.Q_out[(Nt*7):(Nt*8)])
        
        Q_ST_i        = Q_ST      - Q_STbu
        Q_Rm_i        = Q_Rm      - Q_Rmbu
        Q_Gr_i        = Q_Gr      - Q_Grbu
        
        Q_Exud_i      = Q_Exud    - Q_Exudbu
        Q_meso_i      = Q_meso    - Q_mesobu
        
        Q_Rmmax_i     = Q_Rmmax   - Q_Rmmaxbu
        Q_Grmax_i     = Q_Grmax   - Q_Grmaxbu
        Q_Exudmax_i   = Q_Exudmax - Q_Exudmaxbu
        
        Q_out_i       = Q_Rm_i    + Q_Exud_i      + Q_Gr_i
        Q_outmax_i    = Q_Rmmax_i + Q_Exudmax_i   + Q_Grmax_i
          
        assert np.min(Q_ST) >= 0, "np.min(Q_ST) < 0"
        assert np.min(Q_meso) >= 0, "np.min(Q_meso) < 0"
        assert np.min(Q_Rm_i) >= -1e-5, "np.min(Q_Rm_i) < 0"
        assert np.min(Q_Gr_i) >= -1e-5, "np.min(Q_Gr_i) < 0"
        assert np.min(Q_Exud_i) >= -1e-5, "np.min(Q_Exud_i) < 0"
        assert np.min(Q_Rmmax_i) >= -1e-5, "np.min(Q_Rmmax_i) < 0"
        assert np.min(Q_Grmax_i) >= -1e-5, "np.min(Q_Grmax_i) < 0"
        assert np.min(Q_Exudmax_i) >= -1e-5, "np.min(Q_Exudmax_i) < 0"
        
        
        org_ = r.plant.getOrgans(3, False)[0]
        maintstemNodesId = np.array([org_.getNodeIds()])[0]
        OutAuxin = np.array(r.Q_out[(Nt*8):(Nt*9)])
        Q_Auxin =  np.array(r.Q_out[(Nt*9):(Nt*10)])
        Q_Auxin_stem  = Q_Auxin[maintstemNodesId]
        Q_Auxin_other  = np.delete(Q_Auxin,maintstemNodesId)
        InAuxin  += sum(np.array(r.AuxinSource)) * r.auxin_P * dt
        #Q_AuxinInit  = np.array(r.Q_init[(Nt*9):(Nt*10)])
        AuxinDecap   = np.array(r.AuxinLost)
        manualAddAux =  np.array(r.manualAddAux)
        errorAuxin   = sum(Q_Auxin_stem) + sum(Q_Auxin_other) - sum(Q_AuxinInit)- InAuxin +sum(OutAuxin) + sum(AuxinDecap) - sum(manualAddAux)
        #print("sum(Q_Auxin_stem)",sum(Q_Auxin_stem) , sum(Q_Auxin_other) , sum(Q_AuxinInit), InAuxin ,sum(OutAuxin) , sum(AuxinDecap) , sum(manualAddAux))
        #print("errorAuxin",errorAuxin)
        C_Auxin      = np.array(r.C_Auxin)
        AuxinSource  = np.array(r.AuxinSource)
        C_AuxinOut   = np.array(r.C_AuxinOut)
        JAuxin_ST2   = np.array(r.JAuxin_ST2)
        Delta_JA_ST  = np.array(r.Delta_JA_ST)
        """ATT
        
        # C_ST_ = C_ST[1:]
        
        """

        mainStemAux =  np.array(C_Auxin[maintstemNodesId])
        mainStemAux_mean = np.mean(mainStemAux[1:])
        
        if verbosebase :
            print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h", int(((simDuration%1)%24)*(24*60)),"mn", round(r.Qlight *1e6),"mumol m-2 s-1",r.burnInTime ,doDecapitation)
            if Q_in >0.:
                print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,Q_in, 1.)))
            #print("Error in growth:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorGri, relErrorGri))
            print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
            print("water fluxes (cm3/day):\n\ttrans {:5.2e}\tminExud {:5.2e}\tmaxExud {:5.2e}".format(sum(fluxesSoil.values()), min(fluxesSoil.values()), max(fluxesSoil.values())))
            if Q_in >0.:
                print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
                print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
                print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
                print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud)))
                print('init\tST  {:.2e}\tmeso   {:.2e}'.format(sum(Q_ST_init), sum(Q_meso_init)))
                print("aggregated sink satisfaction at last time step (%) :\n\ttot  {:5.1f}\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(
                    sum(Q_out_i)/sum(Q_outmax_i)*100,sum(Q_Rm_i)/sum(Q_Rmmax_i)*100, 
                     div0f(sum(Q_Gr_i), sum(Q_Grmax_i), 1.)*100,div0f(sum(Q_Exud_i),sum(Q_Exudmax_i), 1.)*100))
                print("aggregated sink repartition at last time step (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
                     sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
                print("aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
                     sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
                print("aggregated sink repartition for max (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rmmax_i)/sum(Q_outmax_i)*100, 
                     sum(Q_Grmax_i)/sum(Q_outmax_i)*100,sum(Q_Exudmax_i)/sum(Q_outmax_i)*100))
                print("abs val for max :\n\tRm   {:5.5f}\tGr   {:5.5f}\tExud {:5.5f}".format(sum(Q_Rmmax_i), 
                     sum(Q_Grmax_i),sum(Q_Exudmax_i)))
                print("\tQ_Par {:5.5f}, C_Par {:5.5f}".format(sum(Q_Par), np.mean(C_Par)))
                print("amount Suc (mmol):\n\tAn {:5.2e}\tGr {:5.2e}\tRGr {:5.2e}\n\tRm {:5.2e}\tExud {:5.2e}".format(AnSum, sum(Q_Gr)*r.Gr_Y,sum(Q_Gr)*(1-r.Gr_Y), sum(Q_Rm), sum(Q_Exud))) 
            print("Auxin (mmol):\n\tstem {:5.2e}\totherPl {:5.2e}\tInit {:5.2e}".format(sum(Q_Auxin_stem), sum(Q_Auxin_other), sum(Q_AuxinInit) )) 
            print("\tIn {:5.2e}\tOut {:5.2e}\tAuxCut {:5.2e}".format( InAuxin,sum(OutAuxin), sum(AuxinDecap))) 
            print("Error in Aux_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorAuxin, div0f(errorAuxin,sum(Q_Auxin_stem), 1.)))
            print("C_Auxin (mmol ml-1):\n\tmean {:5.2e}\tmin {:5.2e}\tmax  {:5.2e}".format(mainStemAux_mean, min(mainStemAux[1:]), max(mainStemAux[1:])))
            # print(C_Auxin)
            # print(AuxinSource, max(AuxinSource))
            # print("JAuxin_ST2",JAuxin_ST2)
            # print("Delta_JA_ST",Delta_JA_ST)
            # print("OutAuxin",OutAuxin)
            # print("InAuxin",InAuxin)
            # print("Q_Auxin",Q_Auxin)
            # print("Q_AuxinInit",Q_AuxinInit)
            
        # write_file_array("RhatFhat", RhatFhat)
        # write_file_array("a_STs", a_STs)
        # write_file_array("JW_ST", JW_ST)#cm3/d
        # write_file_array("JS_ST", r.JS_ST)#cm3/d
        # write_file_array("length_ST", length_ST)
        
        ###
        #4UQ
        #
        ###
        
        orgsF = r.plant.getOrgans(3, False)
        orgs = r.plant.getOrgans(3, True)
        toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
        orgs = np.array(orgs)[toKeep]
        #idnodeOfOrg = [org.getNodeIds()[(org.getId()!=1):] for org in orgs]
        #idnodeOfOrg = np.array([item for sublist in idnodeOfOrg for item in sublist])
        # activePhloemv_ = [np.full(len(org.getNodeIds()[(org.getId()!=1):]),org.activePhloem) for org in orgs]
        # activePhloemv = np.array([item for sublist in activePhloemv_ for item in sublist])
        # activePhloemv = [x for _,x in sorted(zip(idnodeOfOrg,activePhloemv))]
        # activeAuxinv_ = [np.full(len(org.getNodeIds()[(org.getId()!=1):]),org.activeAuxin) for org in orgs]
        # activeAuxinv = np.array([item for sublist in activeAuxinv_ for item in sublist])
        # activeAuxinv = [x for _,x in sorted(zip(idnodeOfOrg,activeAuxinv))]
        
        temp_time = time.time()
        #print("num ",thread,temp_time , start_time, "has started. This sampling run took %5.4f seconds." % (temp_time - start_time),'simDuration:',simDuration )
        
        dtreal = ( time.time() - dt_lastWrote)
        dtsimsim = simDuration - dt - dtSIM_lastWrote
        if(not r.burnInTime) and (doDict or ((dtreal >= dt_write) and (dtsimsim >= dtSIM_write) )):
            write_file_array("time_all", np.array([thread, simDuration,simMax,temp_time , start_time,temp_time - start_time]))
            write_file_float("time", simDuration)
            #write_file_array("OutAuxin", OutAuxin)
            #write_file_array("Q_Auxin", Q_Auxin)
            #write_file_array("C_Auxin", C_Auxin)
            write_file_array("C_Auxin_mainStem", mainStemAux[1:])
            #write_file_array("InAuxin", np.array(r.AuxinSource)* r.auxin_P * dt)
            #write_file_array("Delta_JA_ST", Delta_JA_ST)
            #write_file_array("JAuxin_ST2", JAuxin_ST2)


            #write_file_float("doDecapitation", doDecapitation)
            #write_file_float("simMax", simMax)
            #write_file_float("numLNodes",numLNodes)
            #id_orgs = [np.full(org.getNumberOfNodes()-1,org.getId()) for org in orgs]
            #id_orgs = np.array([item for sublist in id_orgs for item in sublist])
            #write_file_array("id_orgsPerNodes", id_orgs)

            #write_file_array("idnodeOfOrg", idnodeOfOrg)
            
            
           
            lengthth_org = np.array([org.getLength(False) for org in orgs])
            length_org = np.array([org.getLength(True) for org in orgs])
            parentOrgId = np.array([org.getParent().getId() for org in orgs])
            ownOrgId = np.array([org.getId() for org in orgs])
            distFromParentBase = np.array([org.getParent().getLength(org.parentNI) for org in orgs])
            #write_file_array("length_org", length_org)
            write_file_array("lengthth_org", lengthth_org)
            #write_file_array("parentOrgId", parentOrgId)
            write_file_array("ownOrgId", ownOrgId)
            write_file_array("distFromParentBase", distFromParentBase)

            cstPerOrg = [np.mean(C_ST[org.getNodeIds()]) for org in orgs]#[1:]
            write_file_array("cstPerOrg", cstPerOrg)#only take shoot C
            #Q_GrmaxPerOrg = [np.mean(Q_Grmax[org.getNodeIds()[1:]]) for org in orgs]
            #write_file_array("Q_GrmaxPerOrg", Q_GrmaxPerOrg)
            #Q_GrPerOrg = [np.mean(Q_Gr[org.getNodeIds()[1:]]) for org in orgs]
            #write_file_array("Q_GrPerOrg", Q_GrPerOrg)
            ot_orgs = np.array([org.organType() for org in orgs])
            st_orgs = np.array([org.getParameter("subType") for org in orgs])
            write_file_array("ot_orgsUQ", ot_orgs)
            write_file_array("st_orgsUQ", st_orgs)

            #activated = np.array([org.activePhloem for org in orgs]) 
            #write_file_array("activatedSUC", activated)
            agesStems = np.array([org.getAge() for org in orgs]) 
            write_file_array("agesStems", agesStems)
            budStage = np.array([org.budStage for org in orgs]) 
            write_file_array("budStage", budStage)
            BerthFact = np.array([org.BerthFact for org in orgs]) 
            write_file_array("BerthFact", BerthFact)
            budStageChange = np.array([org.budStageChange for org in orgs]) 
            write_file_array("bSChange", budStageChange)
            if len(orgs) != len(budStage):
                print(budStage, len(budStage))
                print(ot_orgs,len(ot_orgs))
                print(len(orgs))
                raise Exception
            #activatedaux = np.array([org.activeAuxin for org in orgs]) 
            #write_file_array("activatedAUX", activatedaux)
            
            #ststdistFromParentBase = distFromParentBase[(ot_orgs == 3)* (st_orgs == 2)]
            #ststlactivatedaux = activatedaux[(ot_orgs == 3)* (st_orgs == 2) ]
            #ststactivatedSUC = activated[(ot_orgs == 3)* (st_orgs == 2) ]
            #ststlengthth_org = lengthth_org[(ot_orgs == 3)* (st_orgs == 2) ]
            #ststlengthth_org2 = [x for _,x in sorted(zip(ststdistFromParentBase,ststlengthth_org))]
            #ststlengthth_org[ststdistFromParentBase]
            #mainStemLen = np.array([org_.getLength(nn) for nn in range(org_.getNumberOfNodes())])
            #isbase = np.array([lenlen in ststdistFromParentBase for lenlen in mainStemLen])
            if False:#thread == 100:
                print("Q_ST_stem", len(C_ST))
                # print(sum(mainStemAux[1:] <= r.auxin_threshold))
                print(C_ST[maintstemNodesId])
                print(C_ST[maintstemNodesId] >= r.CSTthreshold)
                
                # print(Q_Auxin_stem)
                # print(manualAddST)
                # print("mainStemAux", len(mainStemAux))
                # print(mainStemAux)
                # print("manualAddAux")
                # print((manualAddAux/volST)[maintstemNodesId])
                # print(mainStemAux <= r.auxin_threshold)
                # print(r.auxin_threshold, simDuration)
                print("mainStemLen")
                print(mainStemLen)
                # print("aux\t",mainStemAux[isbase])
                print("suc\t",C_ST[maintstemNodesId][isbase])
                print("dist\t",mainStemLen[isbase])
                # print("len\t",ststdistFromParentBase)
                print("len\t",ststlengthth_org)
                print("act\t",ststactivatedSUC)#ststlactivatedaux)
                print(maintstemNodesId)
                print("done")
                print()
                if ( (sum(C_ST[maintstemNodesId][isbase] >= r.CSTthreshold ) > sum(ststactivatedSUC))  ):
                    print("sum(C_ST[maintstemNodesId][isbase] >= r.CSTthreshold ) > sum(ststactivatedSUC)")
                    oo += 1
                    #r.doTroubleshooting = True
                    if (oo > 1):
                        raise Exception
            if not doDict:
                dt_lastWrote = time.time()
                dtSIM_lastWrote = simDuration - dt
                #if (sum(mainStemAux[1:] <= r.auxin_threshold)) and not changedSimMax:
                 #   print("(sum(mainStemAux[1:] <= r.auxin_threshold)) and not changedSimMax")
                  #  raise Exception
                #if sum((manualAddAux/volST)[maintstemNodesId]) > 0.:
                 #   print((manualAddAux/volST)[maintstemNodesId])
                    #raise Exception
            
        ###
        #4UQ
        #
        ###
        
        if doVTP:
            ana = pb.SegmentAnalyser(r.plant.mappedSegments())
            
            #print(C_ST)
            #raise Exception
            cutoff = 1e-15 #is get value too small, makes paraview crash
            Delta_JA_ST_p = Delta_JA_ST
            Delta_JA_ST_p[abs(Delta_JA_ST_p) < cutoff] = 0
            C_AuxinOut_p = C_AuxinOut
            C_AuxinOut_p[abs(C_AuxinOut_p) < cutoff] = 0
            JAuxin_ST2_p = JAuxin_ST2
            JAuxin_ST2_p[abs(JAuxin_ST2_p) < cutoff] = 0
            C_Auxin_p = np.full(len(C_Auxin),0.)
            C_Auxin_p[maintstemNodesId] =  np.array(C_Auxin[maintstemNodesId]) 
            C_Auxin_p[abs(C_Auxin_p) < cutoff] = 0
            
            Ag4Phloem_p = np.array(r.Ag4Phloem)
            Ag4Phloem_p[abs(Ag4Phloem_p) < cutoff] = 0
            
            C_ST_p = C_ST
            C_ST_p[abs(C_ST_p) < cutoff] = 0
            fluxes_p = fluxes
            fluxes_p[abs(fluxes_p) < cutoff] = 0
            Q_Exud_i_p = Q_Exud_i
            Q_Exud_i_p[abs(Q_Exud_i_p) < cutoff] = 0
            Q_Rm_i_p = Q_Rm_i
            Q_Rm_i_p[abs(Q_Rm_i_p) < cutoff] = 0
            Q_Gr_i_p = Q_Gr_i
            Q_Gr_i_p[abs(Q_Gr_i_p) < cutoff] = 0
            
            Q_Exudmax_i_p = Q_Exudmax_i
            Q_Exudmax_i_p[abs(Q_Exudmax_i_p) < cutoff] = 0
            Q_Rmmax_i_p = Q_Rmmax_i
            Q_Rmmax_i_p[abs(Q_Rmmax_i_p) < cutoff] = 0
            Q_Grmax_i_p = Q_Grmax_i
            Q_Grmax_i_p[abs(Q_Grmax_i_p) < cutoff] = 0
            
            
            C_Exud_i_p = Q_Exud_i/volST
            C_Exud_i_p[abs(C_Exud_i_p ) < cutoff] = 0
            C_Rm_i_p = Q_Rm_i/volST
            C_Rm_i_p[abs(C_Rm_i_p) < cutoff] = 0
            C_Gr_i_p = Q_Gr_i/volST
            C_Gr_i_p[abs(C_Gr_i_p) < cutoff] = 0
            
            C_Exudmax_i_p = Q_Exudmax_i/volST
            C_Exudmax_i_p[abs(C_Exudmax_i_p) < cutoff] = 0
            C_Rmmax_i_p = Q_Rmmax_i/volST
            C_Rmmax_i_p[abs(C_Rmmax_i_p) < cutoff] = 0
            C_Grmax_i_p = Q_Grmax_i/volST
            C_Grmax_i_p[abs(C_Grmax_i_p) < cutoff] = 0
            
            psiXyl_p = np.array(r.psiXyl)
            psiXyl_p[abs(psiXyl_p) < cutoff] = 0
            
            
            ana.addData("Ag4Phloem_p", Ag4Phloem_p)
            ana.addData("JAuxin_ST2", JAuxin_ST2_p)
            ana.addData("C_AuxinOut", C_AuxinOut_p)
            ana.addData("Delta_JA_ST", Delta_JA_ST_p)
            ana.addData("AuxinSource", AuxinSource)
            ana.addData("isRootTip", np.array(r.isRootTip))
            ana.addData("C_Auxin", C_Auxin_p)
            
            #ana.addData("activePhloem",activePhloemv) 
            #ana.addData("activeAuxin",activeAuxinv)  
            ana.addData("budStage",budStage)             
            
            ana.addData("CST", C_ST_p)
            #do as konz or div per vol or surf?
            #ana.addData("Q_Exud", Q_Exud)  # cut off for vizualisation
            ana.addData("fluxes", fluxes_p)
            ana.addData("Fpsi", np.array(r.Fpsi))
            
            ana.addData("QExud", Q_Exud_i_p)  # cut off for vizualisation
            ana.addData("QRm", Q_Rm_i_p)  # cut off for vizualisation
            ana.addData("QGr", Q_Gr_i_p)  # cut off for vizualisation
            ana.addData("QExudmax", Q_Exudmax_i_p)  # cut off for vizualisation
            ana.addData("QRmmax", Q_Rmmax_i_p)  # cut off for vizualisation
            ana.addData("QGrmax", Q_Grmax_i_p)  # cut off for vizualisation
            
            ana.addData("CExud", C_Exud_i_p)  # cut off for vizualisation
            ana.addData("CRm", C_Rm_i_p)  # cut off for vizualisation
            ana.addData("CGr", C_Gr_i_p)  # cut off for vizualisation
            ana.addData("CExudmax", C_Exudmax_i_p)  # cut off for vizualisation
            ana.addData("CRmmax", C_Rmmax_i_p)  # cut off for vizualisation
            ana.addData("CGrmax", C_Grmax_i_p)  # cut off for vizualisation
            
            ana.addData("psi_Xyl",psiXyl_p)
            ana.write("results"+directoryN+"plotplant_"+dir4allResults+"_"+ str(ö) +".vtp", 
                      ["CST", "fluxes","psi_Xyl",
  "C_Auxin","AuxinSource","Delta_JA_ST","C_AuxinOut","JAuxin_ST2",
     "activePhloem","activeAuxin","Agn4Phloem","isRootTip",
                                "QExud", "QGr", "QRm",
                                "CExud", "CGr", "CRm",
                                "QExudmax", "QGrmax", "QRmmax",
                                "CExudmax", "CGrmax", "CRmmax",
                                "organType", "subType", "Fpsi"]) 
        
        
            fluxes_p = np.insert(fluxes_p,0,0)# "[sucrose]",
            
            #need to adapt plot_plant for this to work
            if(False):#(simDuration > 0):#if(simDuration > (simMax -1)):
                vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","xylem pressure (cm)"],
                                    vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, psiXyl_p], 
                                    filename = "results"+ directoryN +"plotplant_psi_"+dir4allResults+"_"+ str(ö),
                                    range_ = [300,1450])
                vp.plot_plant(r.plant,p_name = [ "fluxes","Exud","Gr","Rm","sucrose concentration (mmol/cm3)"],
                                    vals =[ fluxes_p, Q_Exud_i_p, Q_Gr_i_p, Q_Rm_i_p, C_ST_p], 
                                    filename = "results"+ directoryN +"plotplant_suc_"+dir4allResults+"_"+ str(ö),
                                    range_ = [0,3])   

        if abs(error) > 1e-3:
            print("suc error too high")
            raise Exception    
        if abs(div0f(errorAuxin,InAuxin, 1.)) > 1e-3:
            print("auxin error too high")
            raise Exception    
        if min(C_ST) < 0.0:
            print("min(C_ST) < 0.0", min(C_ST),np.mean(C_ST),max(C_ST))
            raise Exception
            
        if (max(AuxinSource) ==0) and doDecapitation and doStartPM:
            print("no auxin source")
            raise Exception  
            
        if (min(r.Ev) < 0) or (min(r.Jw) < 0) or (abs(errLeuning) > 1e-3) or (min(fluxes_leaves) < 0):
            write_file_array("psiXyl", r.psiXyl)
            write_file_array("trans", r.Ev)
            write_file_array("transrate",r.Jw)
            write_file_array("fluxes", fluxes)#cm3 day-1
            write_file_array("fluxes_leaves", fluxes_leaves)
            r.doTroubleshooting = True
            r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
            verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
            print(toKil)
            print(np.array(r.psiXyl)[toKil])
            print(np.where(np.array(r.psiXyl)<-10000))
            print(np.where(fluxes_leaves<-0))
            print("leaf gaines water", min(r.Ev),min(r.Jw), min(fluxes_leaves))
            
            org_ = r.plant.getOrgans(3, False)[0] #get first (==main) stem
            print("org_.active",org_.active,LL, org._getLength())
            print(np.array([org_.getNodeId(org_.getNumberOfNodes()-2),org_.getNodeId(org_.getNumberOfNodes()-1)]))
            raise Exception
            
        ö +=1
        
        
        if not r.burnInTime:
            r.plant.simulate(dt, False)
            
            
        Ntbu = Nt
        Nt = len(r.plant.nodes)
        
            


        
        Q_Rmbu       =   np.concatenate((Q_Rm, np.full(Nt - Ntbu, 0.)))
        Q_Grbu       =   np.concatenate((Q_Gr, np.full(Nt - Ntbu, 0.))) 
        Q_Exudbu     =   np.concatenate((Q_Exud, np.full(Nt - Ntbu, 0.))) 
        
        Q_Rmmaxbu    =   np.concatenate((Q_Rmmax, np.full(Nt - Ntbu, 0.)))
        Q_Grmaxbu    =   np.concatenate((Q_Grmax, np.full(Nt - Ntbu, 0.))) 
        Q_Exudmaxbu  =   np.concatenate((Q_Exudmax, np.full(Nt - Ntbu, 0.))) 
        
        Q_STbu       =   np.concatenate((Q_ST, np.full(Nt - Ntbu, 0.)))
        Q_mesobu     =   np.concatenate((Q_meso, np.full(Nt - Ntbu, 0.)))
        
        
        if r.burnInTime :
            burninDuration += dt
            org_ = r.plant.getOrgans(3, False)[0]
            maintstemNodesId = np.array([org_.getNodeIds()])
            mainStemAux =  np.array(C_Auxin[maintstemNodesId])[0]
            mainStemAux_std = np.std(mainStemAux[1:])
            mainStemAux_mean = np.mean(mainStemAux[1:])
        else :
            simDuration += dt
            
        # if r.burnInTime and \
        #  ((activeAtThreshold_suc and (min(C_ST)>0.)) or \
        #  (activeAtThreshold_auxin and (mainStemAux_std < 0.01) and (mainStemAux_mean > 1) )):
        if r.burnInTime and (mainStemAux_std < 0.01) and (mainStemAux_mean > 1) : 
            r.burnInTime = False
            
            #r.canStartActivating = True
            #if UseRatiothresholdAux and (r.auxin_threshold ==0):#only set it tyhe first time
             #   r.auxin_threshold =mainStemAux_mean*RatiothresholdAux
            print(r.auxin_init_mean)
            r.auxin_init_mean = mainStemAux_mean
            
            with open('results'+ directoryN+"input"+dir4allResults+ '.txt', 'r+') as f:
                lines = f.read().splitlines()
            lines[0] = lines[0] + ",mainStemAux_mean"
            lines[1] = lines[1] + "," + repr(mainStemAux_mean)
            with open('results'+ directoryN+"input"+dir4allResults+ '.txt', 'w') as f:
                f.write(lines[0])
                f.write("\n")
                f.write(lines[1])
                
            errorBU = error
            Q_outTemp = np.full(len(r.Q_out),0.)
            Q_outTemp[(Ntbu*9):(Ntbu*10)] = Q_Auxin
            Q_outTemp[0:Ntbu] = Q_ST   
            Q_outTemp[Ntbu:(Ntbu*2)] = Q_meso 
            r.Q_out = Q_outTemp
            r.Q_init = Q_outTemp
            InAuxin =0
            Q_in = 0
            Q_Rmbu       =   np.full(Nt , 0.)
            Q_Grbu       =   np.full(Nt , 0.)
            Q_Exudbu     =   np.full(Nt , 0.) 
            Q_Rmmaxbu    =   np.full(Nt , 0.)
            Q_Grmaxbu    =   np.full(Nt , 0.) 
            Q_Exudmaxbu  =   np.full(Nt , 0.) 
            Q_STbu       =   np.full(Nt , 0.)
            Q_mesobu     =   np.full(Nt , 0.)
            #test new balance 
            Q_ST    = np.array(r.Q_out[0:Ntbu])
            Q_meso  = np.array(r.Q_out[Ntbu:(Ntbu*2)])
            Q_Rm    = np.array(r.Q_out[(Ntbu*2):(Ntbu*3)])
            Q_Exud  = np.array(r.Q_out[(Ntbu*3):(Ntbu*4)])
            Q_Gr    = np.array(r.Q_out[(Ntbu*4):(Ntbu*5)])
            #Q_in   += sum(np.array(r.AgPhl)*dt)
            Q_out   = Q_Rm + Q_Exud + Q_Gr
            Q_STLost = np.array(r.SucSTLost)
            Q_MesoLost = np.array(r.SucMesoLost)
            if (len(r.Q_init) >0 ) and (len(Q_ST_init) ==0):
                Q_ST_init = np.array(r.Q_init[0:Ntbu])
                Q_meso_init = np.array(r.Q_init[Ntbu:(Ntbu*2)])
                Q_AuxinInit  = np.array(r.Q_init[(Ntbu*9):(Ntbu*10)])
                
            error   = sum(Q_ST + Q_meso + Q_out + Q_STLost + Q_MesoLost)- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
            if abs(errorBU) < abs(error):
                print("abs(errorBU) < abs(error)", errorBU , error)
                print(len(r.Q_init), sum(Q_ST), sum(Q_meso), sum(Q_out), sum(Q_STLost), sum(Q_MesoLost), Q_in , sum(Q_ST_init)  , sum(Q_meso_init))
                raise Exception
            errorAuxinBU = errorAuxin 
            OutAuxin = np.array(r.Q_out[(Ntbu*8):(Ntbu*9)])
            Q_Auxin  = np.array(r.Q_out[(Ntbu*9):(Ntbu*10)])
            #InAuxin  += sum(np.array(r.AuxinSource)) * r.auxin_P * dt
            
            AuxinDecap   = np.array(r.AuxinLost)
            errorAuxin   = sum(Q_Auxin) - sum(Q_AuxinInit)- InAuxin +sum(OutAuxin) + sum(AuxinDecap)
            if abs(errorAuxinBU) < abs(errorAuxin):
                print("errorAuxinBU < errorAuxin", errorAuxinBU , errorAuxin)
                print(sum(Q_Auxin) , sum(Q_AuxinInit), InAuxin ,sum(OutAuxin) , sum(AuxinDecap))
        
            
        if doDict and (( time.time() - dt_lastWrote) >= dt_write) and (( simDuration - dt - dtSIM_lastWrote ) >= dtSIM_write)  :
            #print("write",time.time())
            for keykey in outputsDict_array.keys():
                name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.txt'
                datas = outputsDict_array[keykey]
                with open(name2, 'a') as log:
                    for data in datas:
                        log.write(','.join([num for num in map(str, data)])  +'\n')
            for keykey in outputsDict_float.keys():
                name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.txt'
                datas = outputsDict_float[keykey]
                with open(name2, 'a') as log:
                    for data in datas:
                        log.write(repr( data)  +'\n')
            outputsDict_array.clear()
            outputsDict_float.clear()
            dt_lastWrote = time.time()
            dtSIM_lastWrote = simDuration - dt
    if doDict:        
        for keykey in outputsDict_array.keys():
            name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.txt'
            datas = outputsDict_array[keykey]
            with open(name2, 'a') as log:
                for data in datas:
                    log.write(','.join([num for num in map(str, data)])  +'\n')
        for keykey in outputsDict_float.keys():
            name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.txt'
            datas = outputsDict_float[keykey]
            with open(name2, 'a') as log:
                for data in datas:
                    log.write(repr( data)  +'\n')
        outputsDict_array.clear()
        outputsDict_float.clear()
        
    dtreal = ( time.time() - dt_lastWrote)
    dtsimsim = simDuration - dt - dtSIM_lastWrote
    if (not doDict) and ("C_Auxin_mainStem" in locals()):
        write_file_array("time_all", np.array([thread, simDuration,simMax,temp_time , start_time,temp_time - start_time]))
        write_file_float("time", simDuration)
        write_file_array("C_Auxin_mainStem", mainStemAux[1:])
        write_file_array("agesStems", agesStems)
        #write_file_float("numLNodes",numLNodes)
        write_file_array("lengthth_org", lengthth_org)
        write_file_array("ownOrgId", ownOrgId)
        write_file_array("distFromParentBase", distFromParentBase)

        write_file_array("cstPerOrg", cstPerOrg)#only take shoot C
        write_file_array("ot_orgsUQ", ot_orgs)
        write_file_array("st_orgsUQ", st_orgs)
        write_file_array("budStage", budStage)
        write_file_array("bSChange", budStageChange)
        #write_file_array("activatedSUC", activated)
        #write_file_array("activatedAUX", activatedaux)
        
    #ana = pb.SegmentAnalyser(r.plant.mappedSegments())
    #ana.write("results"+directoryN+"plotplant_"+dir4allResults +".vtp", 
     #         ["organType", "subType"]) 
    # organType = r.plant.organTypes
    # subType = r.plant.subTypes
    # vp.plot_plant(r.plant,p_name = ["organType", "subType"],
    #                     vals =[ organType, subType], 
    #                     filename = "results"+ directoryN +"plotplant_psi_"+dir4allResults,
    #                     range_ = [0,10])  
    delete_file("running")
    print("finished", thread, time.time() - start_time)
    #os._exit(os.EX_OK)
    del r
    del pl
    gc.collect()
    return 0


if __name__ == '__main__':
    start_time_ = time.time()
    main_dir=os.environ['PWD']#dir of the file
    directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
    results_dir = main_dir +"/results"+directoryN

    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    else:
        test = os.listdir(results_dir)
        for item in test:
            try:
                os.remove(results_dir+item)
            except:
                pass
    
    runSim(
        directoryN_ = directoryN, doVTP = False,verbosebase = True, 
        PRate_ = 6.8e-3,thresholdAux = 0, RatiothresholdAux =0.46,
        UseRatiothresholdAux = True,
        Qmax_ = 450e-6, thresholdSuc = 1.68,  
        GrRatio = 10, maxLBud = 1., budGR = 0.1,L_dead_threshold=2.,
        nodeD = 3, thread = 0,
        testTime=0.5, dtBefore = 1/24, dtAfter= 1/24,
        start_time = start_time_,
        dt_write = 0, dtSIM_write = 10/(60*24), 
        doPrint = True, doDict = False,auxin_D=0.
    )
    end_time_ = time.time()
    print(end_time_ - start_time_ )
    sys.exit(0)