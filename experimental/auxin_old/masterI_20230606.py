""" 
Q:
- if we make kx age dependent better
- Y do we get C_ST > 0 during burn-in time without sun?
- make sure lats have >=2 nodes also when sleeping buds (at least when reach max sleeping size).
- add/check decapita
- why runtime so high?
- remove added C and AAI for added segment of laterals? (or just not add C? because can add fals acitvation of carbon? at least for the first segment?)
- set k_S_ST, C_targ, Q_S_ST init, init Cmeso
- stop respiration if it leads to C_ST below 0
- see where that comes from: corrupted size vs. prev_size in the new version
- check that od and new version give the same restuls
- StopLoss: do only for auxin
- with so much exudation, we should have the roots growing. also why is growth_th changing for the static roots.
- recalibrate for correct carbon-dependent growth:
    - highlight when grown: get the correct growth
    - k_S_ST, C_me => correct growth at night when growth
    - get concentration pick after decapitation
    - reserve of S_ST at the beginning (change of density or An?)
- the computation is slow again. but is it my computer or is it the setup?
"""
import types
import importlib
import os
import sys
SRC_PATH = "../../src/"
sys.path.append("../.."); sys.path.append(SRC_PATH)

# Create a fake plantbox namespace
plantbox = types.SimpleNamespace()

# Automatically import all folders inside src and attach to plantbox
for name in os.listdir(SRC_PATH):
    folder_path = os.path.join(SRC_PATH, name)
    if os.path.isdir(folder_path) and not name.startswith('__'):
        try:
            module = importlib.import_module(name)
            setattr(plantbox, name, module)
            sys.modules[f'plantbox.{name}'] = module
        except ModuleNotFoundError:
            # skip folders that are not importable as modules
            pass



import sys; 
#directoryN = "/"+sys.argv[0].split('.')[0]+"/"
sys.path.append("../.."); sys.path.append("../../src/python_modules")
CPBdir = "../.."
#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs
import gc
from plantbox.functional.phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import psutil

#import matplotlib.pyplot as plt
import time
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.phloem_flux import PhloemFluxPython 
from plantbox.functional.PlantHydraulicParameters import PlantHydraulicParameters 
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
    coef = int((hours)<=(12/24))
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
    ea = es * RH
    pmean =-100#theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,
                    'es':es,'ea':ea,
                    'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
    #print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c = 1.):    
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
                                     
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ,kr_s,kr_s,kr_s,kr_s ,kr_s,kr_s ],[kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l]], kr_length_=l_kr)
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s,kz_s,kz_s,kz_s,kz_s,kz_s,kz_s ],[kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm

    r.psi_air = p_a #*MPa2hPa #used only with xylem
    
def setKrKx_phloem(r): #inC

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
    kroot = 5e-3
    kr_r0 = kroot
    kr_r1 = kroot
    kr_r2 = kroot 
    kr_r3 = kroot
    l_kr = 0.8 #cm
    
    r.cpb_2_pm.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s, kr_s,kr_s,kr_s,kr_s, kr_s,kr_s],[kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l]] , 
                kr_length = l_kr)
    r.cpb_2_pm.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_l,kz_l, kz_l,kz_l,kz_l,kz_l, kz_l,kz_l],
        [kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l]])
    
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
    r.cpb_2_pm.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s,Across_s_s,Across_s_s],[Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])

def doConditionDefault(rinput, timeSinceDecap_,i, simtime, memoryCondition, nodeD, inputdata):
    return False
def successPointsDefault(rinput, timeSinceDecap_,i, simtime, memoryCondition, nodeD, inputdata):
    return memoryCondition

from pathlib import Path

def killChildren(orgToKill):
    toFill = orgToKill.getNodeIds()
    orgToKill.alive = False
    orgToKill.active = False
    orgToKill.budStage = -1
    if orgToKill.getNumberOfChildren() > 0:
        toFill_ = [killChildren(orgToKill.getChild(ni)) for ni in range(orgToKill.getNumberOfChildren())]
        toFill_ = [item for sublist in toFill_ for item in sublist]
        toFill = toFill + toFill_
    return toFill
          
def runSim(startDate, directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, 
           Qmax_, thresholdSuc,useLength,
           maxLBud , maxLBudDormant,#maxLBudDormant_1,
           L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D = 0.,kss=0.2,kaa=1., 
           fileparam ="UQ_1Leaf" , doCondition = doConditionDefault,
           Klight = 0.05, BerthLim = -1,
          simMax_= -1.,
          leafAsIAASource_ = True, limLenActive_ = 0.9, growthUpToNode = 9, successPoints = successPointsDefault):
    outcondition = 0
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
    
    dir4allResults = "_"+str(thread) #"_"+ strPRate_+ "_"+strThA+ "_"+strThRatioA + "_"+ strQ + "_"+strThS+"_"+strDecap+"_"+str(thread)
    dir4allResults=dir4allResults.replace(".", "o")
    def write_file_array(name, data):
        if doPrint :
            name2 = 'results'+ directoryN+ name+ dir4allResults+ '.csv'
            with open(name2, 'a') as log:
                log.write(','.join([num for num in map(str, data)])  +'\n')
        if doDict:
            if name not in outputsDict_array:
                outputsDict_array[name] = list()
            outputsDict_array[name].append(data)
        

    def write_file_float(name, data):
        if doPrint:
            name2 = 'results'+ directoryN+ name+ dir4allResults+ '.csv'
            with open(name2, 'a') as log:
                log.write(repr( data)  +'\n')
        if doDict:
            if name not in outputsDict_float:
                outputsDict_float[name] = np.array([])
            outputsDict_float[name] = np.append(outputsDict_float[name],data )

    def delete_file(name):
        if doPrint:
            name2 = 'results'+ directoryN+ name+dir4allResults+ '.csv'
            my_file = Path(name2)
            if my_file.is_file():
                os.remove(name2)
            else:
                print("file running missing",thread, name2)
    
    namef = 'results'+ directoryN+ "running"+ dir4allResults+ '.csv'
    my_file = Path(namef)
    if my_file.is_file():
        print("input", np.array([directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, 
           Qmax_, thresholdSuc,
           maxLBud , L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D,Klight ]))
        print(namef)        
        print("runningfile already exists")
        raise Exception
    toAdd = 0
    doPrintbu = doPrint
    doDictbu = doDict
    doDict = False
    doPrint = True
    write_file_float("running", thread)  
    allInputs = np.array(["directoryN_","doVTP", "verbosebase",
           "PRate_", "thresholdAux", "RatiothresholdAux", 
           "Qmax_", "thresholdSuc","useLength",
            "maxLBud[0]" , "maxLBudDormant[0]","maxLBudDormant[1]","maxLBudDormant[2]",
                      "L_dead_threshold" ,
           "nodeD", "thread",  
           "testTime", "dtBefore", "dtAfter", "start_time", "dt_write","dtSIM_write",
           "doPrint" ,"doDict", "auxin_D","kss", "kaa" ,
                       "fileparam"  , 
           "Klight" , "BerthLim" ])
    allInputsData = np.array([directoryN_,doVTP, verbosebase,
           PRate_, thresholdAux, RatiothresholdAux, 
           Qmax_, thresholdSuc,useLength,
           maxLBud[0] , maxLBudDormant[0],maxLBudDormant[1],maxLBudDormant[2],
           L_dead_threshold ,
           nodeD, thread,  
           testTime, dtBefore, dtAfter, start_time, dt_write,dtSIM_write,
           doPrint ,doDict, auxin_D ,kss,kaa, 
           fileparam  , 
           Klight , BerthLim ])
    write_file_array("input", allInputs) 
    write_file_array("input", allInputsData) 
    
    allInputsdict = {}
    allInputsDataL = allInputsData.tolist()
    for key in allInputs.tolist():
        for value in allInputsDataL:
            allInputsdict[key] = value
            allInputsDataL.remove(value)
            break

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
    name = fileparam

    pl.readParameters("./modelparameter/UQ_1LeafRS.xml") #path + name + ".xml")


    #raise Exception
    sdf = pb.SDF_PlantBox(np.inf, np.inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    #pl.activeAtThreshold_auxin = activeAtThreshold_auxin
    #pl.activeAtThreshold = activeAtThreshold_suc
    pl.initialize(verbose = False)#, stochastic = False)
    pl.maxLBud = maxLBud
    pl.maxLBudDormant = maxLBudDormant
    pl.budGR = 1.8
    leafArea = np.array([])
    budlength = np.array([])
    while (len(pl.getOrgans(3, False)) ==0 ) or (pl.getOrgans(3, False)[0].getNumberOfLinkingNodes() < (3)) or (len(leafArea)<3) or (sum(leafArea) < (0.69)):
        pl.simulate(startDate, False) #20
        scalLeaves = pl.getOrgans(4, True)
        leafArea = np.array([org.getLength(True) * org.getParameter("Width_blade") for org in scalLeaves])
        stems = pl.getOrgans(3, True)[1:]
        budlength = np.array([org.getLength(True) for org in stems])
        if __name__ == '__main__':
            print('leafArea', leafArea)
            print('budlength', budlength)
        if((len(pl.getOrgans(3, False)) >0 ) and (__name__ == '__main__')):
            print('numofnodes',pl.getOrgans(3, False)[0].getNumberOfLinkingNodes(), 
            'numkid',pl.getOrgans(3, False)[0].getNumberOfChildren(),
            'length',pl.getOrgans(3, False)[0].getLength(True))
        leafArea[:2]
        simDuration += startDate
    
    
    stems = np.array(pl.getOrgans(3, True))
    toKeep = np.array([org.getParameter("subType") <= 2 for org in stems])
    stems  = stems[toKeep]
    stemlengths = np.array([org.getLength(False) for org in stems])
    if(min(stemlengths) == 0) and (min(maxLBudDormant)>0):
        print("min(stemlengths) == 0)",stemlengths ,maxLBudDormant)
        raise Exception
        
    if simMax_ < 0:
        simMax = 100#simDuration + testTime #if !doDecapitation
    else:
        simMax = simMax_

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

    pl.setSoilGrid(picker)  # maps segment


    """ Parameters phloem and photosynthesis """
    global r
    
    params = PlantHydraulicParameters()  # |\label{l52:hydraulic}|
    #setKrKx_xylem(20.,params)
    params.read_parameters("./modelparameter/pea_toredo")  # |\label{l52:hydraulic_end}|
    
    r = PhloemFluxPython(pl,params, psiXylInit = min(sx),ciInit = weatherInit["cs"]*0.5) #XylemFluxPython(pl)#
    r.initialize()
    setKrKx_phloem(r)
    r.g0 = 8e-2
    r.VcmaxrefChl1 =2
    r.VcmaxrefChl2 =7
    r.a1 = 2
    r.a3 = 2.2
    r.alpha = 0.27
    r.theta = 0.51
    r.k_meso = 1e-3
    r.cpb_2_pm.setKrm2([[2e-5 ]])
    r.cpb_2_pm.setKrm1([[1.3e-2 ]])
    rho_org = [[0.,1.3,1.3,1.3],[0.,1.4,1.4,1.4],[0.,1.5,1.5,1.5,1.5]]#g glucose/gDW
    
    density = 0.17 #g DW/cm3
    g_2_mmolSuc = 2.92 
    #density and rho in mmol suc
    rho_org = np.array([np.array(xi) for xi in rho_org],dtype=object)*density/2 * g_2_mmolSuc #/2 => glucose to sucrose
    r.cpb_2_pm.setRhoSucrose(rho_org)
    grRoot = 1.3
    grRate = [[0.,2.*grRoot,0.7*grRoot,0.5*grRoot,2.*grRoot],[0.,1.,1.,1.],[0.,1.,1.,1.,1.]]
    grRate = np.array([np.array(xi) for xi in grRate],dtype=object)
    
    r.cpb_2_pm.setRmax_st(grRate)
    r.KMfu = 0.2
    #r.exud_k = np.array([2.4e-4])#*10#*(1e-1)
    #r.k_gr = 1#0
    r.sameVolume_meso_st = False
    r.sameVolume_meso_seg = True
    r.withInitVal = True
    r.initValST = 0.3
    if startDate < 10.:
        r.initVal_S_ST = 1e10
    r.initValMeso = 0.#0.2
    r.beta_loading = 10
    r.Vmaxloading = 0.3 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 3.3*1e-3#0.2#
    r.Gr_Y = 1#0.75
    r.CSTimin = 0.25#0.3#    
    r.k_S_ST = 5/25 *1000  #daudet2002 * 1000 
    r.k_S_meso = 0.# 5/25 *1000 
    print("change for the cmeso day/night")
    r.C_targ = r.CSTimin  * 6
    r.C_targ_meso = r.CSTimin 
    #r.surfMeso=0.0025
    #r.cs = weatherInit["cs"]

    r.expression = 6
    r.update_viscosity = False
    r.solver = 1
    r.atol = 1e-10# 1e-14
    r.rtol = 1e-6#1e-10
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 60.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4

    #turn off effect of water limitation
    r.cpb_2_pm.psiMin = -10000000000*(1/0.9806806)

    """ for post processing """
    structSumInit = 0
    orgs_all = r.plant.getOrgans(-1, True)

    for org in orgs_all:
        if org.organType() < 2:
            raise Exception("ot < 3")
        #structSumInit += org.orgVolume(-1,False) * r.rhoSucrose_f(
         #   int(org.getParameter("subType")-1),
          #                 org.organType())

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
    r.thread = thread

    r.auxin_threshold = thresholdAux
    r.auxin_D = auxin_D #e-6#muM /d  #e-6#3e-7
    r.auxin_P = PRate_ #muM /d  #*1e-6# 6e-7
    r.auxin_alpha =0.00239336#0.00270427# 0.002727303 #0.002704272#0.00245669#0.02 #5#e-6# 2.4e-3
    r.initValAuxin = 0 #0.2
    r.deleteAtRootTip = True
    r.plant.useCWGr = useCWGr 
    InAuxin = 0
    #r.canStartActivating = False
    r.cpb_2_pm.CSTthreshold = thresholdSuc 
    r.StopLoss = True #// when add a new segment, create an initial concentraiton to avoid empty puls
    
    r.cpb_2_pm.L_dead_threshold = L_dead_threshold
    orgs_stems = r.plant.getOrgans(3, True)

    r.computeBerth = lambda ss, aa, org:  ((aa+kaa)/(ss+kss))#*(1-(0.15/(ss+kss)))*2 
    #0.5/(1+np.exp(-kaa*(aa-0.4/2)))+0.5/(1+np.exp(kss*(ss+1)))
    #10/(1+np.exp(-kaa*(aa-2.4/2))+np.exp(kss*(ss+1)))
    #(((aa_*100)**kaa)/100**kaa/(ss_*10+kss))*(1-(0.15/(ss_*10+kss)))*2 
    #teset
    if __name__ == '__main__':
        print(r.computeBerth(0.1,2.49,  r.plant.getOrgans(3, False)[0]))
    
    r.cpb_2_pm.BerthLim = BerthLim
    r.useLength = useLength
    r.leafAsIAASource = leafAsIAASource_ #only young expending leaf production is relevant
    r.cpb_2_pm.limLenActive = limLenActive_

    
    
    
    burninDuration = simDuration
    r.burnInTime = True#activeAtThreshold_auxin

    changedSimMax = False
    oo = 0

    stem = r.plant.getOrgans(3, False)[0]
    leaves = np.array([stem.getChild(np) for np in range(stem.getNumberOfChildren()) if stem.getChild(np).organType() == pb.leaf])
    roots = r.plant.getOrgans(2, True)
    lln = stem.getLlocalId_linking_nodes()
    lln = np.diff([stem.getLength(lln_i) for lln_i in lln]) + np.array(stem.epsilonDxPerPhyto) 
    ll = np.array([ll.getLength(False) for ll in leaves])

    lr = np.array([ll_.getLength(False) for ll_ in roots])
    idRoots = np.array([ll_.getId() for ll_ in roots])
    lr = dict(zip(idRoots, lr))
    grr_th = np.array([np.nan])
    while  simDuration < simMax:#
        temp_time = time.time()
        
            
        if r.burnInTime:
            weatherX = weather(burninDuration, Qmax_)
        else:
            weatherX = weather(simDuration+toAdd, Qmax_)
            
            
        # print('solve water')
        # r.solve(sim_time = simDuration, 
                # rsx=sx, 
                 # ea= weatherX["ea"], es= weatherX["es"], 
                 # PAR = weatherX["Qlight"] * (24 * 3600) / 1e4 , # [mol photons m-2 s-1] to [mol photons cm-2 d-1] 
                 # TairC = weatherX["TairC"])  
        # print('ok water')
        
        r.seg_leaves_idx = r.rs.getSegmentIds(4)
        leafBladeSurface = np.array(r.rs.leafBladeSurface)[np.array( r.seg_leaves_idx)]
        An = np.zeros(len(leafBladeSurface))
        An[leafBladeSurface > 0.] = 15e-6 * weatherX["Qlight"]/(1000e-6 ) /2.
        r.An = An
        assert((weatherX["Qlight"]/(1000e-6 )) <= 1.)
        """ dumux """    
        r.getAg4Phloem() # 0.88  mmol Suc d-1
        # growth shoot, theory = 0.88*0.2 = 0.177  mmol Suc d-1
        # real growth use: 
        #       stem: 0.5 * 2* np.pi * (0.5 ** 2) * 1.4 * 0.17 = 0.187  mmol Suc d-1
        #       leaf: 0.5 * 4 * 0.05 * 6.8 * 1.5 * 0.17
        
        
        AnSum += np.sum(r.Ag4Phloem)*dt
        startphloem= simDuration
        endphloem = startphloem + dt
        stepphloem = 1
        
        # save_stdout = sys.stdout
        # sys.stdout = open('trash', 'w')
        verbose_phloem = True
        
        filename = "results"+ directoryN +"inPM_"+str(ö)+dir4allResults+ '.txt'
        
        
        r.useStemTip = True
        doStartPM = True
        r.doTroubleshooting = False
        r.stopAt = -1 
        if doStartPM:
            try:
                print('solve sucrose')
                r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)                
                print('ok sucrose')            
            except:
                print("error while running r.startPM",thread, simDuration)
                r.doTroubleshooting = True
                #r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
                raise Exception("error while running r.startPM")
        if ((r.stopAt >= 0) and (r.stopAt < 13)) or (not doStartPM): 
            r.Q_out = np.full(Nt*10,0.)
            r.C_ST = np.full(Nt,0.)
            r.vol_ST = np.full(Nt,1.)
            r.Fl = np.full(Nt,0.)
            r.vol_Meso = np.full(Nt,1.)
            r.SucSTLost= np.full(Nt,0.)
            r.SucMesoLost= np.full(Nt,0.)
            r.manualAddST= np.full(Nt,0.)
            r.manualAddMeso= np.full(Nt,0.)
            r.cpb_2_pm.AuxinSource= [int(0) for i in range(Nt)]
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
        Q_S_ST  = np.array(r.Q_out[(Nt*7):(Nt*8)])
        
        
        C_ST    = np.array(r.C_ST)
        Q_Par   = np.full(len(C_ST),0.)#r.Q_out[(Nt*8):(Nt*9)])
        Fl      = np.array(r.Fl)
        volST   = np.array(r.vol_ST)
        volMeso   = np.array(r.vol_Meso)
        C_Par   = Q_Par/volST
        C_S_ST   = Q_S_ST/volST
        C_meso  = Q_meso/volMeso
        Q_in   += sum(np.array(r.AgPhl)*dt)
        Q_out   = Q_Rm + Q_Exud + Q_Gr
        Q_STLost = np.array(r.SucSTLost)
        Q_MesoLost = np.array(r.SucMesoLost)
        manualAddST =  np.array(r.manualAddST)
        manualAddMeso =  np.array(r.manualAddMeso)
            
        if (r.withInitVal) and (len(Q_ST_init) == 0) :
            Q_ST_init = r.initValST * volST
            Q_S_ST_init = r.initVal_S_ST * volST
            Q_meso_init = r.initValMeso * volMeso
            Q_meso_init[volMeso < 0.] = 0. 
        
        # if (len(r.Q_init) >0 ) and (len(Q_ST_init) ==0): # somehow not working
            # Q_ST_init = np.array(r.Q_init[0:Ntbu])
            # Q_meso_init = np.array(r.Q_init[Ntbu:(Ntbu*2)])
            # Q_AuxinInit  = np.array(r.Q_init[(Ntbu*9):(Ntbu*10)])
            
        error   = sum(Q_ST +Q_meso + Q_out + Q_STLost + Q_MesoLost + Q_S_ST)- Q_in - sum(Q_ST_init) - sum(Q_S_ST_init)  - sum(Q_meso_init) - sum(manualAddST) -sum(manualAddMeso)
        Q_labile = Q_ST + Q_meso + Q_S_ST
        Q_Rmmax       = np.array(r.Q_out[(Nt*5):(Nt*6)])
        Q_Grmax       = np.array(r.Q_out[(Nt*6):(Nt*7)])
        Q_GrmaxBis       = r.Q_Grmaxv
        # print(Q_Grmax)
        # print(Q_GrmaxBis)
        # print(dt)
        # print(np.array(Q_Grmax) - np.array(Q_GrmaxBis)*dt)
        #raise Exception
        Q_Exudmax     = Q_Exud.copy()
        
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
        InAuxin  += sum(np.array(r.cpb_2_pm.AuxinSource)) * r.auxin_P * dt
        #Q_AuxinInit  = np.array(r.Q_init[(Nt*9):(Nt*10)])
        AuxinDecap   = np.array(r.AuxinLost)
        manualAddAux =  np.array(r.manualAddAux)
        errorAuxin   = sum(Q_Auxin_stem) + sum(Q_Auxin_other) - sum(Q_AuxinInit)- InAuxin +sum(OutAuxin) + sum(AuxinDecap) - sum(manualAddAux)
        #print('Q_Auxin_stem',Q_Auxin_stem, 'Q_Auxin_other', Q_Auxin_other, 'Q_AuxinInit', Q_AuxinInit, 
        #        'InAuxin', InAuxin, 'OutAuxin', OutAuxin, 'AuxinDecap', AuxinDecap, 'manualAddAux', manualAddAux)
        
        C_Auxin      = np.array(r.C_Auxin)
        AuxinSource  = np.array(r.cpb_2_pm.AuxinSource)
        C_AuxinOut   = np.array(r.C_AuxinOut)
        JAuxin_ST2   = np.array(r.JAuxin_ST2)
        Delta_JA_ST  = np.array(r.Delta_JA_ST)
        """ATT
        
        # C_ST_ = C_ST[1:]
        
        """

        mainStemAux =  np.array(C_Auxin[maintstemNodesId])
        mainStemAux_mean = np.mean(mainStemAux[1:])
        if(ö == 0):
            mainStemAuxBU = np.array(C_Auxin[maintstemNodesId][1:])
            mainStemAuxBU[:] = -1
            
            
            
        if __name__ == '__main__':
            #int(((simDuration%1)%24)*(24*60)),"mn",
            print("\n\n\n\t\tat ", int(np.floor(simDuration)),"d", int((simDuration%1)*24),"h",round(weatherX["Qlight"] *1e6),"mumol m-2 s-1",
                    'burnInTime?',r.burnInTime,'doDecapitation?',doDecapitation)
            if Q_in >0.:
                print("Error in Suc_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(error, div0f(error,Q_in, 1.)))
            #print("Error in growth:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(errorGri, relErrorGri))
            #print("Error in photos:\n\tabs (cm3/day) {:5.2e}".format(errLeuning))
            #print("water fluxes (cm3/day):\n\ttrans {:5.2e}\tmin {:5.2e}\tmax {:5.2e}".format(sum(fluxesSoil.values()), min(fluxesSoil.values()), max(fluxesSoil.values()))) Q_S_ST
            if Q_in >0.:
                print("tot balance",'leafBladeSurface',sum(leafBladeSurface),'an',sum(r.An), sum(r.Ag4Phloem),sum(r.AgPhl),
                            'Q_S_ST', sum(Q_S_ST), 'Q_meso', sum(Q_meso), 'Q_ST', sum(Q_ST), 'Q_Gr', sum(Q_Gr),
                                'Q_Rm', sum(Q_Rm), 'Q_Exud', sum(Q_Exud))
                print("C_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} at {:d} segs \tmax  {:5.2e}".format(np.mean(C_ST), min(C_ST), len(np.where(C_ST == min(C_ST) )[0]), max(C_ST)))        
                print("C_S_ST (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e} \tmax  {:5.2e}".format(np.mean(C_S_ST), min(C_S_ST), max(C_S_ST)))            
                print("C_me (mmol ml-1):\n\tmean {:.2e}\tmin  {:5.2e}\tmax  {:5.2e}".format(np.mean(C_meso), min(C_meso), max(C_meso)))        
                #print('Q_X (mmol Suc): \n\tST   {:.2e}\tmeso {:5.2e}\tin   {:5.2e}'.format(sum(Q_ST), sum(Q_meso), Q_in))
                #print('\tRm   {:.2e}\tGr   {:.2e}\tExud {:5.2e}'.format(sum(Q_Rm), sum(Q_Gr), sum(Q_Exud))) Q_labile
                print('init\tST  {:.2e}\tmeso   {:.2e}'.format(sum(Q_ST_init), sum(Q_meso_init)))
                print("aggregated sink satisfaction at last time step (%)\n\ttot {:5.1f}\tRm {:5.1f}\tGr {:5.1f}".format(
                     div0f(sum(Q_Rm_i+ Q_Gr_i),sum(Q_Rmmax_i + Q_Grmax_i))*100,div0f(sum(Q_Rm_i),sum(Q_Rmmax_i))*100, 
                     div0f(sum(Q_Gr_i), sum(Q_Grmax_i), 1.)*100))
                print("aggregated sink repartition at last time step (%):\n\tRm {:5.1f}\tGr {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm_i)/sum(Q_out_i)*100, 
                     sum(Q_Gr_i)/sum(Q_out_i)*100,sum(Q_Exud_i)/sum(Q_out_i)*100))
                print("aggregated sink repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rm)/sum(Q_out)*100, 
                     sum(Q_Gr)/sum(Q_out)*100,sum(Q_Exud)/sum(Q_out)*100))
                print("aggregated sink repartition for max (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}".format(sum(Q_Rmmax_i)/sum(Q_outmax_i)*100, 
                     sum(Q_Grmax_i)/sum(Q_outmax_i)*100,sum(Q_Exudmax_i)/sum(Q_outmax_i)*100))
                print("aggregated C repartition (%) :\n\tRm   {:5.1f}\tGr   {:5.1f}\tExud {:5.1f}\tLabile {:5.1f}".format(sum(Q_Rm)/sum(Q_out + Q_labile)*100, 
                     sum(Q_Gr)/sum(Q_out + Q_labile)*100,sum(Q_Exud)/sum(Q_out + Q_labile)*100,sum(Q_labile)/sum(Q_out + Q_labile)*100))
                print("aggregated Q_labile repartition (%) :\n\tQ_ST   {:5.1f}\tQ_S_ST   {:5.1f}\tQ_meso {:5.1f}".format(sum(Q_ST)/sum(Q_labile)*100, 
                     sum(Q_S_ST)/sum(Q_labile)*100,sum(Q_meso)/sum(Q_labile)*100))
            print("Error in Aux_balance:\n\tabs (mmol) {:5.2e}\trel (-) {:5.2e}".format(
                    errorAuxin, div0f(errorAuxin,sum(Q_Auxin_stem), 1.)))
        
        if (max(C_ST) > 3):
            stringError = "ERROR max(C_ST) > 3 for thread "+ str(thread)
            print(stringError)
            write_file_array("stringError", np.array([stringError])) 
            raise Exception(stringError)
            
        ###
        #4UQ
        #
        ###
        Mstem = r.plant.getOrgans(3, True)[0]#main shoot
        
        orgsF = r.plant.getOrgans(3, False)
        orgs = np.array([Mstem.getChild(i) for i in range(Mstem.getNumberOfChildren())])#r.plant.getOrgans(3, True)
        toKeep = np.array([org.organType() == 3 for org in orgs])
        orgs = np.array(orgs)[toKeep]
        toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
        orgs = np.array(orgs)[toKeep]
        
        temp_time = time.time()
        
        dtreal = ( time.time() - dt_lastWrote)
        dtsimsim = simDuration - dt - dtSIM_lastWrote
        if(not r.burnInTime) and (doDict or ((dtreal >= dt_write) and (dtsimsim >= dtSIM_write) )):
            write_file_array("time_all", np.array([thread, simDuration,simMax,temp_time , start_time,temp_time - start_time]))
            write_file_float("time", simDuration)
            #write_file_array("OutAuxin", OutAuxin)
            #write_file_array("Q_Auxin", Q_Auxin)
            write_file_array("C_Auxin", C_Auxin)
            auxTested= np.array([org.auxTested for org in orgs])
            write_file_array("auxTested", auxTested)
            write_file_array("InAuxinAll", np.array(r.cpb_2_pm.AuxinSource)* r.auxin_P * dt)
            #write_file_array("Delta_JA_ST", Delta_JA_ST)
            #write_file_array("JAuxin_ST2", JAuxin_ST2)


            #write_file_float("doDecapitation", doDecapitation)
            #write_file_float("simMax", simMax)
            #write_file_float("numLNodes",numLNodes)
            #id_orgs = [np.full(org.getNumberOfNodes()-1,org.getId()) for org in orgs]
            #id_orgs = np.array([item for sublist in id_orgs for item in sublist])
            #write_file_array("id_orgsPerNodes", id_orgs)

            #write_file_array("idnodeOfOrg", idnodeOfOrg)
            
                
            #write_file_array("length_org", length_org)
            write_file_array("lengthth_org", lengthth_org)
            
            if((min(lengthth_org[:2]) == 0) and (min(maxLBudDormant)>0)):
                print("min(stemlengths) == 0)",lengthth_org )
                raise Exception
                
            #raise Exception
            #write_file_array("parentOrgId", parentOrgId)
            write_file_array("ownOrgId", ownOrgId)
            write_file_array("distFromParentBase", distFromParentBase)

            write_file_array("sucTested", sucTested)#only take shoot C
            write_file_array("sucTestedWeighted",sucTested*lengthth_org)
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
            parentLinkingNode = np.array([org.parentLinkingNode for org in orgs]) 
            write_file_array("parentLinkingNode", parentLinkingNode)
            orgsLeaves = r.plant.getOrgans(4, True)
            
            stemNodes = np.array(Mstem.getNodeIds())
            write_file_array("Grmax_i_s1", Q_Grmax_i[stemNodes])
            write_file_array("Gr_i_s1",Q_Gr_i[stemNodes])
            write_file_array("Rm_i_s1",Q_Rm_i[stemNodes])
            write_file_array("C_ST_s1",C_ST[stemNodes])
            write_file_float("Qexudi",sum(Q_Exud_i))
            write_file_float("Qrmi",sum(Q_Rm_i))
            
            if __name__ == '__main__':
                print("numLNodes",numLNodes,simMax)
                print("sucTested",sucTested)
                print("auxTested",auxTested)
                print("lengthth_buds",lengthth_org)
                print("buds_age",buds_age)
                print("BerthFact",BerthFact)
                print("budStage",budStage)
                print("sum aux source", sum(AuxinSource))
                print("stem length",lln)
                print("growth rate",(lln - lln_old)/dt)
                print("leaf_th length",ll)
                #print('true', ll_true)
                print("growth rate",(ll - ll_old)/dt)
                print('leaf ages',leaf_ages)
                #print("get ids",idRoots)
                print("root_th length",sum([value.item() for key, value in lr.items()]))
                gr_real = []
                for key, value in lr.items():
                    if key in lr_old:
                        gr_real.append((value.item() - lr_old[key].item())/dt)
                    else:
                        gr_real.append(value.item()/dt)
                print("growth rate",sum(gr_real))
                print("growth_th",sum(grr_th))
                #print("leaf_obs length",np.array([ll.getLength(True) for ll in leaves]))
                print("leaf",nodeD,growthUpToNode,len(leaves),leafArea,maxLeafArea)
                print("do decapitate?",numLNodes, growthUpToNode, (not r.burnInTime), nodeD  , (not changedSimMax))
            if(changedSimMax):
                timeSinceDecap = simDuration - (simMax - testTime)
                if((not (budStage[(nodeD):] ==-1).all()) and (nodeD>0)):
                    print(thread, "not (arr[(nodeD+1):] ==-1).all()")
                    print(budStage,nodeD,budStage[(nodeD+1):] ,(budStage[(nodeD):] ==-1),(not (budStage[(nodeD):] ==-1).all()))
                    errorMessage = str(thread)+" not (arr[(nodeD+1):] ==-1).all() "
                    erM2 = str(budStage) +" "+str(nodeD) +" "+str(budStage[(nodeD+1):]) \
                        +" "+str((budStage[(nodeD):] ==-1)) +" "+str((not (budStage[(nodeD):] ==-1).all()))
                    raise Exception(errorMessage + erM2)
            else:
                timeSinceDecap = -1
            
            outcondition = doCondition(r,timeSinceDecap, thread,(temp_time - start_time)/(60*60*24), outcondition, nodeD,allInputsdict)
            ctoohigh = (budStage[0]==2)and (max(C_ST)>3)
            if (((temp_time - start_time)/(60*60*24) > 2) or (outcondition != 0) or ctoohigh): #success or falur
                simMax = -1
                
            if len(orgs) != len(budStage):
                print(budStage, len(budStage))
                print(ot_orgs,len(ot_orgs))
                print(len(orgs))
                raise Exception
            if not doDict:
                dt_lastWrote = time.time()
                dtSIM_lastWrote = simDuration - dt
                
            stem = r.plant.getOrgans(3, False)[0]
        
        ### size of phytomeres
        

        Mstem = r.plant.getOrgans(3, True)[0]#main shoot
        orgs = np.array([Mstem.getChild(i) for i in range(Mstem.getNumberOfChildren()) if Mstem.getChild(i).organType() == pb.stem])#r.plant.getOrgans(3, True)
        lengthth_org = np.array([org.getLength(False) for org in orgs])
        buds_age = np.array([org.getAge() for org in orgs])
        length_org = np.array([org.getLength(True) for org in orgs])
        parentOrgId = np.array([org.getParent().getId() for org in orgs])
        ownOrgId = np.array([org.getId() for org in orgs])
        distFromParentBase = np.array([org.getParent().getLength(org.parentNI) for org in orgs])
        sucTested = np.array([org.sucTested for org in orgs])
        
        lln_old = lln
        epsilonDxPerPhyto = np.array(stem.epsilonDxPerPhyto)
        lln = stem.getLlocalId_linking_nodes()
        lln = np.diff([stem.getLength(lln_i) for lln_i in lln]) + epsilonDxPerPhyto 
        numLNodes = np.sum(lln > 2.6)
        leaves = np.array([stem.getChild(np) for np in range(stem.getNumberOfChildren()) if stem.getChild(np).organType() == pb.leaf])
        ll_old = ll
        ll = np.array([ll_.getLength(False) for ll_ in leaves])    
        ll_true = np.array([ll_.getLength(False) for ll_ in leaves])             
        roots = r.plant.getOrgans(2, True)
        lr_old = lr
        lr = np.array([ll_.getLength(False) for ll_ in roots])
        idRoots = np.array([ll_.getId() for ll_ in roots])
        lr = dict(zip(idRoots, lr))
        leaf_ages = np.array([ll_.getAge() for ll_ in leaves]) 
        
        
        leafArea = leaves[-1].getLength(False) * leaves[-1].getParameter("Width_blade")  
        maxLeafArea = leaves[-1].getParameter("k") * leaves[-1].getParameter("Width_blade") 
        
        forDecapitate = (numLNodes >= growthUpToNode) and (not r.burnInTime) and (nodeD !=0)  and (not changedSimMax)
        
        if  (numLNodes >= growthUpToNode) and (not r.burnInTime) and (nodeD !=0) and (not changedSimMax): 
            print('do decapitate')
            raise Exception
            org_ = r.plant.getOrgans(3, False)[0] #get first (==main) stem
            org_.active = False
            org_.budStage = -1
            lastNode = leaves[nodeD -1].parentNI-1
            toKil = np.array([org_.getNodeId(nnid) for nnid in range(lastNode,org_.getNumberOfNodes())])#-1
            kid_pni = np.array([org_.getChild(kkiidd).getNodeId(0) for kkiidd in range(org_.getNumberOfChildren())])
            kid_id = np.array([org_.getChild(kkiidd).getId() for kkiidd in range(org_.getNumberOfChildren())])
            
            if __name__ == '__main__':
                print(toKil,org_.getNumberOfChildren())
                print("kid_pni",kid_pni,kid_id)
                
            if len(kid_pni) > 0:
                selectKids =np.concatenate(([np.where( kid_pni ==toKill_)[0] for toKill_ in toKil]))#,np.where(kid_pni==toKil[1])[0]
                if __name__ == '__main__':
                    print("selectKids",selectKids)
                dying_kids = np.array([org_.getChild(kkiidd) for kkiidd in selectKids])
                k_ =[killChildren(ni) for ni in dying_kids]
                toFill_ = [item for sublist in k_ for item in sublist]
                toFill_= np.concatenate((toFill_,toKil))
                toKil  = np.array(np.unique(toFill_),dtype = np.int32)
            r.plant.node_Decapitate = toKil[1:]
            doDecapitation = False
            if __name__ == '__main__':
                print("toKil end",toKil)
            if not changedSimMax:
                if simMax_ < 0:
                    simMax = simDuration + testTime #end 7 days after decapitation
                changedSimMax = True
                dt = dtAfter #1MIN
                toAdd = np.ceil(simDuration) - simDuration #to have start of day after decapitation
        #assert numLNodes <= growthUpToNode
        #if doDecapitation and (numLNodes > nodeD): 
         #   raise Exception("too many linking nodes")
        
        if (not doDecapitation) and (not changedSimMax) and (nodeD ==0) and (numLNodes >=growthUpToNode): 
            if simMax_ <0:
                simMax = simDuration + testTime
            changedSimMax = True
            dt = dtAfter #1MIN
            #raise Exception
            
        ###
        #
        #4UQ
        #
        ###
        
        if doVTP == 2:
            ana = pb.SegmentAnalyser(r.plant.mappedSegments())
            
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
            C_meso
            
            C_meso_p = C_meso
            C_meso_p[abs(C_meso_p) < cutoff] = 0
            
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
            
            #AuxinSource_p = np.array(r.AuxinSource)
            
            ana.addData("C_meso", C_meso_p)
            ana.addData("JAuxin_ST2", JAuxin_ST2_p)
            ana.addData("C_AuxinOut", C_AuxinOut_p)
            ana.addData("Delta_JA_ST", Delta_JA_ST_p)
            ana.addData("AuxinSource", AuxinSource)
            ana.addData("isRootTip", np.array(r.isRootTip))
            ana.addData("C_Auxin", C_Auxin_p)
            
            #ana.addData("activePhloem",activePhloemv) 
            #ana.addData("activeAuxin",activeAuxinv)  
            #ana.addData("budStage",budStage)             
            
            ana.addData("CST", C_ST_p)
            #do as konz or div per vol or surf?
            #ana.addData("Q_Exud", Q_Exud)  # cut off for vizualisation
            ana.addData("fluxes", fluxes_p)
            #ana.addData("Fpsi", np.array(r.Fpsi))
            
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
            #write_file_array("Ag4Phloem_p", r.Ag4Phloem_p)
            
            #try: 
            #    ana.addData("Ag4Phloem", Ag4Phloem_p)
            #except:
            #    errorMsg = str(thread)+" "+str(len(Ag4Phloem_p))+" "+str(len(Q_Exud_i))+" "+str(len(psiXyl_p))
            #    raise Exception(errorMsg)
            # 
            
            ana.write("results"+directoryN+"plotplant_"+dir4allResults+"_"+ str(ö) +".vtp", 
                      ["CST", "fluxes","psi_Xyl",
                        "C_Auxin","AuxinSource","Delta_JA_ST",
                        "C_AuxinOut","JAuxin_ST2",
                        #"activePhloem",
                        #"activeAuxin",#"Ag4Phloem",
                        "isRootTip","C_meso",
                        "QExud", "QGr", "QRm",
                        "CExud", "CGr", "CRm",
                        "QExudmax", "QGrmax", "QRmmax",
                        "CExudmax", "CGrmax", "CRmmax",
                        "organType", "subType"#, "Fpsi"
                        ]) 
        
        
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
 
        ö +=1
        
        
        if not r.burnInTime:
            grr_th = np.array([ll_.getTheoreticalGrowth(dt) for ll_ in roots])/dt
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
            mainStemAux_std = np.max(np.abs((mainStemAux[1:]-mainStemAuxBU)))
            mainStemAuxBU = mainStemAux[1:]
            mainStemAux_mean = np.mean(mainStemAux[1:])
            if __name__ == '__main__':
                print("mainStemAux_std",mainStemAux_std, mainStemAux_mean, mainStemAux)
        else :
            simDuration += dt
            
        # if r.burnInTime and \
        #  ((activeAtThreshold_suc and (min(C_ST)>0.)) or \
        #  (activeAtThreshold_auxin and (mainStemAux_std < 0.01) and (mainStemAux_mean > 1) )):
        r.burnInTime = False
        if r.burnInTime and (mainStemAux_std < 0.01) and (mainStemAux_mean > 1) : 
            r.burnInTime = False
            r.auxin_init_mean = mainStemAux_mean
            
            with open('results'+ directoryN+"input"+dir4allResults+ '.csv', 'r+') as f:
                lines = f.read().splitlines()
            lines[0] = lines[0] + ",mainStemAux_mean"
            lines[1] = lines[1] + "," + repr(mainStemAux_mean)
            with open('results'+ directoryN+"input"+dir4allResults+ '.csv', 'w') as f:
                f.write(lines[0])
                f.write("\n")
                f.write(lines[1])
                
            errorBU = error
            Q_outTemp = np.full(len(r.Q_out),0.)
            Q_outTemp[(Ntbu*9):(Ntbu*10)] = Q_Auxin
            Q_outTemp[0:Ntbu] = Q_ST   
            Q_outTemp[Ntbu:(Ntbu*2)] = Q_meso 
            r.Q_out = Q_outTemp
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
            
            orgs = r.plant.getOrgans(3, True)
            toKeep = np.array([org.getParameter("subType") <= 2 for org in orgs])
            orgs = np.array(orgs)[toKeep]
            #for org in orgs:
             #   org.auxTested = 2.49330566
            #print(np.array([org.auxTested for org in orgs]))
            #raise Exception
        
            
        if doDict and (( time.time() - dt_lastWrote) >= dt_write) and (( simDuration - dt - dtSIM_lastWrote ) >= dtSIM_write)  :
            #print("write",time.time())
            for keykey in outputsDict_array.keys():
                name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.csv'
                datas = outputsDict_array[keykey]
                with open(name2, 'a') as log:
                    for data in datas:
                        log.write(','.join([num for num in map(str, data)])  +'\n')
            for keykey in outputsDict_float.keys():
                name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.csv'
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
            name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.csv'
            datas = outputsDict_array[keykey]
            with open(name2, 'a') as log:
                for data in datas:
                    log.write(','.join([num for num in map(str, data)])  +'\n')
        for keykey in outputsDict_float.keys():
            name2 = 'results'+ directoryN+ keykey+ dir4allResults+ '.csv'
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
        
    delete_file("running")
    outId = -1
    mySuccessPoints = successPoints(r,timeSinceDecap, thread,(temp_time - start_time)/(60*60*24), outcondition, nodeD,allInputsdict)
    write_file_float("successPoints", mySuccessPoints)
    if (outcondition >= 0) and (changedSimMax):
        outId = thread
        print(thread, outId,nodeD,"success", budStage,simMax, lengthth_org,mySuccessPoints)
    if not changedSimMax:
        print(thread, outId,nodeD,"fail, not reached changedSimMax", budStage, simMax,mySuccessPoints)
    #print("finished", thread, time.time() - start_time)
    #os._exit(os.EX_OK)
    del r
    del pl
    gc.collect()
    #import psutil
    #print("intern pid",psutil.Process().pid)#, r.burnInTime)
    return outId


if __name__ == '__main__':
    startDate = int(sys.argv[1])
    
    start_time_ = time.time()
    print("intern pid_start",psutil.Process().memory_info())
    from CalibP1Database import toTry
    #from CalibP1Database import doCondition_
    params = toTry()
    Qsv=params['Qsv']
    MulimSucv=params['MulimSucv']
    nodeDv=0
    GrRatiov= params['GrRatiov']
    CarbonCostv= params['CarbonCostv']
    Klightv= params['Klightv']
    i = 1
    totrun = 255*2
    
    try:
        nameSim = sys.argv[1]
    except:
        nameSim = "WTD"
    try:
        lightLevel = sys.argv[2]
    except:
        lightLevel = "high"
    main_dir=os.environ['PWD']#dir of the file
    directoryN = "/"+os.path.basename(__file__)[:-3]+"_"+str(nameSim)+"/"#+"_"+str(lightLevel)"/a_threshold/"
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
    
    assert lightLevel in np.array(["low","medium","high"])
    Qinit = 1000*1e-6 
    if lightLevel == "medium":
        Qinit *= 0.86 
    elif lightLevel == "low":
        Qinit *= 0.64
    fileparam_ = "UQ_1LeafRS"
        
    BerthLim_ =  3.
    thresholdSuc_ = 1#3e-2 #*2*2
    simMax__ = -1
    runSim(startDate, directoryN_ = directoryN, doVTP = 0, verbosebase = False,
    PRate_ = 6.8e-3, thresholdAux = 0, 
     RatiothresholdAux = 1,useLength = 2,limLenActive_ = 1.5,
     Qmax_ = Qinit, Klight = 0.,
     thresholdSuc = thresholdSuc_, 
     maxLBud = np.array([1.,1.,1.,1.]),  maxLBudDormant = np.array([0.15,0.15,0.15,0.15]), #([0.1,0.15,0.05])
     L_dead_threshold=2000.,
     kss=0.2,kaa=1,
    BerthLim =BerthLim_,
     nodeD = 7, thread = 1,
     testTime=9, dtBefore = 1/24, dtAfter= 30/(60*24),
    start_time = start_time_,
     doPrint = True, doDict = False,
     dt_write = 0, dtSIM_write = 30/(60*24),auxin_D=0.,
    #doCondition = doCondition_
           simMax_ = simMax__,
           fileparam =fileparam_,
           leafAsIAASource_ = True, growthUpToNode =9
          )
    end_time_ = time.time()
    print(end_time_ - start_time_ )
    sys.exit(0)
    
def runAllSim(i, arg1,arg2):
    start_time_ = time.time()   
    
    nameSim = arg1
    lightLevel = arg2
    main_dir=os.environ['PWD']#dir of the file
    directoryN = "/"+os.path.basename(__file__)[:-3]+"_"+str(nameSim)+"_"+str(lightLevel)+"/"#"/a_threshold/"
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
    
    assert nameSim in np.array(["DW","SLM","WT","WTD"])
    assert lightLevel in np.array(["low","medium","high"])
    Qinit = 500*1e-6 
    GrRatio_ = 10
    GrRatioLats_ = 5
    GrRatioLeaf = 10
    if lightLevel == "medium":
        Qinit *= 0.86 
    elif lightLevel == "low":
        Qinit *= 0.64
        
    if nameSim == "DW":
        fileparam_ = "UQ_1Leaf_DW"
        Qmax__ = Qinit*0.5 #- 25e-6 #175*1e-6
        GrRatio_ = GrRatio_/3
        GrRatioLats_ = GrRatioLats_/3
        GrRatioLeaf = GrRatioLeaf/3
    else:
        fileparam_ = "UQ_1Leaf"
        Qmax__ = Qinit
    if nameSim == "SLM":
        BerthLim_ = 1000
        thresholdSuc_ = 0
    else:
        BerthLim_ =  3.
        thresholdSuc_ = 3e-2 #*2*2
    if nameSim == "WTD":
        nodeD_ = 7
        simMax__ = -1
    else:
        nodeD_ = 0
        simMax__ = 25 + 10/24#1.2708333333343
    runSim(directoryN_ = directoryN, doVTP = 0, verbosebase = False,
    PRate_ = 6.8e-3,  thresholdAux = 0, 
     RatiothresholdAux = 1,useLength = 2,limLenActive_ = 0.5,
     Qmax_ = Qmax__, Klight = 0.,
     thresholdSuc = thresholdSuc_, 
           maxLBud = np.array([1.]),  
           maxLBudDormant = np.array([0.05,0.15,0.05]), #([0.1,0.15,0.05])
     budGR = 0.1,L_dead_threshold=2000.,
     kss=0.2,kaa=1,
    BerthLim =BerthLim_,
     nodeD = nodeD_, thread = i,
     testTime=7, dtBefore = 1/24, dtAfter= 30/(60*24),
    start_time = start_time_,
     doPrint = True, doDict = False,
     dt_write = 0, dtSIM_write = 30/(60*24),auxin_D=0.,
    #doCondition = doCondition_
           simMax_ = simMax__,
           fileparam =fileparam_,
           leafAsIAASource_ = True
          )
    end_time_ = time.time()
    #print(end_time_ - start_time_ )