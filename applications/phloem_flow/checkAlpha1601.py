""" water movement within the root (static soil) """


import sys; 
#directoryN = "/"+sys.argv[0].split('.')[0]+"/"
sys.path.append("../.."); sys.path.append("../../src/python_modules")
CPBdir = "../.."
#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

from phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import vtk_plot as vp
#import matplotlib.pyplot as plt
from datetime import datetime, timedelta

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
                                     
    r.setKr([[kr_r0],[kr_s],[kr_l]], kr_length_=l_kr)
    r.setKx([[kz_r0],[kz_s],[kz_l]])
    
    
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
    kroot = 5e-4
    kr_r0 = kroot
    kr_r1 = kroot
    kr_r2 = kroot
    kr_r3 = kroot
    l_kr = 0.8 #cm
    
    r.setKr_st([[kr_r0],[kr_s],[kr_l]] , kr_length_= l_kr)
    r.setKx_st([[kz_r0],[kz_l],[kz_l]])
    
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
    r.setAcross_st([[Across_s_r0],[Across_s_s],[Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])
    

def runSim(directoryN_,auxin_alpha =0.0027, thread = 0,verbosebase = True, Array2d_result = np.array([]), doLogalpha = True,
          parameterFile = "setAlpha", dt_ = 1/(24*60*60),
          dt_write = 1/24):
    PRate_ = 1
    threshold = -40
    doVTP = False
    nodeD = 6
    simInit= 3
    testTime = 10/24 
    doDecapitation = (nodeD > 0)
    directoryN = directoryN_
    strPRate_ = str(int(np.round(PRate_)))
    strTh = str(int(threshold*10))
    strDecap = str(nodeD)#[0]
    
    #print("in runSim")
    
    if len(Array2d_result)==0:
        CSVData = open('fig2bRenton2012.txt')
        Array2d_result_ = np.loadtxt(CSVData, delimiter="\t")
    
    def write_file_array(name, data, doLogalpha_=doLogalpha):
        if doLogalpha_:
            name2 = 'results'+ directoryN+ name+ "_"+ str(thread)+'.txt'
            with open(name2, 'a') as log:
                log.write(','.join([num for num in map(str, data)])  +'\n')
    def open_file_array(name):
        name2 = 'results'+ directoryN+ name+ "_"+ str(thread) +'.txt'
        return np.loadtxt(open(name2, 'r'), delimiter=",") 

    def write_file_float(name, data, doLogalpha_=doLogalpha):
        if doLogalpha_:
            name2 = 'results' + directoryN+  name+ "_"+  str(thread)+ '.txt'
            with open(name2, 'a') as log:
                log.write(repr( data)  +'\n')
            
    Qmax_ = 1000e-6        
    weatherInit = weather(0, Qmax_)
    
    simDuration = 0#simInit # [day] init simtime
    depth = 60
    dt = 1#dt
    verbose = False

    # plant system 
    pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
    path = CPBdir+"/modelparameter/plant/"
    name = parameterFile

    pl.readParameters(path + name + ".xml")


    #raise Exception
    sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

    pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


    #pl.activeAtThreshold_auxin = True
    #pl.activeAtThreshold = False
    pl.initialize(verbose = False)#, stochastic = False)
    pl.simulate(1, False)
    simDuration +=1
    mainstem = pl.getOrgans(3, False)[0]
    #print(mainstem.getLength(), mainstem.getNumberOfLinkingNodes())
    #print([mainstem.getLength(i) for i in range(mainstem.getNumberOfNodes())])
    assert(mainstem.getLength()==17.)
    maintstemLen = np.array([mainstem.getLength(i) for i in range(mainstem.getNumberOfNodes())])
    lenlen2 = np.array([0.5*i for i in range(int(17.5*2))])
    print(maintstemLen, lenlen2)
    assert(sum(maintstemLen == lenlen2)==len(maintstemLen))
    
    # while (len(pl.getOrgans(3, False)) ==0 ) or (pl.getOrgans(3, False)[0].getNumberOfLinkingNodes() <6):
    #     pl.simulate(dt, False)#, "outputpm15.txt")
    #     simDuration += dt
    #     if len(pl.getOrgans(3, False)) !=0:
    #         print("pl.getOrgans(3, False)[0].getNumberOfLinkingNodes()",pl.getOrgans(3, False)[0].getNumberOfLinkingNodes())
    simMax = 100#simDuration + testTime #if !doDecapitation

    """ Coupling to soil """
    min_b = [-3./2, -12./2, -61.]#distance between wheat plants
    max_b = [3./2, 12./2, 0.]
    cell_number = [6, 24, 61]#1cm3? 
    layers = depth; soilvolume = (depth / layers) * 3 * 12
    k_soil = []
    initial = weatherInit["p_mean"]#mean matric potential [cm] pressure head


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
    grRate = [[4.],[2.],[3.]]
    grRate = np.array([np.array(xi) for xi in grRate],dtype=object)
    grRate[1:] *= 10 #[1:]

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
    r.CSTimin = 0.05#4
    #r.surfMeso=0.0025
    r.cs = weatherInit["cs"]

    r.expression = 6
    r.update_viscosity = True
    r.solver = 1
    r.atol = 1e-14
    r.rtol = 1e-10
    #r.doNewtonRaphson = False;r.doOldEq = False
    SPAD= 31.0
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
        structSumInit += org.orgVolume(-1,False) * r.rhoSucrose_f(int(org.getParameter("subType")-1),org.organType())

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


    beginning = datetime.now()
    #1h for 1d when dxMin = 0.3

    AnSum = 0
    Q_ST_init = np.array([])
    Q_meso_init  = np.array([])
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
    Q_Parbu      = np.array([0.])
    Q_mesobu    = np.array([0.])


    r.doTroubleshooting = False

    #print("r.plant.activeAtThreshold_auxin,r.auxin_threshold,r.auxin_alpha:",r.plant.activeAtThreshold_auxin,r.auxin_threshold,r.auxin_alpha)
    #print("r.initValAuxin",r.initValAuxin)
    r.auxin_threshold = threshold
    r.auxin_D =0# 0.5 #e-6#muM /d  #e-6#3e-7
    r.auxin_P = PRate_ #muM /d  #*1e-6# 6e-7
    minAlhpa = 1e-4
    maxAlpha = 0.075
    r.auxin_alpha =auxin_alpha#5#e-6# 2.4e-3
    r.initValAuxin = 0 #0.2
    r.deleteAtRootTip = True
    r.plant.useCWGr = False
    InAuxin = 0
    #print("r.plant.activeAtThreshold_auxin,r.auxin_threshold,r.auxin_alpha:",r.plant.activeAtThreshold_auxin,r.auxin_threshold,r.auxin_alpha)
    #print("r.plant.activeAtThreshold,r.CSTthreshold,r.canStartActivating:",r.plant.activeAtThreshold,r.CSTthreshold,r.canStartActivating)
    #print("r.initValAuxin",r.initValAuxin)
    r.canStartActivating = False
    r.CSTthreshold = threshold #0.8#[[-1],[-1,0.3,-1],[-1]]
    orgs_stems = r.plant.getOrgans(3, True)
   
    def killChildren(orgToKill):
        toFill = orgToKill.getNodeIds()
        orgToKill.alive = False
        if orgToKill.getNumberOfChildren() > 0:
            toFill_ = [killChildren(orgToKill.getChild(ni)) for ni in range(orgToKill.getNumberOfChildren())]
            toFill_ = [item for sublist in toFill_ for item in sublist]
            toFill = toFill + toFill_
        return toFill
    
    #todo: adapt distance selctionso that does not matter when change resolution. 
    #test for different time step and check that get same result.
    
        
    day4burnin = simMax
    r.burnInTime = True
    LL = -1
    changedSimMax = False
    calibrateD = np.array([2*i for i in range(9)]) +1 #cm
    calibrateT = np.array([i/24 for i in range(11)])*24*60
    timeAlpha = 0
    while  timeAlpha < 10/24:
        #print('simDuration:',simDuration )
        
        numLNodes = r.plant.getOrgans(3, False)[0].getNumberOfLinkingNodes()
        if  not r.burnInTime: #(numLNodes >= nodeD) 
            org_ = r.plant.getOrgans(3, False)[0] #get first (==main) stem
            org_.active = False
            LL = org_.getLength()
            toKil = np.array([org_.getNodeId(org_.getNumberOfNodes()-1)])
            kid_pni = np.array([oo.getNodeId(0) for oo in r.plant.getOrgans(-1, True)])
            #print("kil1 ", toKil, kid_pni)
            #print(np.where( kid_pni ==toKil[0] )[0])
            if len(kid_pni) > 0:
                selectKids =np.where( kid_pni ==toKil[0] )[0]
                #print("selectKids ", selectKids)
                k_ =[killChildren(np.array(r.plant.getOrgans(-1, True))[ni]) for ni in selectKids]
                toFill_ = [item for sublist in k_ for item in sublist]
                toFill_= np.concatenate((toFill_,toKil))
                toKil  = np.unique(toFill_)
            r.plant.node_Decapitate = toKil
            #print(toKil)
            #raise Exception
            doDecapitation = False
            if not changedSimMax:
                simMax = simDuration + testTime #end 7 days after decapitation
                changedSimMax = True
                dt = 1/(24)#1hr
                #mainStemAuxRatio = (mainStemAux/mainStemAux_mean)[idToKeep]
                #write_file_array("mainStemAuxRatio", mainStemAuxRatio,doLogalpha_= True)
                timeAlphaInit = simDuration
                timeAlpha = simDuration - timeAlphaInit
            #print("toKil",toKil)
            if __name__ == '__main__':
                print(mainStemAux)
                print(mainStemAux_std,mainStemAux_mean)

        if doDecapitation and (numLNodes > nodeD): 
            raise Exception("too many linking nodes")
        
        if (not doDecapitation) and (not changedSimMax) and (nodeD ==0) and (numLNodes ==8): 
            simMax = simDuration + testTime
            changedSimMax = True
            
        weatherX = weather(simDuration, Qmax_)

        r.Qlight = weatherX["Qlight"] #; TairC = weatherX["TairC"] ; text = "night"
            
            
        setKrKx_xylem(weatherX["TairC"], weatherX["RH"])
        
        #r.maxLoop = 3
        r.solve_photosynthesis(sim_time_ = simDuration, sxx_=sx, cells_ = True,RH_ = weatherX["RH"],
            verbose_ = False, doLog_ = False,TairC_= weatherX["TairC"] )
        ##print(r.psiXyl)
        #raise Exception
        ots = np.concatenate((np.array([0]), r.get_organ_types()))#per node
        leavesSegs = np.where(ots[1:] ==4)
        fluxes = np.array(r.outputFlux)
        fluxes_leaves = fluxes[leavesSegs]
        errLeuning = sum(r.outputFlux)
            
        """ dumux """    
        fluxesSoil = r.soilFluxes(simDuration, r.psiXyl, sx, approx=False)
        r.Ag4Phloem = np.full(len(r.Ag4Phloem),0)
        AnSum += np.sum(r.Ag4Phloem)*dt
        startphloem= simDuration
        endphloem = startphloem + dt
        stepphloem = 1
        filename = "results/pmincpb_" + str(simDuration) + "_15pm.txt" 
        
        # save_stdout = sys.stdout
        # sys.stdout = open('trash', 'w')
        verbose_phloem = False
        
        filename = "results"+ directoryN +"inPM_"+str(ö)+ "_"+ strPRate_ + "_"+strTh+"_"+strDecap+".txt"
        #print("startpm")
        
        r.useStemTip = True
        r.startPM(startphloem, endphloem, stepphloem, ( weatherX["TairC"]  +273.15) , verbose_phloem, filename)
        if r.withInitVal and (len(Q_ST_init) ==0) :
            Q_ST_init = np.array(r.Q_init[0:Nt])
            Q_meso_init = np.array(r.Q_init[Nt:(Nt*2)])
        # print(r.deltaSucOrgNode)
        # orgs_test = r.plant.getOrgans(-1, True)
        # for org_ in orgs_test:
            # print(org_.organType(), org_.getNodeIds())    
        # raise Exception
        
        Q_ST    = np.array(r.Q_out[0:Nt])
        Q_meso  = np.array(r.Q_out[Nt:(Nt*2)])
        Q_Rm    = np.array(r.Q_out[(Nt*2):(Nt*3)])
        Q_Exud  = np.array(r.Q_out[(Nt*3):(Nt*4)])
        Q_Gr    = np.array(r.Q_out[(Nt*4):(Nt*5)])
        #Q_Gr4       = Q_Gr[np.where(ot_4phloem==4)[0]]#Q_Gr4     - Q_Gr4bu
        #Q_Gr3       = Q_Gr[np.where(ot_4phloem==3)[0]]#Q_Gr3     - Q_Gr3bu
        #Q_Gr2       = Q_Gr[np.where(ot_4phloem==2)[0]]#Q_Gr2     - Q_Gr2bu
        
        
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
        
        error   = sum(Q_ST +Q_Par+ Q_meso + Q_out + Q_STLost + Q_MesoLost)- Q_in - sum(Q_ST_init)  - sum(Q_meso_init)
        
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
        
        OutAuxin = np.array(r.Q_out[(Nt*8):(Nt*9)])
        Q_Auxin  = np.array(r.Q_out[(Nt*9):(Nt*10)])
        InAuxin  += sum(np.array(r.AuxinSource)) * r.auxin_P * dt
        Q_AuxinInit  = np.array(r.Q_init[(Nt*9):(Nt*10)])
        AuxinDecap   = np.array(r.AuxinLost)
        errorAuxin   = sum(Q_Auxin) - sum(Q_AuxinInit)- InAuxin +sum(OutAuxin) + sum(AuxinDecap)
        C_Auxin      = np.array(r.C_Auxin)
        AuxinSource  = np.array(r.AuxinSource)
        C_AuxinOut   = np.array(r.C_AuxinOut)
        JAuxin_ST2   = np.array(r.JAuxin_ST2)
        Delta_JA_ST  = np.array(r.Delta_JA_ST)
        """ATT
        
        # C_ST_ = C_ST[1:]
        
        """
        
        ###
        #4UQ
        #
        ###
        org_ = r.plant.getOrgans(3, False)[0]
        maintstemNodesId = np.array([org_.getNodeIds()])
        mainStemAux =  np.array(C_Auxin[maintstemNodesId])[0]   
        idToKeep = np.array([ms in calibrateD for ms in maintstemLen])
        maintstemLen = np.round(np.array([org_.getLength(i) for i in range(org_.getNumberOfNodes())])*10)/10
        
        if not r.burnInTime:
            timeAlpha  += dt
            #print(timeAlpha*24*60,(round(timeAlpha*24*60)  in  calibrateT))
            #if (timeAlpha - dt_lastwrote >= dt_write):
            if(round(timeAlpha*24*60)  in  calibrateT):
                mainStemAuxRatio = (mainStemAux/mainStemAux_mean)[idToKeep]
                mainStemAuxRatioAll = (mainStemAux/mainStemAux_mean)[1:]
                write_file_array("mainStemAuxRatioAll", mainStemAuxRatioAll,doLogalpha_= True)
                write_file_array("maintstemLenAll", maintstemLen[1:],doLogalpha_= True)
                write_file_array("mainStemAuxRatio",mainStemAuxRatio,doLogalpha_= True)
                write_file_array("maintstemLen", maintstemLen[idToKeep],doLogalpha_= True)
                #print(timeAlpha*24*60)
                write_file_float("timeAlpha", timeAlpha*24*60,doLogalpha_= True)
                dt_lastwrote  = timeAlpha
                
                
        ###
        #4UQ
        #
        ###
        
        if abs(error) > 1e-3:
            print("suc error too high")
            raise Exception    
        if abs(div0f(errorAuxin,InAuxin, 1.)) > 1e-3:
            print("auxin error too high")
            raise Exception    
        if min(C_ST) < 0.0:
            print("min(C_ST) < 0.0", min(C_ST),np.mean(C_ST),max(C_ST))
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
            org_ = r.plant.getOrgans(3, False)[0] 
            raise Exception
            
        ö +=1
        
        
        verbose_simulate = False
        #if not r.burnInTime:
         #   r.plant.simulate(dt, False)
            
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
            upperB = np.where(np.round(maintstemLen) == 17)[0][0].astype(int) +1
            lowerB = np.where(np.round(maintstemLen) == 1)[0][0].astype(int)
            mainStemAux_std = np.std(mainStemAux[lowerB:upperB])
            mainStemAux_mean = np.mean(mainStemAux[lowerB:upperB])
           
        else :
            simDuration += dt
            
        if (mainStemAux_std < 0.01) and (mainStemAux_mean > 1) and ( r.burnInTime):
            r.burnInTime = False
            dt = dt_
            r.canStartActivating = True
            Q_outTemp = np.array(r.Q_out)
            Q_outTemp[:] = 0
            Q_outTemp[(Nt*9):(Nt*10)] = Q_Auxin
            Q_outTemp[0:Nt] = Q_ST   
            Q_outTemp[Nt:(Nt*2)] = Q_meso 
            r.Qinit = Q_outTemp
            
            Q_Rmbu       =   np.full(Nt , 0.)
            Q_Grbu       =   np.full(Nt , 0.)
            Q_Exudbu     =   np.full(Nt , 0.) 
            Q_Rmmaxbu    =   np.full(Nt , 0.)
            Q_Grmaxbu    =   np.full(Nt , 0.) 
            Q_Exudmaxbu  =   np.full(Nt , 0.) 
            Q_STbu       =   np.full(Nt , 0.)
            Q_mesobu     =   np.full(Nt , 0.)
            mainStemAuxRatio = (mainStemAux/mainStemAux_mean)[idToKeep]
            mainStemAuxRatioAll = (mainStemAux/mainStemAux_mean)
            write_file_array("mainStemAuxRatioAll", mainStemAuxRatioAll[1:],doLogalpha_= True)
            write_file_array("maintstemLenAll", maintstemLen[1:],doLogalpha_= True)
            write_file_array("mainStemAuxRatio",mainStemAuxRatio,doLogalpha_= True)
            write_file_array("maintstemLen", maintstemLen[idToKeep],doLogalpha_= True)
            write_file_float("timeAlpha", timeAlpha*24*60)
            #print(timeAlpha*24*60)
            dt_lastwrote = timeAlpha
            
    if len(Array2d_result) == 0:
        CSVData = open('fig2bRenton2012.txt')
        Array2d_result = np.loadtxt(CSVData, delimiter="\t")
    
    outputSim = open_file_array('mainStemAuxRatio')
    test2 = outputSim - Array2d_result
    test2[np.isnan(test2)] = 0
    MSE = np.average(test2**2)
    return(MSE)

if __name__ == '__main__':
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
            
    MSE =runSim(directoryN,auxin_alpha =  0.002727303 , thread = 0,verbosebase = True, Array2d_result = np.array([]), doLogalpha = True,
          parameterFile = "setAlpha", dt_ = 1/24, dt_write = 1/(24*60))#0.00260063
    print(MSE)