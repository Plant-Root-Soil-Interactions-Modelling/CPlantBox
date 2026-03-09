""" water movement within the root (static soil) """


import sys; 
import gc
import math
import os
import numpy as np
import psutil

import time



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
    ea = specificHumidity * Pair / (0.378 * specificHumidity + 0.622)
    RH = qair2rh(specificHumidity, es, Pair)
    
    pmean =-100#theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,
                    'es':es,'ea':ea,
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
        
def doVTP(r):
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
                
    
    ana.addData("CST", C_ST_p)
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
                "C_Auxin","AuxinSource","Delta_JA_ST",
                "C_AuxinOut","JAuxin_ST2",
                "activePhloem","activeAuxin",#"Ag4Phloem",
                "isRootTip","C_meso"
                "QExud", "QGr", "QRm",
                "CExud", "CGr", "CRm",
                "QExudmax", "QGrmax", "QRmmax",
                "CExudmax", "CGrmax", "CRmmax",
                "organType", "subType", "Fpsi"]) 


    fluxes_p = np.insert(fluxes_p,0,0)# "[sucrose]",
    

def setPhloemflow_data(r): 

    r.C_targMesophyll = 2e-5 # todo: use data from new paper. same for starch
    r.k_S_Mesophyll = 1
    
    r.k_meso = 1e-3#1e-4
    r.cpb_2_pm.setKrm2([[2e-4]])#2e-4
    r.cpb_2_pm.setKrm1([[1.3e-1]])
    GrRatioLeaf = 10
    GrRatioRoot = 1
    CarbonCostRoot = 1/100
    CarbonCostLeaf = 1/10
    CarbonCost = 2
    rho_org = [[0.,1.34,1.34,1.34],[0.,1.44*CarbonCost,1.44*CarbonCost,1.44*CarbonCost],[0.,1.56,1.56,1.56,1.56]]#g C/gDW?
    
    density = 0.17 #g DW/cm3?
    #density and rho in mmol suc or carbon?
    rho_org = np.array([np.array(xi) for xi in rho_org],dtype=object)*density/2 #/2 => glucose to sucrose
    rho_org[2] *= CarbonCostLeaf
    r.cpb_2_pm.setRhoSucrose(rho_org)
    GrRatio = 10
    GrRatioLats = 5
    grRate = [[0.,4.,2.,1.0,4.],[0.,1.8*GrRatio,1.8*GrRatioLats,1.8*GrRatioLats],[0.,3.,3.,3.,3.]]
    grRate = np.array([np.array(xi) for xi in grRate],dtype=object)
    grRate[0] *= GrRatioRoot

    grRate[2] *= GrRatioLeaf

    r.cpb_2_pm.setRmax_st(grRate, verbose = True)
    r.KMfu = 0.2
    r.sameVolume_meso_st = True
    r.sameVolume_meso_seg = False
    r.withInitVal =True
    r.initValST = 0.#0.15
    r.initValMeso = 0.#0.2
    r.beta_loading = 0#3
    r.Vmaxloading = 0.15#0.3 #mmol/d, needed mean loading rate:  0.3788921068507634
    r.Mloading = 3.3*1e-3#0.2#
    r.Gr_Y = 1#0.75
    r.CSTimin = 0.25#0.3#
    r.update_viscosity = True
    r.atol = 1e-10# 1e-14
    r.rtol = 1e-6#1e-10

    #turn off effect of water limitation, TODO : reset again?
    r.psiMin = -10000000000*(1/0.9806806)
    r.psiMax = -1000000000*(1/0.9806806)
    
def setProtosynthesis_data(r): 
    # todo: take from paper?
    r.g0 = 8e-2
    r.VcmaxrefChl1 =2
    r.VcmaxrefChl2 =7
    r.a1 = 2#0.8/0.2
    r.a3 = 2.2
    #r.Rd_ref = 0#2e-6
    r.alpha = 0.27
    r.theta = 0.51
    SPAD= 60.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    r.Csoil = 1e-4

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh

def setKrKx_xylem(TairC,  r): #inC
    raise Excpetion
    #mg/cm3
    specificHumidity = 0.0097
    Pair = 1010.00 #hPa
    hPa2cm = 1/0.9806806
    es =  6.112 * np.exp((17.67 * TairC)/(TairC + 243.5))
    RH = qair2rh(specificHumidity, es, Pair)
    
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
    
    r.set_kr_distance_dependent([l_kr], [kr_l], subType = [1,2,3,4], organType = 4)
    r.set_kr_distance_dependent([l_kr], [kr_l], subType = [1,2,3,4], organType = 2)
    r.set_kr_distance_dependent([l_kr], [kr_l], subType = [1,2,3,4], organType = 3)
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ,kr_s,kr_s,kr_s,kr_s ,kr_s,kr_s ],[kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l]], kr_length_=l_kr)
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s,kz_s,kz_s,kz_s,kz_s,kz_s,kz_s ],[kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
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
    kroot = 5e-2
    kr_r0 = kroot
    kr_r1 = kroot
    kr_r2 = kroot
    kr_r3 = kroot
    l_kr = 0.8 #cm
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s, kr_s,kr_s,kr_s,kr_s, kr_s,kr_s],[kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l,kr_l]] , 
                kr_length = l_kr, verbose = False)
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_l,kz_l, kz_l,kz_l,kz_l,kz_l, kz_l,kz_l],[kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l,kz_l]])
    
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
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s,Across_s_s,Across_s_s],[Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l,Across_s_l]])
    #actually, don t use perimeter currently
    #r.setPerimeter_st([[0,Perimeter_s_r0,Perimeter_s_r12,Perimeter_s_r12,Perimeter_s_r0],[0,Perimeter_s_s,Perimeter_s_s],[0,Perimeter_s_l]])

def doConditionDefault(rinput, timeSinceDecap_,i, simtime, memoryCondition, nodeD):
    '''
    function which return True if the calibration has failed => did not respect the condition
    '''
    return False

def write_file_array(name, data, directoryN="/"):
    name2 = 'results'+ directoryN+ name+ '.csv'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data, directoryN="/"):
    name2 = 'results'+ directoryN+ name+ '.csv'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

def delete_file(name):
    name2 = 'results'+ directoryN+ name+ '.csv'
    my_file = Path(name2)
    if my_file.is_file():
        os.remove(name2)
    else:
        print("file running missing",thread, name2)



def killChildren(orgToKill):
    toFill = orgToKill.getNodeIds()
    orgToKill.alive = False
    orgToKill.budStage = -1
    if orgToKill.getNumberOfChildren() > 0:
        toFill_ = [killChildren(orgToKill.getChild(ni)) for ni in range(orgToKill.getNumberOfChildren())]
        toFill_ = [item for sublist in toFill_ for item in sublist]
        toFill = toFill + toFill_
    return toFill

    
    