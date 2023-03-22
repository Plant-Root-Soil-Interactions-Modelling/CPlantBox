import numpy as np

"""avoid division per 0 during post processing"""
def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c

"""compute soil water potential from water content and VanGenouchten parameters"""
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0] 
    thetas = vg[1] 
    alpha = vg[2] 
    n = vg[3] 
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2


"""get environmental conditions"""
def weather(simDuration):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1] #van Gemuchten parameters
    nightPart = 0.3
    Qnight = 0; Qday = 1000e-6 /nightPart #absorbed light
    Tnight = 15.8; Tday = 22 #temperature
    RHnight = 0.8; RHday = 0.5
    Pair = 1010.00 #hPa
    thetaInit = 30/100 #soil water content

    coefhours = sinusoidal(simDuration) #sinusoidal coeficient to get daily variation
    RH = RHnight + (RHday - RHnight) * coefhours
    TairC_ = Tnight + (Tday - Tnight) * coefhours
    Q_ = Qnight + (Qday - Qnight) *max(0, coefhours -nightPart)
    cs = 850e-6 #co2 paartial pressure at leaf surface (mol mol-1)
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5)) 
    ea = es * RH
    
    pmean = theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,'ea':ea,'es':es,
                    'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
    return weatherVar

def setKrKx_xylem(TairC, RH, r): #inC
    #mg/cm3
    hPa2cm = 1.0197
    dEauPure = (999.83952 + TairC * (16.952577 + TairC * 
        (- 0.0079905127 + TairC * (- 0.000046241757 + TairC * 
        (0.00000010584601 + TairC * (- 0.00000000028103006)))))) /  (1 + 0.016887236 * TairC)
    siPhi = (30 - TairC) / (91 + TairC)
    siEnne=0
    mu =  pow(10, (- 0.114 + (siPhi * (1.1 + 43.1 * pow(siEnne, 1.25) )))) 
    mu = mu /(24*60*60)/100/1000; #//mPa s to hPa d, 1.11837e-10 hPa d for pure water at 293.15K
    #water viscosity
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
    l_kr = 0.8 #cm
    r.setKr([[kr_r0,kr_r1,kr_r2,kr_r0],[kr_s,kr_s ],[kr_l]]) 
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    #log(-) * (cm^3*MPa/K/mol) * (K) *(g/cm3)/ (g/mol) * (hPa/MPa) * (cm/hPa) =  cm                      
    p_a = np.log(RH) * Rgaz * rho_h2o * (TairC + 273.15)/Mh2o * MPa2hPa * hPa2cm
    
    #air water potential
    r.psi_air = p_a #*MPa2hPa #used only with xylem
    return r

    
def setKrKx_phloem(r): 

    #number of vascular bundles
    VascBundle_leaf = 32
    VascBundle_stem = 52
    VascBundle_root = 1 #valid for all root type
            
    #number of sieve tubes per bundle
    numL = 18
    numS = 21
    numr0 = 33
    numr1 = 25
    numr2 = 25
    numr3 = 1 

    # axial conductivity [cm^3/day] , mu is added later as it evolves with CST  
    beta = 0.9 #Thompson 2003a
    kz_l   = VascBundle_leaf * numL* (0.00025 **4)   * np.pi /8 * beta  
    kz_s   = VascBundle_stem * numS *(0.00019 **4)    * np.pi /8 * beta
    kz_r0  = VascBundle_root *numr0 *(0.00039 **4)  * np.pi /8 * beta
    kz_r12 = VascBundle_root * numr1*(0.00035**4) * np.pi /8 * beta
    kz_r3  = VascBundle_root * numr3 *(0.00068**4)  * np.pi /8 * beta
    
    
    #radial conductivity [1/day], for carbon exudation
    kr_l  = 0. #always 0 for leaves
    kr_s  = 0. #always 0 for stem
    kr_r0 = 5e-2
    kr_r1 = 5e-2
    kr_r2 = 5e-2
    kr_r3 = 5e-2
    l_kr = 0.8 #cm zone from root tip at which exudation occures
    
    r.setKr_st([[kr_r0,kr_r1 ,kr_r2 ,kr_r0],[kr_s,kr_s ],[kr_l]])
    r.setKx_st([[kz_r0,kz_r12,kz_r12,kz_r0],[kz_s,kz_s ],[kz_l]])
    
    #cross-sectional area, to switch between length and volume
    a_ST = [[0.00039,0.00035,0.00035,0.00039 ],[0.00019,0.00019],[0.00025]]
    Across_s_l   = numL*VascBundle_leaf *(a_ST[2][0]**2)*np.pi
    Across_s_s   = numS *VascBundle_stem * (a_ST[1][0]**2)*np.pi     
    Across_s_r0  = numr0 *VascBundle_root * (a_ST[0][0]**2)*np.pi   
    Across_s_r12 = numr1*VascBundle_root * (a_ST[0][1]**2)*np.pi
    Across_s_r3  =  numr3 *VascBundle_root *(a_ST[0][2]**2)*np.pi   
    r.setAcross_st([[Across_s_r0,Across_s_r12,Across_s_r12,Across_s_r0],[Across_s_s,Across_s_s],[Across_s_l]])
    return(r)

def setPhotosynthesisParameters(r,weatherInit):
    r.g0 = 8e-3              #residual stomatal opening at night 
    r.VcmaxrefChl1 =1.28     #influence of leaf chlorophyl content on carboxylation rate
    r.VcmaxrefChl2 = 8.33    #influence of leaf chlorophyl content on carboxylation rate
    r.a1 = 0.5               #ci/(cs - ci)
    r.a3 = 1.5               # VcrefMax to VjrefMax ratio
    r.alpha = 0.4            #effect of light on photon flux rate
    r.theta = 0.6            #effect of light on photon flux rate
    r.cs = weatherInit["cs"] #external CO2 partial rpessure
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) #leaf chlorophyle content (mean value or defined per leaf segment)
    return r