import numpy as np

def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

    
def qair2rh(qair, press,TK):
    #pressPa = press *100
    #eaPa = qair* pressPa/(0.622+0.378*qair)
    #ea = eaPa/100
    T0 = 273.16
    RH = 26.3 * press * qair  /(exp((17.67*(TK - T0))/(TK- 29.65)))
    return RH
    
def weather( hp:float=1):
    Qnigh = 0; Qday = 960e-6 #458*2.1
    loam = [0.08, 0.43, 0.04, 1.6, 50]

    Tnigh = 15.8; Tday = 22
    #Tnigh = 13; Tday = 20.7
    #specificHumidity = 0.0097
    RHday = 0.6; RHnigh = 0.88
    Pair = 1010.00 #hPa
    #thetaInit = 0.4#
    pmean = -100.
    cs = 350e-6

    coefhours = 1
    RH_ = RHnigh + (RHday - RHnigh) * coefhours
    TairC_ = Tnigh + (Tday - Tnigh) * coefhours
    Q_ = Qnigh + (Qday - Qnigh) * coefhours
    
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    ea = es*RH_#qair2ea(specificHumidity,  Pair)
    assert ea < es
    #RH = ea/es
    assert ((RH_ > 0) and(RH_ < 1))
    bl_thickness = 1/1000 #1mm * m_per_mm
    diffusivity= 2.5e-5#m2/sfor 25*C
    rbl =bl_thickness/diffusivity #s/m 13
    #cs = 350e-6
    Kcanopymean = 1e-1 # m2/s
    meanCanopyL = (2/3) * hp /2
    rcanopy = meanCanopyL/Kcanopymean
    windSpeed = 2 #m/s
    zmzh = 2 #m
    karman = 0.41 #[-]

    rair = 1
    if hp > 0:
        rair = np.log((zmzh - (2/3)*hp)/(0.123*hp)) * np.log((zmzh - (2/3)*hp)/(0.1*hp)) / (karman*karman*windSpeed)
        #print()
        #raise Exception


    #pmean = min(theta2H(loam, thetaInit),-100.)

    weatherVar = {'TairC' : TairC_,'TairK' : TairC_ + 273.15,'Pair':Pair,"es":es,
                    'Qlight': Q_,'rbl':rbl,'rcanopy':rcanopy,'rair':rair,"ea":ea,
                    'cs':cs, 'RH':RH_, 'p_mean':pmean#,'theta':thetaInit
                 }
    #print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar   

def resistance2conductance(resistance,r, weatherX):
    resistance = resistance* (1/100) #[s/m] * [m/cm] = [s/cm]
    resistance = resistance * r.R_ph * weatherX["TairK"] / r.Patm # [s/cm] * [K] * [hPa cm3 K-1 mmol-1] * [hPa] = [s] * [cm2 mmol-1]
    resistance = resistance * (1000) * (1/10000)# [s cm2 mmol-1] * [mmol/mol] * [m2/cm2] = [s m2 mol-1]
    return 1/resistance

def init_conductivities(r, TairC:float = 20):
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
    betaXylX = 1#0.1      
    kz_l  = VascBundle_leaf *(rad_x_l_1 + rad_x_l_2)    *np.pi /(mu * 8)  * betaXylX
    kz_s  = VascBundle_stem *(rad_x_s_1 + rad_x_s_2)    *np.pi /(mu * 8) * betaXylX 
    kz_r0 = VascBundle_root * rad_x_r0_1                *np.pi /(mu * 8) * betaXylX  
    kz_r1 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX
    kz_r2 = VascBundle_root *(rad_x_r12_1 + rad_x_r12_2)*np.pi /(mu * 8)  * betaXylX 
    kz_r3 = VascBundle_root * rad_x_r3_1                *np.pi /(mu * 8) * betaXylX

    #radial conductivity [1/day],0.00014 #
    betaXyl = 1#0.1#0.1
    kr_l  = 3.83e-5 * hPa2cm * betaXyl# init: 3.83e-4 cm/d/hPa
    kr_s  = 0.#1.e-20  * hPa2cm # set to almost 0
    kr_r0 =6.37e-5 * hPa2cm * betaXyl
    kr_r1 =7.9e-5  * hPa2cm * betaXyl
    kr_r2 =7.9e-5  * hPa2cm * betaXyl
    kr_r3 =6.8e-5  * hPa2cm * betaXyl
    
    ratio_decrease = 1.5/100
    
    r.setKrTables([[[kr_r0,kr_r0,kr_r0*ratio_decrease,kr_r0*ratio_decrease],
                    [kr_r1,kr_r1,kr_r1*ratio_decrease,kr_r1*ratio_decrease],
                    [kr_r2,kr_r2,kr_r2*ratio_decrease,kr_r2*ratio_decrease],
                    [kr_r0,kr_r0,kr_r0*ratio_decrease,kr_r0*ratio_decrease]],
                    [[kr_s],[kr_s],[kr_s] ],[[kr_l]]],
            [[[0,0.8,1,10000],[0,0.8,1,10000],[0,0.8,1,10000],[0,0.8,1,10000]],
            [[0],[0],[0]],[[0]]],verbose = False, ageBased = False) 
            
    r.setKx([[kz_r0,kz_r1,kz_r2,kz_r0],[kz_s,kz_s,kz_s ],[kz_l]])
    
    
    Rgaz=8.314 #J K-1 mol-1 = cm^3*MPa/K/mol
    rho_h2o = dEauPure/1000#g/cm3
    Mh2o = 18.05 #g/mol
    MPa2hPa = 10000
    hPa2cm = 1/0.9806806
    return r

def photosynthesisParam(r):
    
    r.weatherX = weather()
    r = init_conductivities(r,r.weatherX['TairC'])
    r.g0 = 8e-6
    r.VcmaxrefChl1 =1.28#/2
    r.VcmaxrefChl2 = 8.33#/2
    r.a1 = 0.6/0.4#0.7/0.3#0.6/0.4 #ci/(cs - ci) for ci = 0.6*cs
    r.a3 = 1.5
    r.alpha = 0.4#0.2#/2
    r.theta = 0.6#0.9#/2
    r.k_meso = 1e-3#1e-4
    
    r.fwr = 1e-16
    r.fw_cutoff =   0.04072 # 0.09497583
    r.sh =5e-4# 4.655e-4# 4e-4
    r.gm=0.01
    r.p_lcrit = -8375#-8847# from calibration using data of corso2020 #-15000*0.6
    
    r.limMaxErr = 1/100
    r.maxLoop = 10000
    r.minLoop=900
    r.cs = r.weatherX["cs"]

    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6)/10
    r.Chl = np.array( [chl_]) 
    

    r.Patm = r.weatherX["Pair"]

    ##resistances
    r.g_bl = resistance2conductance(r.weatherX["rbl"],r, r.weatherX) / r.a2_bl
    r.g_canopy = resistance2conductance(r.weatherX["rcanopy"],r, r.weatherX) / r.a2_canopy
    r.g_air = resistance2conductance(r.weatherX["rair"],r, r.weatherX) / r.a2_air
    # intercepted PAR
    Qlightv = np.linspace(r.weatherX["Qlight"]/100, r.weatherX["Qlight"],num = len(r.get_segments_index(4)))
    r.vQlight = Qlightv # 
    r.Qlight = r.weatherX["Qlight"]

    
    return r