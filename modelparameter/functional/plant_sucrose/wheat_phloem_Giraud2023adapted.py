import numpy as np

    
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
