import numpy as np


def setPhotosynthesisParameters(r,weatherInit):
    r.g0 = 8e-3              #residual stomatal opening at night 
    r.VcmaxrefChl1 =1.28     #influence of leaf chlorophyl content on carboxylation rate
    r.VcmaxrefChl2 = 8.33    #influence of leaf chlorophyl content on carboxylation rate
    r.a1 = 0.5               #ci/(cs - ci)
    r.a3 = 1.5               # VcrefMax to VjrefMax ratio
    r.alpha = 0.4            #effect of light on photon flux rate
    r.theta = 0.6            #effect of light on photon flux rate
    r.pCO2 = weatherInit["cs"] #external CO2 partial rpessure
    SPAD= 41.0
    chl_ = (0.114 *(SPAD**2)+ 7.39 *SPAD+ 10.6) # Quian 2021, doi: 10.1029/2020JG006076
    r.Chl = np.array( [chl_]) #leaf chlorophyle content (mean value or defined per leaf segment)
    return r