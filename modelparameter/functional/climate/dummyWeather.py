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
