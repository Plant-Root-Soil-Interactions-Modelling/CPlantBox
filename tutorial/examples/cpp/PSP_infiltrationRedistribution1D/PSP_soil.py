#PSP_soil.py
from PSP_readDataFile import *
from math import sqrt, log

g = 9.8065                 

NODATA = -9999.

CAMPBELL = 1
RESTRICTED_VG = 2
IPPISCH_VG = 3
VAN_GENUCHTEN = 4

CELL_CENT_FIN_VOL = 1
NEWTON_RAPHSON_MP = 2
NEWTON_RAPHSON_MFP = 3

LOGARITHMIC = 0
HARMONIC = 1
GEOMETRIC = 2

class Csoil:
    upperDepth = NODATA        
    lowerDepth = NODATA        
    Campbell_he = NODATA      
    Campbell_b = NODATA        
    CampbellMFP_he = NODATA   
    Campbell_b3 = NODATA 
    VG_alpha = NODATA         
    VG_n = NODATA             
    VG_m = NODATA
    VG_he = NODATA             
    VG_alpha_mod = NODATA      
    VG_n_mod = NODATA          
    VG_m_mod = NODATA          
    VG_Sc = NODATA             
    VG_thetaR = NODATA        
    Mualem_L = NODATA         
    thetaS = NODATA            
    Ks = NODATA                
        
def readSoil(soilFileName):
    mySoil = []
    A, isFileOk = readDataFile(soilFileName, 1, ',', False)
    if ((not isFileOk) or (len(A[0]) < 12)):
        print("error: wrong soil file.")
        return False, mySoil
    
    for i in range(len(A)):
        horizon = Csoil()
        horizon.upperDepth = A[i,0]
        horizon.lowerDepth = A[i,1]
        horizon.Campbell_he = A[i,2]
        horizon.Campbell_b = A[i,3]
        horizon.Campbell_n = 2.0 + (3.0 / horizon.Campbell_b)
        horizon.VG_he = A[i,4]
        horizon.VG_alpha = A[i,5]
        horizon.VG_n = A[i,6]
        horizon.VG_m =  1. - (1. / horizon.VG_n)
        horizon.VG_alpha_mod = A[i,7]
        horizon.VG_n_mod = A[i,8]
        horizon.VG_m_mod =  1. - (1. / horizon.VG_n_mod)
        horizon.VG_Sc =((1.+
        (horizon.VG_alpha_mod*abs(horizon.VG_he))**horizon.VG_n_mod)**(-horizon.VG_m_mod))
        horizon.VG_thetaR = A[i,9]
        horizon.thetaS = A[i,10]
        horizon.Ks = A[i,11]
        horizon.Mualem_L = 0.5
        horizon.CampbellMFP_he = horizon.Ks * horizon.Campbell_he / (1.0 - horizon.Campbell_n) 
        horizon.Campbell_b3 = (2.0 * horizon.Campbell_b + 3.0) / (horizon.Campbell_b + 3.0)
        mySoil.append(horizon)
    return True, mySoil

def getHorizonIndex(soil, depth):
    for index in range(len(soil)):
        if ((depth >= soil[index].upperDepth) 
        and (depth < soil[index].lowerDepth)): 
            return(index)
    lastHorizon = len(soil)-1
    if depth >= soil[lastHorizon].lowerDepth:
        return lastHorizon
    else:
        return(-1)

def airEntryPotential(funcType, soil): 
    if (funcType == CAMPBELL):
        return(soil.Campbell_he)
    elif (funcType == IPPISCH_VG):
        return(soil.VG_he)
    elif (funcType == RESTRICTED_VG):
        return(0)
    else:
        return(NODATA)
     
def waterPotential(funcType, soil, theta):
    psi = NODATA
    Se = SeFromTheta(funcType, soil, theta)
    if (funcType == RESTRICTED_VG):
        psi = -(1./soil.VG_alpha)*((1./Se)**(1./soil.VG_m) - 1.)**(1./soil.VG_n)
    elif (funcType == IPPISCH_VG):
        psi = -((1./soil.VG_alpha_mod)*
                ((1./(Se*soil.VG_Sc))**(1./soil.VG_m_mod)-1.)**(1./soil.VG_n_mod))
    elif (funcType == CAMPBELL):
        psi = soil.Campbell_he * Se**(-soil.Campbell_b)
    return(psi)
    
def SeFromTheta(funcType, soil, theta):
    if (theta >= soil.thetaS): return(1.)
    if (funcType == CAMPBELL):
        Se = theta / soil.thetaS
    else:
        Se = (theta - soil.VG_thetaR) / (soil.thetaS - soil.VG_thetaR)
    return (Se)

def thetaFromSe(funcType, soil, Se):
    if (funcType == RESTRICTED_VG) or (funcType == IPPISCH_VG):
        theta = (Se * (soil.thetaS - soil.VG_thetaR) + soil.VG_thetaR)
    elif (funcType == CAMPBELL):
        return(Se * soil.thetaS) 
    return(theta)

def degreeOfSaturation(funcType, soil, psi):
    if (psi >= 0.): return(1.)
    Se = NODATA
    if (funcType == IPPISCH_VG):
        if (psi >= soil.VG_he): Se = 1.
        else: 
            Se = (1./soil.VG_Sc) * pow(1.+pow(soil.VG_alpha_mod 
                                    * abs(psi), soil.VG_n_mod), -soil.VG_m_mod)
    elif (funcType == RESTRICTED_VG):
        Se = 1 / pow(1 + pow(soil.VG_alpha * abs(psi), soil.VG_n), soil.VG_m)  
    elif (funcType == CAMPBELL):
        if psi >= soil.Campbell_he: Se = 1.
        else: Se = pow(psi / soil.Campbell_he, -1. / soil.Campbell_b)
    return(Se) 

def thetaFromPsi(funcType, soil, psi):
    Se = degreeOfSaturation(funcType, soil, psi)
    theta = thetaFromSe(funcType, soil, Se)
    return(theta)
          
def hydraulicConductivityFromTheta(funcType, soil, theta): 
    k = NODATA      
    if (funcType == RESTRICTED_VG):
        Se = SeFromTheta(funcType, soil, theta)
        k = (soil.Ks * pow(Se, soil.Mualem_L) * 
             (1. -pow(1. -pow(Se, 1./soil.VG_m), soil.VG_m))**2)
    elif (funcType == IPPISCH_VG):
        Se = SeFromTheta(funcType, soil, theta)
        num   = 1. - pow(1. - pow(Se * soil.VG_Sc, 1./ soil.VG_m_mod), soil.VG_m_mod);
        denom = 1. - pow(1. - pow(soil.VG_Sc, 1./ soil.VG_m_mod), soil.VG_m_mod);
        k = soil.Ks * pow(Se, soil.Mualem_L) * pow((num / denom), 2.)
    elif (funcType == CAMPBELL):
        psi = waterPotential(funcType, soil, theta)
        k = soil.Ks * (soil.Campbell_he / psi)**soil.Campbell_n 
    return(k)

#---------------------------------------------
# dTheta/dH = dSe/dH (Theta_s - Theta_r)
#---------------------------------------------
def dTheta_dPsi(funcType, soil, psi):
    airEntry = airEntryPotential(funcType, soil)
    if (psi > airEntry): return 0.0
     
    if (funcType == RESTRICTED_VG):
        dSe_dpsi = soil.VG_alpha * soil.VG_n * (soil.VG_m 
                * pow(1. + pow(soil.VG_alpha * abs(psi), soil.VG_n), 
                -(soil.VG_m + 1.)) * pow(soil.VG_alpha * abs(psi), soil.VG_n - 1.))      
        return dSe_dpsi * (soil.thetaS - soil.VG_thetaR)
    elif (funcType == IPPISCH_VG):
        dSe_dpsi = soil.VG_alpha_mod * soil.VG_n_mod * (soil.VG_m_mod 
                * pow(1. + pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod), 
                -(soil.VG_m_mod + 1.)) * pow(soil.VG_alpha_mod * abs(psi), soil.VG_n_mod - 1.))      
        dSe_dpsi *= (1. / soil.VG_Sc)
        return dSe_dpsi * (soil.thetaS - soil.VG_thetaR)
    elif (funcType == CAMPBELL):
        theta = soil.thetaS * degreeOfSaturation(funcType, soil, psi) 
        return -theta / (soil.Campbell_b * psi)

def MFPFromTheta(soil, theta): 
    return (soil.CampbellMFP_he * (theta / soil.thetaS)**(soil.Campbell_b + 3.0)) 

def MFPFromPsi(soil, psi):
    return (soil.CampbellMFP_he * (psi / soil.Campbell_he)**(1.0 - soil.Campbell_n)) 

def thetaFromMFP(soil, MFP):
    if (MFP > soil.CampbellMFP_he):
        return(soil.thetaS) 
    else:
        return(soil.thetaS * (MFP / soil.CampbellMFP_he)**(1.0/(soil.Campbell_b + 3.0)))  

def hydraulicConductivityFromMFP(soil, MFP):
    k = soil.Ks * (MFP / soil.CampbellMFP_he)**soil.Campbell_b3
    return(k)
	
def dTheta_dH(funcType, soil, H0, H1, z): 
    psi0 = H0 + g*z
    psi1 = H1 + g*z
    if (abs(psi1-psi0) < 1E-5):
        return dTheta_dPsi(funcType, soil, psi0)
    else:
        theta0 = thetaFromPsi(funcType, soil, psi0)
        theta1 = thetaFromPsi(funcType, soil, psi1)
        return (theta1 - theta0) / (psi1 - psi0)

def meanK(meanType, k1, k2):
    if (meanType == LOGARITHMIC):
        if (k1 != k2):
            k = (k1-k2) / log(k1/k2)
        else:
            k = k1
    elif (meanType == HARMONIC): 
        k = 2.0 / (1.0 / k1 + 1.0 / k2)
    elif (meanType == GEOMETRIC): 
        k = sqrt(k1 * k2)
    return k