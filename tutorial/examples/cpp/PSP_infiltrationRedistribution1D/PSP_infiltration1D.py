#PSP_infiltration1D
from __future__ import division

import PSP_grid as grid
from PSP_ThomasAlgorithm import ThomasBoundaryCondition
from PSP_soil import *

waterDensity = 1000.        
area = 1                    
maxNrIterations = 100
tolerance = 1e-6
n = 100
 
hor = np.zeros(n+2, int)                      
z = np.zeros(n+2, float) 
zCenter = np.zeros(n+2, float) 
vol = np.zeros(n+2, float)     
a = np.zeros(n+2, float)       
b = np.zeros(n+2, float)       
c = np.zeros(n+2, float)       
d = np.zeros(n+2, float)       
dz = np.zeros(n+2, float)      
psi = np.zeros(n+2, float)     
dpsi = np.zeros(n+2, float)    
theta = np.zeros(n+2, float)   
oldTheta = np.zeros(n+2, float)  
C = np.zeros(n+2, float)       
k = np.zeros(n+2, float)       
u = np.zeros(n+2, float)       
du = np.zeros(n+2, float)      
f = np.zeros(n+2, float)       
H = np.zeros(n+2, float)       
H0 = np.zeros(n+2, float)      
        
def initializeWater(funcType, soil, se_0, solver):
    global z
    # vector depth [m]
    lastHorizon = len(soil)-1
    z = grid.linear(n, soil[lastHorizon].lowerDepth)    
    vol[0] = 0
    for i in range(n+1): 
        dz[i] = z[i+1]-z[i]
        if (i > 0): vol[i] = area * dz[i]
    for i in range(n+2): 
        zCenter[i] = z[i] + dz[i]*0.5
        
    if (solver == CELL_CENT_FIN_VOL):
        for i in range(n+1): 
            dz[i] = zCenter[i+1]-zCenter[i]
                  
    #initial conditions   
    psi[0] = 0
    for i in range(1, n+2):
        hor[i] = getHorizonIndex(soil, zCenter[i])
        theta[i] = thetaFromSe(funcType, soil[hor[i]], se_0)
        oldTheta[i] = theta[i]
        if (solver == NEWTON_RAPHSON_MFP):
            psi[i] = MFPFromTheta(soil[hor[i]], theta[i])
            k[i] = hydraulicConductivityFromMFP(soil[hor[i]], psi[i])
        else:
            psi[i] = waterPotential(funcType, soil[hor[i]], theta[i])
            k[i] = hydraulicConductivityFromTheta(funcType, soil[hor[i]], theta[i])
        H[i] = psi[i] - zCenter[i]*g
                    
def NewtonRapsonMP(funcType, soil, dt, ubPotential, isFreeDrainage):
    #apply upper boundary condition
    airEntry = airEntryPotential(funcType, soil[0])
    psi[1] = min(ubPotential, airEntry) 
    oldTheta[1] = thetaFromPsi(funcType, soil[0], psi[1])
    theta[1] = oldTheta[1]
    
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    
    nrIterations = 0
    massBalance = 1
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        massBalance = 0
        for i in range(1, n+1):
            k[i] = hydraulicConductivityFromTheta(funcType, soil[hor[i]], theta[i])
            u[i] = g * k[i]
            du[i] = -u[i] * soil[hor[i]].Campbell_n / psi[i]
            capacity = dTheta_dPsi(funcType, soil[hor[i]], psi[i])
            C[i] = (waterDensity * vol[i] * capacity) / dt
        
        for i in range (1, n+1):
            f[i] = ((psi[i+1] * k[i+1] - psi[i] * k[i]) 
                    / (dz[i] * (1 - soil[hor[i]].Campbell_n))) - u[i]
            if (i == 1): 
                a[i] = 0
                c[i] = 0
                b[i] = k[i] / dz[i] + C[i] + du[i]
                d[i] = 0
            else:
                a[i] = -k[i-1] / dz[i-1] - du[i-1]
                c[i] = -k[i+1] / dz[i]
                b[i] = k[i] / dz[i-1] + k[i] / dz[i] + C[i] + du[i]
                d[i] = f[i-1] - f[i] + (waterDensity * vol[i] 
                                        * (theta[i] - oldTheta[i])) /dt
                massBalance += abs(d[i])
    
        ThomasBoundaryCondition(a, b, c, d, dpsi, 1, n)
        
        for i in range(1, n+1):
            psi[i] -= dpsi[i]
            psi[i] = min(psi[i], airEntry)        
            theta[i] = thetaFromPsi(funcType, soil[hor[i]], psi[i])
        nrIterations += 1
        
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
    
    if (massBalance < tolerance):
        flux = -f[1]
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0
    
def cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage, meanType):
    #apply upper boundary condition
    airEntry = airEntryPotential(funcType, soil[0])
    psi[0] = min(ubPotential, airEntry)
    theta[0] = thetaFromPsi(funcType, soil[0], psi[0])
    theta[1] = theta[0]
    psi[1] = psi[0]
    
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        H[n+1] = psi[n+1] - zCenter[n+1]*g
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
        
    sum0 = 0
    for i in range(1, n+1):
        H0[i] = psi[i] - zCenter[i]*g
        H[i] = H0[i]
        sum0 += waterDensity * vol[i] * theta[i]
    
    massBalance = sum0
    nrIterations = 0
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        for i in range(1, n+1):
            k[i] = hydraulicConductivityFromTheta(funcType, soil[hor[i]], theta[i])
            capacity = dTheta_dH(funcType, soil[hor[i]], H0[i], H[i], zCenter[i])
            C[i] = (waterDensity * vol[i] * capacity) / dt
            
        f[0] = 0
        for i in range(1, n+1):
            f[i] = area * meanK(meanType, k[i], k[i+1]) / dz[i]
        
        for i in range(1, n+1):    
            a[i] = -f[i-1]
            if (i == 1):
                b[i] = 1
                c[i] = 0
                d[i] = H0[i]  
            elif (i < n):
                b[i] = C[i] + f[i-1] + f[i]
                c[i] = -f[i]
                d[i] = C[i] * H0[i]
            elif (i == n):
                b[n] = C[n] + f[n-1]
                c[n] = 0
                if (isFreeDrainage):
                    d[n] = C[n] * H0[n] - area * k[n] * g
                else:
                    d[n] = C[n] * H0[n] - f[n]* (H[n]-H[n+1])
              
        ThomasBoundaryCondition(a, b, c, d, H, 1, n)
        
        newSum = 0
        for i in range(1, n+1):
            psi[i] = H[i] + g*zCenter[i]
            theta[i] = thetaFromPsi(funcType, soil[hor[i]], psi[i])
            newSum += waterDensity * vol[i] * theta[i]
            
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]

        if (isFreeDrainage):    
            massBalance = abs(newSum - (sum0 + f[1]*(H[1]-H[2])*dt 
                                        - area*k[n]*g*dt))
        else:
            massBalance = abs(newSum - (sum0 + f[1]*(H[1]-H[2])*dt 
                                        - f[n]*(H[n]-H[n+1])*dt))
        nrIterations += 1
        
    if (massBalance < tolerance):
        flux = f[1]*(H[1]-H[2])
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0
    
# Infiltration simulation using Matric Flux Potential
def NewtonRapsonMFP(funcType, soil, dt, ubPotential, isFreeDrainage):
    #apply upper boundary condition
    airEntry = airEntryPotential(funcType, soil[0])
    ubPotential = min(ubPotential, airEntry) 
    psi[1] = MFPFromPsi(soil[0], ubPotential)
    oldTheta[1] = thetaFromMFP(soil[0], psi[1])
    theta[1] = oldTheta[1]
    k[1] = hydraulicConductivityFromMFP(soil[0], psi[1])
    psi[0] = psi[1]
    k[0] = 0.0
    
    if (isFreeDrainage):
        psi[n+1] = psi[n]
        theta[n+1] = theta[n] 
        k[n+1] = k[n]
    
    nrIterations = 0
    massBalance = 1
    while ((massBalance > tolerance) and (nrIterations < maxNrIterations)):
        massBalance = 0
        for i in range(1, n+1):
            k[i] = hydraulicConductivityFromMFP(soil[hor[i]], psi[i])
            capacity = theta[i] / ((soil[hor[i]].Campbell_b + 3.0) * psi[i])
            C[i] = waterDensity * vol[i] * capacity / dt
            u[i] = g * k[i]
            f[i] = (psi[i+1] - psi[i]) / dz[i] - u[i]
            if (i == 1): 
                a[i] = 0
                c[i] = 0
                b[i] = 1.0/dz[i] + C[i] + g * soil[hor[i]].Campbell_b3 * k[i] / psi[i]
                d[i] = 0
            else:
                a[i] = -1.0/dz[i-1] -g * soil[hor[i-1]].Campbell_b3 * k[i-1] / psi[i-1]
                c[i] = -1.0/dz[i]
                b[i] = 1.0/dz[i-1] + 1.0/dz[i] + C[i] + g * soil[hor[i]].Campbell_b3 * k[i] / psi[i]
                d[i] = f[i-1] - f[i] + (waterDensity * vol[i] 
                                        * (theta[i] - oldTheta[i]) / dt)
                massBalance += abs(d[i])
    
        ThomasBoundaryCondition(a, b, c, d, dpsi, 1, n)
        
        for i in range(1, n+1):
            psi[i] -= dpsi[i]
            psi[i] = min(psi[i], soil[hor[i]].CampbellMFP_he)
            theta[i] = thetaFromMFP(soil[hor[i]], psi[i])
            
        if (isFreeDrainage):
            psi[n+1] = psi[n]
            theta[n+1] = theta[n] 
            k[n+1] = k[n]
            
        nrIterations += 1
        
    if (massBalance < tolerance):
        flux = -f[1]
        return True, nrIterations, flux
    else:
        return False, nrIterations, 0