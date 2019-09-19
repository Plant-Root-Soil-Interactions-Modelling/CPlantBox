import math
import time as timer
import os # add the search path for py_rootbox.so (probably there is a nicer way to do it?)
import sys
cwd = os.getcwd()
i = cwd.index("CRootBox"+os.sep)
sys.path.append(cwd[0:i+8])

import numpy as np
import scipy.sparse.linalg as LA
from scipy import sparse

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import PSP_infiltration1D as inf

import py_rootbox as rb    
from rb_tools import *
import xylem_flux 
        
        
        
#
# look up function for radial xylem fluxes for soil matric potential
#    
def soil_p(x,y,z,psi,depth):
    z = min(z,0.)
    z = max(z,-depth+1.e-9) # z \in (-depth, 0]
    i = math.floor(float(-z/depth*(inf.n))) # i \in [0, n-1]
    return psi[i+1] # (i+1) \in [1,n], boundaries are psi[0] and psi[n+1]
 
 
             
#
# Initialize soil domain (from Soil Physics with Python)
#    
isSuccess, soil = inf.readSoil("soilUniform.txt") # clay.txt
if not isSuccess: 
    print("warning: wrong soil file.")
    quit()             
funcType = inf.CAMPBELL         
Se = 0.6 # initial degree of saturation ]0-1]
inf.initializeWater(funcType, soil, Se, inf.CELL_CENT_FIN_VOL) # for simplicity we use a linear grid
for i in range(inf.n+2):
    inf.oldTheta[i] = inf.theta[i]              
ubPotential = inf.waterPotential(funcType, soil[0],inf.thetaFromSe(funcType, soil[0], Se))    
print("Soil upper boundary potential ", ubPotential)
isFreeDrainage = True

sumInfiltration = 0
totalIterationNr = 0



#
# Initialize root domain
#
rsname = "Sorghum_bicolor_NA_NA" # "Lupinus_albus_Leitner_2014"
rs = rb.RootSystem()
rs.openFile(rsname,parameterPath())
for i in range(0,10): # set axial resolution   
    rs.getRootTypeParameter(i+1).dx = 0.1
rs.initialize() 
rs_Kr = np.array([ 2.e-10, 2.e-10, 2.e-10, 2.e-10, 2.e-10, 2.e-11, 2.e-11 ]) # s/m; root hydraulic radial conductivity per root type 
rs_Kz = np.array([ 5.e-14, 5.e-14, 5.e-14, 5.e-14, 5e-14, 5e-14, 5e-14 ]) # m²*s; root hydraulic axial conductivity per root type 

dirichlet = False



#
# Simulation parameters
#                   
simTime = 28*24*3600  #28*24*3600   
maxTimeStep = 3600           
dt = 36.              
time = 0              

out_delay = 2*3600 # figure output delay
out_next = out_delay             
out_img = np.zeros((inf.n+2,round(simTime/out_delay)))
out_c = 0
rs_out_next = 7*24*3600 # root system output times (vtp)

rho = 1e3 # kg / m^3      
g = 1.e-3*9.8065 # m / s^2   

pot_trans = np.array([2*-2e-9]) # # m^3 s^-1 potential transpiration [-5.e-8]
top_pot = -1500 # J kg^-1 top potential (wilting point)

ctflux = 0 # cumulative transpiration
  
plt.ion()
fig1 = plt.figure() # Prepare output 
plt.rc('xtick', labelsize=24) # bigger
plt.rc('ytick', labelsize=24) 


#
# Simulation loop
#
while (time < simTime):
    
    dt = min(dt, out_next - time)     
    
    #
    # 1. Richards Code (from Soil Physics with Python)
    #     
    t1 = timer.time()    
    ubPotential = inf.psi[2] # should be order 2 over time (centered around 1)
    success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage, inf.LOGARITHMIC)  
    while (not success):
        print ("dt =", dt, "no convergence")
        dt = max(dt/2., 1e-9) # go slower
        for i in range(inf.n+2): # reset to old time step
            inf.theta[i] = inf.oldTheta[i]
            inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
        ubPotential = inf.psi[2]
        success, nrIterations, flux = inf.cellCentFiniteVolWater(funcType, soil, dt, ubPotential, isFreeDrainage,inf.LOGARITHMIC)

    for i in range(inf.n+2):
        inf.oldTheta[i] = inf.theta[i]  
            
    sumInfiltration += flux * dt 
    time += dt    
    
    #
    # 2. Root System  
    #
    t2 = timer.time()
    rs.simulate(dt/3600/24) # seconds to days
  
    seg = seg2a(rs.getSegments())
    nodes = vv2a(rs.getNodes())/100. # convert to meter
    rs_ana = rb.SegmentAnalyser(rs) 
    type = v2a(rs_ana.getScalar(rb.ScalarType.type))
    radius = v2a(rs_ana.getScalar(rb.ScalarType.radius))/100. # convert to meter 
    time_ = v2a(rs_ana.getScalar(rb.ScalarType.time))*3600*24 # convert to seconds   
        
    #
    # 3. Root System fluxes
    #    
    t3 = timer.time()     
    kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
    kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))    
          
    soil_p2 = lambda x,y,z : soil_p(x,y,z, inf.psi,soil[-1].lowerDepth) # J/kg
                    
    if dirichlet == False: # Neumann
        Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p2) 
        Q, b = xylem_flux.bc_neumann(Q, b, np.array([0]), np.array([pot_trans]))
        x = LA.spsolve(Q, b) # direct
        eff_trans = xylem_flux.axial_flux0(x, seg, nodes, kz, rho, g) # verify that eff_trans == pot_trans        
        print("using neumann ( Effective Transpiration = " +str(eff_trans)+" m^3 s^-1)", "top potential", x[0] )        
         
    if dirichlet or (x[0]<top_pot):        
        Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p2) 
        dirichlet = True
        Q, b = xylem_flux.bc_dirichlet(Q, b, np.array([0]), np.array([top_pot])) 
        x = LA.spsolve(Q, b) # direct
        eff_trans = xylem_flux.axial_flux0(x, seg, nodes, kz, rho, g)
        dirichlet = eff_trans<pot_trans
        print("using dirichlet ( Effective Transpiration = "+str(eff_trans)+" m^3 s^-1)")
         
    radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p2) # used to calculate the sink term             
    
    #
    # 4. Sink term
    #    
    t4 = timer.time()
    sink = np.zeros(inf.n)
    
    for i in range(0,len(seg)): # calculate sink
        z1 = nodes[seg[i,0],2]
        z2 = nodes[seg[i,1],2]
        z = 0.5*(z1+z2)
        if z>-soil[-1].lowerDepth: # -1 meter
            ind = math.floor(-z/soil[-1].lowerDepth*inf.n)
            sink[ind] += ( 40 * radial_flux[i] * dt / 0.01 ) # per volume, i.e. per layer (40 plants per m^2)
    
    for i in range(0,inf.n): # apply sink
        inf.theta[i+1] += sink[i] 
        if inf.theta[i+1]>1:
            print("WARNING: cut max, layer ", i+1)
            inf.theta[i+1] = min(inf.theta[i+1],1)
        if inf.theta[i+1]<0:
            print("WARNING: cut min, layer ", i+1)
            inf.theta[i+1] = max(inf.theta[i+1],0)          
          
    for i in range(1,inf.n+1): # update psi
        inf.psi[i] = inf.waterPotential(funcType, soil[inf.hor[i]], inf.theta[i])
     
    ctflux += -float(np.sum(sink)) # cumulative flux  
     
    #
    # Output
    #
    t = timer.time()    
    t0 = t-t1    
    s0 = str(format(t0, '.3f'))
    s1 = str(format((t2-t1)/t0, '.3f'))
    s2 = str(format((t3-t2)/t0, '.3f'))
    s3 = str(format((t4-t3)/t0, '.3f'))
    s4 = str(format((t -t4)/t0, '.3f'))
    print("simtime =", format(time/3600./24., '.3f'), "days, spent ="+ s0 +"s ("+s1+", "+s2+", "+s3+", "+s4+")", rs.getNumberOfNodes(), "nodes,", "dt =", dt, 
          " iterations =", int(nrIterations), "cum trans", ctflux, "trans", -float(np.sum(0.01*sink/40))/dt, "summed infiltration:", format(sumInfiltration, '.3f'))
    
    if time>=rs_out_next: 
        rs_out_next += (7*24*3600)  
              
        segP = nodes2seg(nodes,seg,x)# save vtp 
        axial_flux = xylem_flux.axial_flux(x, seg, nodes, kz, rho, g)
        radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p2)
        net_flux = axial_flux+radial_flux
        rs_ana.addUserData(a2v(segP),"pressure")
        rs_ana.addUserData(a2v(axial_flux),"axial_flux")
        rs_ana.addUserData(a2v(radial_flux),"radial_flux")
        rs_ana.addUserData(a2v(net_flux),"net_flux")
        rs_ana.write(rsname+"_"+str(round(time/(24*3600)))+".vtp")                
              
        fig1a = plt.figure(figsize=(0.6*14.,0.6*12.)) # plot root length distribution          
        ax1 = fig1a.gca()
        ax2 = ax1.twiny()          
        rsl = v2a(rs_ana.distribution(rb.ScalarType.length,0.,100.,inf.n,False))
        ax2.set_ylim(-0.25, 0)
        ax2.set_xlim(0, 1.2)
        # ax2.set_xlabel("Root system length (m)",color='green',fontsize=32)
        # ax2.set_ylabel("Depth (m)",fontsize=32)                 
        l2 = ax2.plot(rsl/100.,-inf.z[1:len(inf.z)-1],color='green') # for greyscale...     
        plt.setp(l2, linewidth=2)     
#         xt_ = [0.,0.3,0.6,0.9,1.2,1.5]        
#         ax2.set_xticks([100*x for x in xt_] , [str(s) for s in xt_] )        
                
        ax1.set_ylim(-0.25, 0) # plot sink
        ax1.set_xlim(0.,2.75e-10)
        # ax1.set_xlabel("Sink (1)",color='blue',fontsize=32)
        # ax1.set_ylabel("Depth (m)",fontsize=32)         
        l1 = ax1.plot(-0.01*sink/40./dt,-inf.z[1:len(inf.z)-1],'-',color='blue')   # 0.01 m is the layer width, 40 plants
        plt.setp(l1, linewidth=2)
                               
        fig1a.show()                        
    
    if time>=out_next:
        out_next += out_delay
        
        out_img[:,out_c] =  inf.theta[:]/soil[0].thetaS # plot effective saturation            
        out_c += 1       
        fig1.clf()
        ax = fig1.gca()
        ax.set_xlabel("Time",fontsize=16)             
        ax.set_ylabel("Depth",fontsize=16)
        ax.set_title("Effective saturation",fontsize=16)     
        # ax2.set_xticks( )
        cax = ax.imshow(out_img[0:50,:], interpolation="nearest", cmap="hot")             
        fig1.colorbar(cax,ticks=np.linspace(0,Se,7))       
        # fig1.colorbar(cax)       
        plt.pause(0.0001) 
            
    if (float(nrIterations/inf.maxNrIterations) < 0.1): 
        dt = min(dt*2, maxTimeStep) # go faster        
           
print("finished without crashing!")
plt.ioff()
plt.show()
# fig1.show()







