import math
import numpy as np
from scipy import sparse
from numpy.linalg.linalg import norm

#
# Creates the linear system describing the pressure inside a xylem network
#
# in:
# seg     numpy array (Ns,2) of segment indices [1]
# nodes   numpy array (N,3) of the node coordinates [L]
# radius  segment radii [L]
# kr      radial conductivity for each segment [L2 T M−1]
# kz      axial conductivity for each segment [L5 T M-1]
# rho     density of soil water [M L-3]
# g       gravitational acceleration [L T−2]
# soil_p  lambda function returning the soil matric potential at a given location, p=soil_p(x,y,z) [M L−1 T−2]
#
# out: 
# Q,b     The equations are represented by the linear system Qx=b
# 
def linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p):

    Ns = seg.shape[0]
    N = nodes.shape[0]
    
    I = np.zeros(4*Ns)
    J = np.zeros(4*Ns)    
    V = np.zeros(4*Ns)
    b = np.zeros(N)
    
    k = 0     

    for c in range(0,Ns):
        
        i = seg[c,0]
        j = seg[c,1]
        
        n1 = nodes[i,:]
        n2 = nodes[j,:]
       
        mid = 0.5*(n1+n2)        
        p_s = soil_p(mid[0],mid[1],mid[2]) # evaluate soil matric potential

        v = n2-n1
        l = norm(v)        
        vz = v[2] / l # normed direction        
        
        a = radius[c]
        
        cii = a*math.pi*l*kr[c]/2 + kz[c]/l # Eqn (10)
        cij = a*math.pi*l*kr[c]/2 - kz[c]/l # Eqn (11)
        bi = a*math.pi*l*kr[c]*p_s # first term of Eqn (12) & (13)            
        
        # edge ij
        b[i] +=  ( bi - kz[c]*rho*g*vz )  # Eqn (12)     
        
        I[k] = i
        J[k] = i     
        V[k] += cii
        k += 1                
        
        I[k] = i
        J[k] = j        
        V[k] += cij
        k += 1 
        
        # edge ji
        i, j = j, i
        b[i] += ( bi + kz[c]*rho*g*vz ) # Eqn (13) 
  
        I[k] = i
        J[k] = i  
        V[k] += cii    
        k += 1                
          
        I[k] = i
        J[k] = j        
        V[k] += cij        
        k += 1 
         
    Q = sparse.coo_matrix((V,(I,J)))    
    Q = sparse.csr_matrix(Q) # Sparse row matrix seems the most reasonable to solve Qx = b iteratively
    
    return (Q, b)

#
# Modifies the linear system to describe Diriclet BC at the node indices n0
#
# in:
# Q, b    the linear system
# n0      node indices where to apply the Dirichlet BC
# d       fixed potential at n0, i.e. len(d)==len(n0)
#
# out:
# Q, b    the updated linear system
#
def bc_dirichlet(Q, b, n0, d):
    c = 0
    for c in range(0, len(n0)):
        i = n0[c]# print("Dirichlet BC at node "+str(i))             
        e0 = np.zeros((1,Q.shape[1])) # build zero vector
        Q[i,:] = sparse.csr_matrix(e0) # replace row i with ei
        Q[i,i] = 1
        b[i] = d[c]    

    return Q, b 

#
# Modifies the linear system to describe a Neumann BC at the segments seg0
#
# in:
# Q, b    the linear system
# n0      node indices where to apply the Dirichlet BC
# f       flux at n0, i.e. len(d)==len(n0)
#
# out:
# Q, b    the updated linear system
#
def bc_neumann(Q, b, n0, f):
    c = 0
    for c in range(0, len(n0)):                
        i = n0[c]  # print("Neumann BC at node "+str(i))       
        b[i] += f[c]        

    return Q, b 

#
# Calculates the axial flux for each segment
#
# in: 
# p       xylem pressure (i.e. solution vector)
# seg     numpy array (Ns,2) of segment indices [1]
# nodes   numpy array (N,3) of the node coordinates [L]
# kz      axial conductivity for each segment [L5 T]
# rho     density of soil water [M L-3]
# g       gravitational acceleration [L T−2]
#
# out:
# the axial flux
def axial_flux(p, seg, nodes, kz, rho, g):
    af = np.zeros(seg.shape[0])
    c = 0
    for s in seg:
        i = s[0]
        j = s[1]        
        v = nodes[j,:]-nodes[i,:] # segment direction
        l = norm(v) # length        
        v = v / l # normed direction           
        af[c] = -kz[c]*((p[j]-p[i])/l+rho*g*v[2])  # Eqn (6)
        c += 1
        
    return af

#
# Calculates the axial flux for the top segment
#
# in:
# p       xylem pressure (i.e. solution vector) 
# seg     numpy array (Ns,2) of segment indices [1]
# nodes   numpy array (N,3) of the node coordinates [L]1]
# kz      axial conductivity for each segment [L5 T]
# rho     density of soil water [M L-3]
# g       gravitational acceleration [L T−2]
#
# out:
# the axial flux
def axial_flux0(p, seg, nodes, kz, rho, g):
    s = seg[0]
    i = s[0]
    j = s[1]        
    v = nodes[j,:]-nodes[i,:] # segment direction
    l = norm(v) # length        
    v = v / l # normed direction           
    af = -kz[0]*((p[j]-p[i])/l+rho*g*v[2])  # Eqn (6)
    return af

#
# Calculates the radial flux for each segment
#
# in:
# p       xylem pressure (i.e. solution vector) 
# seg     numpy array (Ns,2) of segment indices [1]
# nodes   numpy array (N,3) of the node coordinates [L]
# radius  segment radii [L]
# kr      radial conductivity for each segment [L2 T M−1]
# soil_p  lambda funciton returning the soil matric potential at a given location, p=soil_p(x,y,z) [M L−1 T−2]
#
# out:
# the radial flux [L3 T-1] its radial flow
def radial_flux(p, seg, nodes, radius, kr, soil_p):
    rf = np.zeros(seg.shape[0])
    c = 0
    for s in seg:
        i = s[0]
        j = s[1] 
        n1 = nodes[i,:]
        n2 = nodes[j,:]       
        l = norm(n2-n1) # length
        a = radius[c]
        mid = (n1+n2)/2 # segment mid point 
        ps = soil_p(mid[0],mid[1],mid[2])
        rf[c] = -2*a*math.pi*l*kr[c]*(ps-(p[j]+p[i])/2) # Eqn (7)
        c += 1
                
    return rf

#
# Calculates the radial net for each segment
#
# in:
# p       xylem pressure (i.e. solution vector) 
# seg     numpy array (Ns,2) of segment indices [1]
# nodes   numpy array (N,3) of the node coordinates [L]
# radius  segment radii [L]
# kr      radial conductivity for each segment [L2 T M−1]
# kz      axial conductivity for each segment [L5 T]
# rho     density of soil water [M L-3]
# g       gravitational acceleration [L T−2]
# soil_p  lambda funciton returning the soil matric potential at a given location, p=soil_p(x,y,z) [M L−1 T−2]
#
# out:
# the net flux
def net_flux(p,seg, nodes, radius, kr, kz, rho, g, soil_p):
    return axial_flux(p, seg, nodes, radius, kr, kz, rho, g, soil_p) + radial_flux(p,seg, nodes, radius, kr, kz, rho, g, soil_p)


