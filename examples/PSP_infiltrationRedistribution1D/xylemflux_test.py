import math
import time
import os # add the search path for py_rootbox.so (probably there is a nicer way to do it?)
import sys
cwd = os.getcwd()
i = cwd.index("CRootBox"+os.sep)
sys.path.append(cwd[0:i+8]) 

import numpy as np
from scipy import sparse
import scipy.sparse.linalg as LA

import matplotlib.pyplot as plt

import py_rootbox as rb    
from rb_tools import *

import xylem_flux 

#
# initialize (root system)
#
rs = rb.RootSystem()
rsname = "Lupinus_albus_Leitner_2014" 

# path = parameterPath()
# rs.openFile(rsname,path) 

p0 = rb.RootTypeParameter()
p1 = rb.RootTypeParameter()
# Taproot
p0.name = "taproot"
p0.type = 1
p0.lb = 1
p0.la = 10
p0.nob = 20
p0.ln = 89./19.
p0.r = 1
p0.dx = 0.5
# p0.k = maxRootLength(p0.la,p0.lb,p0.ln,p0.nob)
# print(p0)
# 1st order lateral
p1.name = "lateral"
p1.type = 2
p1.la = 25
p1.ln = 0
p1.r = 2
p1.k = 25
p1.dx = 0.1
rs.setRootTypeParameter(p0)
rs.setRootTypeParameter(p1)

rs.initialize() # hydrotropism is not set right now, link to soil is missing


#
# simulate
#
rs.simulate(30)

#
# results 
#
seg = seg2a(rs.getSegments())
nodes = vv2a(rs.getNodes())/100. # convert to meter
rs_ana = rb.SegmentAnalyser(rs) 
type = v2a(rs_ana.getScalar(rb.ScalarType.type))
radius = v2a(rs_ana.getScalar(rb.ScalarType.radius))/100. # convert to meter 
time_ = v2a(rs_ana.getScalar(rb.ScalarType.time))*3600*24 # convert to seconds 

#
# initialize (xylem_flux)
#
rs_Kr = np.array([ 1.16e-6, 1.74e-5, 1.74e-5, 1.74e-5, 11.74e-5, 1.74e-5, 1.74e-5 ]) # s/m; root hydraulic radial conductivity per root type 
rs_Kz = np.array([ 2.3e-8, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11, 1.16e-11 ]) # m²*s; root hydraulic axial conductivity per root type  
kr = np.array(list(map(lambda t: rs_Kr[int(t)-1], type))) # convert from 'per type' to 'per segment'
kz = np.array(list(map(lambda t: rs_Kz[int(t)-1], type)))
# kr = kr * (10*time_ + 1)
# kz = kz / (10*time_ + 1)
# kr = kr / 100
# kz = kz * 100

rho = 1e3 # kg / m^3      
g = 9.8  # m / s^2   

soil_p = lambda x,y,z : -100 # 

#
# create linear system
#
Q, b = xylem_flux.linear_system(seg, nodes, radius, kr, kz, rho, g, soil_p)  

#
# apply BC
#
n0= np.array([0]) # node indices
seg0 = np.array([0]) # segment indices
potT = np.array([-5.79e-4]) # potential Transpiration
topP = np.array([-1500]) # top potential (J/kg)

# Q, b = xylem_flux.bc_neumann(Q, b, seg0, potT, seg, nodes) # Neumann
Q, b = xylem_flux.bc_dirichlet(Q, b, n0, topP) # Dirichlet

#
# solve LS
#
t = time.time()
print("solve")
# x0 = d*np.ones(Q.shape[0]) # empirically proofen to be the best
# x, info = LA.cg(Q,b,x0=x0, tol=1e-12)  # tested with CG, CGS, GMRES, BICG. CG by far the best
x = LA.spsolve(Q, b) # direct
print("fin")
# print("CG: " + str(time.time()-t) +" sec" ) 
print("spsolve: " + str(time.time()-t) +" sec" )

#
# output
# 
f0 = xylem_flux.axial_flux0(x, seg, nodes, kz, rho, g)
print("Effective transpiration.: ", f0, " [?]")


segP = nodes2seg(nodes,seg,x) 
axial_flux = xylem_flux.axial_flux(x, seg, nodes, kz, rho, g)
radial_flux = xylem_flux.radial_flux(x, seg, nodes, radius, kr, soil_p)
net_flux = axial_flux+radial_flux

rs_ana.addUserData(a2v(segP),"pressure")
rs_ana.addUserData(a2v(axial_flux),"axial_flux")
rs_ana.addUserData(a2v(radial_flux),"radial_flux")
rs_ana.addUserData(a2v(net_flux),"net_flux")
rs_ana.write(rsname+".vtp")




