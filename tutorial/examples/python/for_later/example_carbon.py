import py_rootbox as rb
from rb_tools import *
import math

# accuracy = 0.1 # cm
# maxiter = 10
#
# # was added to the c++ code
# def elongate(rs, inc, dt, se):
# 
#     ol = np.sum(v2a(rs.getScalar(rb.ScalarType.length)))
#     i = 0
#         
#     rs_ = rb.RootSystem(rs) # copy
#     se.setScale(1.)
#     rs_.simulate(dt, True)
#     inc_ = np.sum(v2a(rs_.getScalar(rb.ScalarType.length))) - ol
#     
#     if inc_>inc and abs(inc_-inc)>accuracy: # check if we have to perform a binary search  
#         
#         sl = 0. # left           
#         sr = 1. # right
#                         
#         while abs(inc_-inc)>accuracy and i<maxiter: # binary search 
#                                     
#             m = (sl+sr)/2. # mid
#             rs_ = rb.RootSystem(rs) # copy        
#             se.setScale(m)  
#             rs_.simulate(dt, True) 
#             inc_ = np.sum(v2a(rs_.getScalar(rb.ScalarType.length))) - ol
#             print("\tsl, mid, sr ", sl, m, sr, inc_)
#             
#             if inc_>inc: # concatenate
#                 sr = m
#             else:
#                 sl = m                  
#             
#             i += 1            
#             
#         return rs_                        
#         
#     else:
#         return rs_
    
# Parameter
simtime = 300. # days
dt = 1
N = round(simtime/dt) # steps
maxinc = 20; # maximal length increment (cm/day), TODO base this value on some fancy model 

# Initialize root system
rs = rb.RootSystem()
name = "Zea_mays_4_Leitner_2014" 
rs.openFile(name) 

# Set up depth dependent elongation scaling function
scale_elongation = rb.EquidistantGrid1D(0,-50, 100) # for root elongation from 0 cm to -50 cm, 100 nodes = 99 layers
soil_strength = np.ones((99,))*0.5 # some data for the 99 layers, e.g. np.linspace(0.1, 1., 99)          
scales = np.exp(-0.4*soil_strength) # scales from some equation (scale = function(soil_strength) ), where scale in (0,1)
scale_elongation.data = a2v(scales) # set proportionality factors
  
# Proportionally scale this function
se = rb.ProportionalElongation()
se.setBaseLookUp(scale_elongation)
  
# Manually set scaling function 
for i in range(0,10):  
    p = rs.getRootTypeParameter(i+1)
    p.se = se

rs.initialize() 

ol = 0
 
# Simulation loop
for i in range(0,N):    
    
    print("\nSimulation step", i)

    # if maxinc is dynamic: set maxinc (cm/day) according to some model
    
    # if soil_strength is dynamic: update soil_strength according to some model (update like in L58-L60)
    
    rs.simulate(dt, maxinc, se, True) # True = disable debug messages, False = enable debug messages
        
    l = np.sum(v2a(rs.getScalar(rb.ScalarType.length)))
    inc =  l - ol
    ol = l
        
    print("elongated for", inc, " cm")    
        
rs.write("results/example_carbon.vtp")



# 
# print("root system copy test")
# 
# rs = rb.RootSystem()
# name = "Anagallis_femina_Leitner_2010" 
# rs.openFile(name)     
# rs.initialize() 
# rs.simulate(20) # for a bit
# 
# rs2 = rb.RootSystem(rs) # copy the root system
# 
# nodes = vv2a(rs.getNodes())
# nodes2 = vv2a(rs2.getNodes())
# print(nodes.shape, nodes2.shape)
# 
# nodes2 = vv2a(rs2.getNodes())
# 
# rs.simulate(10)
# rs2.simulate(10)
# 
# nodes = vv2a(rs.getNodes())
# nodes2 = vv2a(rs2.getNodes())
# print(nodes.shape, nodes2.shape)
# 
# uneq = np.sum(nodes!=nodes2)
# print("Unequal nodes", uneq)


