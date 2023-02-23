import py_rootbox as rb
from rb_tools import *
import math
import matplotlib.pyplot as plt
from scipy import stats



def analyseTropism(N,sigma,dx):

    print("Analysis with N="+str(N)+", sigma="+str(sigma)+", dx="+str(dx),"\n")
    
    # Open plant and root parameter from a file
    rs = rb.RootSystem()
    name = "Triticum_aestivum_a_Bingham_2011" 
    rs.openFile(name) 

    for i in range(0,10):
        p = rs.getRootTypeParameter(i+1)
        p.tropismT = rb.TropismType.gravi
        p.tropismN = N
        p.tropismS = sigma
        p.dx = dx

    # Simulate
    rs.initialize() 
    rs.simulate(30) 

    # Analysis of angles
    poly = rs.getPolylines()

    e3 = np.array([0.,0.,-1.])
    angles = []
    dangles = []
    for line in poly:
        line_ = vv2a(line)    
        for i in range(1,line_.shape[0]-1):
        
            v1 = line_[i,:]-line_[i-1,:]
            v2 = line_[i+1,:]-line_[i,:]
            v1 = v1/np.linalg.norm(v1)
            v2 = v2/np.linalg.norm(v2)
        
            a1 = math.acos(float(np.inner(e3,v1)))/math.pi*180.
            angles.append(a1)
            a2 = math.acos(np.inner(e3,v2))/math.pi*180.
            dangles.append(a2-a1)
    return angles, dangles
            

sigma_ = np.array([5.,10.,15.,20.])/180*math.pi # angular change per cm
N_ = np.array([0.,1.,2.,5.,10.]) # number of trials 
dx = 0.5 # cm

for sigma in sigma_:
    for N in N_: 
        angles, dangles = analyseTropism(N,sigma,dx)
        slope, intercept, r_value, p_value, std_err = stats.linregress(angles, dangles)
        print("Slope ", slope)
        print("Intercept", intercept)      
        print("\n")     
           
# # Plot results (i do not see much at those plots)
# plt.title("Tropism behaviour N="+str(N)+", sigma="+str(sigma)+", dx="+str(dx))
# plt.xlabel("Angle between segment and ez")
# plt.ylabel("Angular change")
# plt.plot(angles, dangles, 'ro')
# # plt.savefig("results/example_tropisms.png")
# plt.show()
