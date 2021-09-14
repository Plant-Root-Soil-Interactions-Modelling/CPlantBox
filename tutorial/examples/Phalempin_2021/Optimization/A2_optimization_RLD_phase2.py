""" updated soil core example """

import sys
sys.path.append("../../../..")
import plantbox as pb
import vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

path = "rootsystem_params/"
name = "Optimization_loam"

def get_result(rs, time :float):
    """ 
    Retrieves a the state of the root systems at a certain time 
    in a single SegmentAnalyser object
    
    @param allRS     list of root systems 
    @param time      of the simulation result (days)
    @return SegmentAnalyser object conaining the segments of 
    all root sytstems at the specific time 
    """
    a = pb.SegmentAnalyser(rs)
    a.filter("creationTime", 0., time)
    a.pack()
    return a


def soil_cores(x :list, y :list, r :float,  h :list, up: list):
    """
    A lsit of soil core geometries with a fixed location in the field  
 
    @param x     x coordinates of the soil cores (cm)
    @param y     y coordinates of the soil cores (cm)
    @param r     radius of the soil core (cm)
    @param h     height of the soil core (cm)
    """
    assert len(x) == len(y), "coordinate length must be equal"
    core = []
    for i in range(0,len(h)): 
        core.append(pb.SDF_PlantContainer(r, r, h[i], False))
    cores = []
    for i in range(0, len(x)):
        cores.append(pb.SDF_RotateTranslate(core[i], 0., pb.SDF_Axis.xaxis, pb.Vector3d(x[i], y[i], up[i])))  # just translate
    return cores;

def set_all_sd(rs):
    for p in rs.getRootRandomParameter():
        p.las = 0
        p.lbs = 0
        p.rs = 0
        p.lmaxs = 0
        p.thetas = 0

def err(fitparams):
    r = fitparams[0]; ln1=fitparams[1]; ln1s=fitparams[2]; ln2 = fitparams[3];
    ln2s = fitparams[4]; maxB=fitparams[5]; maxBs=fitparams[6]; 
    
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap root type
    p1.r = r;
    p1.ln = ln1;
    p1.lns = ln1s;
    rs.setOrganRandomParameter(p1)
    
    p4 = rs.getRootRandomParameter(4)  # basal root type
    p4.r = r;
    p4.ln = ln1;
    p4.lns = ln1s;
    rs.setOrganRandomParameter(p4)
    
    p2 = rs.getRootRandomParameter(2)
    p2.ln = ln2
    p2.lns = ln2s
    rs.setOrganRandomParameter(p2)
    
    srp = rs.getRootSystemParameter()
    srp.maxB = maxB
    srp.maxBs = maxBs
    rs.setRootSystemParameter(srp)

    #set the standard deviation of all parameters zero to facilitate optimization
    set_all_sd(rs)

    #define the cores 
    A = np.zeros((3,4))
    A = [[-2.5,   0,  -3.5,  -6.5],
        [-2.5,    0,  -8.5,  -11.5],
        [-2.5,    0,  -13.5, -16.5]]
    A_ = np.array(A)

    x = A_[:,0]
    y = A_[:,1]
    h = np.subtract(A_[:,2],A_[:,3])
    up = A_[:,2]
    dow = A_[:,3]

    cores = soil_cores(x, y, r1, h, up)
    
    #make simulations
    rs.setGeometry(bigcyl) #set the geometry of the big cylinder

    for i in range(0,sims):
        rs.setSeed(i+1) 
        rs.initialize()
        rs.simulate(times)

        for j in range(0, len(cores)):
            core_analyser = get_result(rs, times)
            core_analyser.crop(cores[j]);
            core_analyser.pack()
            tl1 = core_analyser.distribution("length", up[j], dow[j], 1, True)  # vertical length distribution
            tl1 = np.array(tl1) / ( r1 * r1 * math.pi * h[j])  #RLD within the ind cores in cm/cm³
            rm_[j,i] = tl1 #matrix that contains the root length density

    rm = np.mean(rm_,axis = 1)
    std = np.std(rm_,axis = 1)
    rmstd = np.concatenate((rm, std))
    np.savetxt("RLDs/p2_virt_rm_std_"+mat+"_"+den+".txt",rmstd,fmt='%.2f')

    #compare measured and simulated core RLDs 
    with open("RLDs/p2_measured_RLD_std_"+mat+"_"+den+".txt") as f: #read in measured RLDs
        real_RLD_ = [[float(xx) for xx in line.split()] for line in f]

    #flatten the vectors
    RLD = np.reshape(rmstd, (corenum*2,1))
    real_RLD = np.reshape(np.array(real_RLD_), (corenum*2,1))

    NRMSE = math.sqrt(sum(((np.subtract(RLD,real_RLD)))**2)/len(RLD))
    err = NRMSE
    print('error= ', err)
    return err

######################################################################################
# Set geometry of soil cylinder,  top radius 5 cm, bot radius 5 cm, height 18 cm
bigcyl = pb.SDF_PlantContainer(5, 5, 18, False)

#Parameters to be fitted and initial values
r = 2.6 #cm/d
ln1 = 0.6
ln2 = 0.6 #cm
maxB = 15
ln1std = 0.05
ln2std = 0.05 #cm
maxBstd = 10

p0=[r,ln1,ln1std, ln2, ln2std,maxB,maxBstd] #initial guess
mat = 'loam' #or sand
den = 'low' #or high
sims = 20

corenum = 3 #number of soil cores taken 
times = 21 #days after planting
r1 = 1.5  # core radius
exportVTP = False  # export single root system cores (for debugging)
rm_ = np.zeros((corenum,sims)) #preallocation of result matrices
rm = np.zeros((corenum)) #preallocation of result matrices
std = np.zeros((corenum)) #preallocation of result matrices
rmstd = np.zeros((corenum*2)) #preallocation of result matrices


res = scipy.optimize.minimize(err, p0, method='L-BFGS-B')

print("initial guess= ", p0)
print("optimized params = ", res.x)

#save the parameters of the final root system 
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
p1 = rs.getRootRandomParameter(1)
p1.r = res.x[0]
p1.ln = res.x[1]
p1.lns = res.x[2]
rs.setOrganRandomParameter(p1)

p4 = rs.getRootRandomParameter(4)
p4.r = res.x[0]
p4.ln = res.x[1]
p4.lns = res.x[2]
rs.setOrganRandomParameter(p4)

p2 = rs.getRootRandomParameter(2)
p2.ln = res.x[3]
p2.lns = res.x[4]
rs.setOrganRandomParameter(p2)

srp = rs.getRootSystemParameter()
srp.maxB = res.x[5]
srp.maxBs = res.x[6]
rs.setRootSystemParameter(srp)

set_all_sd(rs) #set the std zero
rs.writeParameters("rootsystem_params/p2_Optimization_"+mat+".xml")








