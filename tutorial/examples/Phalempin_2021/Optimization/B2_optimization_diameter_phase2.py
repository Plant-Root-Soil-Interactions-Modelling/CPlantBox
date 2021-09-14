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
name = "p2_Optimization_loam"

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
    a1 = fitparams[0]; a2=fitparams[1]; a3=fitparams[2]; a1s = fitparams[3];
    a2s = fitparams[4]; a3s=fitparams[5];  
    
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap root type
    p1.a = a1;
    p1.a_s = a1s;
    rs.setOrganRandomParameter(p1)
    p4 = rs.getRootRandomParameter(4)  # basal root type
    p4.a = a1;
    p4.a_s = a1s;
    rs.setOrganRandomParameter(p4)
    p2 = rs.getRootRandomParameter(2)
    p2.a = a2;
    p2.a_s = a2s;
    rs.setOrganRandomParameter(p2)
    p3 = rs.getRootRandomParameter(3)
    p3.a = a3;
    p3.a_s = a3s;
    rs.setOrganRandomParameter(p3)

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

            
            radius = np.asarray(core_analyser.getParameter("radius"))
            length = np.asarray(core_analyser.getParameter("length"))

            prod = 0
            for pp in range(len(radius)):
                if radius[pp]>0.09:
                    radius[pp] = 0.09
                if radius[pp] == 0:
                    radius[pp] = 0.003
                prod = prod+radius[pp]*length[pp]

            mean_diam = prod/np.sum(length)*2
            rm_[j,i] = mean_diam #matrix that contains the root length density


    rm = np.nanmean(rm_,axis = 1)
    std = np.nanstd(rm_,axis = 1)
    rmstd = np.concatenate((rm, std))
    np.savetxt("diams/p2_virt_rm_std_"+mat+"_"+den+".txt",rmstd,fmt='%.4f')

    #compare measured and simulated core RLDs 
    with open("diams/p2_measured_diam_std_"+mat+"_"+den+".txt") as f: #read in measured RLDs
        real_diam_ = [[float(xx) for xx in line.split()] for line in f]

    #flatten the vectors
    diam = np.reshape(rmstd, (corenum*2,1))
    real_diam = np.reshape(np.array(real_diam_), (corenum*2,1))
    np.set_printoptions(suppress=True)
    print(diam)
    print(real_diam)

    #momentary
    #ana = pb.SegmentAnalyser(rs)
    #ana.write("results_diam/RS_phase2.vtp") # Write all RS of last sims into one single file


    NRMSE = math.sqrt(sum(((np.subtract(diam,real_diam)))**2)/len(diam))
    err = NRMSE
    print('error= ', err)
    return err

######################################################################################
# Set geometry of soil cylinder,  top radius 5 cm, bot radius 5 cm, height 18 cm
bigcyl = pb.SDF_PlantContainer(5, 5, 18, False)

#Parameters to be fitted and initial values
a1 = 0.1322/2 #cm/d
a2 = 0.0630/2
a3 = 0.0019/2 #cm
a1s = 0.0021
a2s = 0.0038
a3s = 0.0076

bnds = ((a2, 0.066), (0.003, a2), (0.003, a2), (0, 1), (0, 1), (0, 1))

p0=[a1, a2, a3, a1s, a2s, a3s] #initial guess
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

res = scipy.optimize.minimize(err, p0, method='L-BFGS-B',bounds = bnds) #optimize params, the tolerance is set smaller than the default to speed up the optim process

print("initial guess= ", p0)
print("optimized params = ", res.x)

#save the final root system parameters
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
p1 = rs.getRootRandomParameter(1)
p1.a = res.x[0]
p1.a_s = res.x[3]
rs.setOrganRandomParameter(p1)

p4 = rs.getRootRandomParameter(4)
p4.a = res.x[0]
p4.a_s = res.x[3]
rs.setOrganRandomParameter(p4)

p2 = rs.getRootRandomParameter(2)
p2.a = res.x[1]
p2.a_s = res.x[4]
rs.setOrganRandomParameter(p2)

p3 = rs.getRootRandomParameter(3)
p3.a = res.x[2]
p3.a_s = res.x[5]
rs.setOrganRandomParameter(p3)

rs.writeParameters("rootsystem_params/p2_Optimization_diam.xml")





