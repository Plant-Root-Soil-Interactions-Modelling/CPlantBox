""" updated soil core example """

import sys
sys.path.append("../../..")
import plantbox as pb
import vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

path = "data/"
name = "Zeamays_synMRI_modified"


# User tropism 1: print input arguments to command line
class My_Info_Tropism(pb.Tropism):

    def tropismObjective(self, pos, old, a, b, dx, root):
        print("Postion \t", pos)
        print("Heading \t", old.column(0))
        print("Test for angle alpha = \t", a)
        print("Test for angle beta = \t", b)
        print("Eesolution of next segment \t", dx)
        print("Root id", root.getId())
        print()
        return 0.

class My_Age_Tropism(pb.Tropism):

    def __init__(self, rs, n, sigma, age):
        super(My_Age_Tropism, self).__init__(rs)
        self.plagio = pb.Plagiotropism(rs, 0., 0.)
        self.gravi = pb.Gravitropism(rs, 0., 0.)
        self.setTropismParameter(n, sigma)
        self.age = age

    def tropismObjective(self, pos, old, a, b, dx, root):
        age = root.getAge()
        if age < self.age:
            d = self.plagio.tropismObjective(pos, old, a, b, dx, root)
            return d
        else:
            return self.gravi.tropismObjective(pos, old, a, b, dx, root)


def simulate_rs(times :list, rs):
    """ 
    Simulates all root systems for 

    @param times     simulation times (days)
    @param allRS     list of root systems to simulate
    """
    t = times.copy()
    t.insert(0, 0.)
    dt_ = np.diff(np.array(t))
    for dt in dt_:
        rs.simulate(dt)

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
        p.lns = 0
        p.rs = 0
        p.lmaxs = 0
        p.thetas = 0

def err(fitparams):
    r = fitparams[0]; ln2=fitparams[1];  ln3 = fitparams[2];
    delayRC = fitparams[3]; delayB = fitparams[4];
    tropismN=fitparams[5]; tropismS=fitparams[6]; agegravi = fitparams[7]
    
    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap and basal root type
    #p1.lmax = lmax
    p1.r = r;
    p1.ln = ln2;
    rs.setOrganRandomParameter(p1)
    p2 = rs.getRootRandomParameter(2)
    p2.ln = ln3
    rs.setOrganRandomParameter(p2)
    srp = rs.getRootSystemParameter()
    srp.delayRC = delayRC
    srp.delayB = delayB
    rs.setRootSystemParameter(srp)

    #set the standard deviation of all parameters zero to facilitate optimization
    set_all_sd(rs)

    #define the cores 
    A = np.zeros((10,4))
    A = [[10,    0,     0,    -20],
        [10,     0,   -20,    -40],
        [10,     0,   -40,    -60]]
    A_ = np.array(A)

    x = A_[:,0]
    y = A_[:,1]
    h = np.subtract(A_[:,2],A_[:,3])
    up = A_[:,2]
    dow = A_[:,3]

    cores = soil_cores(x, y, r1, h, up)
    
    #make simulations
    rs.setSeed(0)
    rs.initialize()
    
    #define age dependent tropism
    mytropism1 = My_Info_Tropism(rs)
    mytropism1.setTropismParameter(2., 0.2)
    mytropism2 = My_Age_Tropism(rs, tropismN, tropismS, agegravi)
    rs.setTropism(mytropism2,1)
    
    simulate_rs(times, rs)

    # Export results as single vtp files (as polylines)
    if exportVTP: #only the last root system is exported as vtp
        ana = pb.SegmentAnalyser(rs)
        ana.mapPeriodic(interrow, row)
        ana.write("results/RS_tot_Leitner_Streuber.vtp") # Write all RS of last sims into one single file

    for i, t in enumerate(times):

        for j in range(0, len(cores)): 
            core_analyser = get_result(rs, t)
            core_analyser.mapPeriodic(interrow, row)
            core_analyser.crop(cores[j]);
            core_analyser.pack()
            tl1 = core_analyser.distribution("length", up[j], dow[j], 1, True)  # vertical length distribution
            tl1 = np.array(tl1) / ( r1 * r1 * math.pi * h[j])  #RLD within the ind cores in cm/cm³
            rm[j,i] = tl1 #matrix that contains the root length density

            if exportVTP: #export only the last root systems 
                vtp_name = "results/"+ name + "_core_cropped" + str(j) + "_time"+str(i)+".vtp";  # export cropped core for vizualisaten
                core_analyser.write(vtp_name);


    #compare measured and simulated core RLDs 
    with open("RLD_measured/measured_RLD_"+soil+"_"+genotype+".txt") as f: #read in measured RLDs 
        real_RLD_ = [[float(xx) for xx in line.split()] for line in f]

    #flatten the vectors
    RLD = np.reshape(rm, (corenum*len(times),1))
    real_RLD = np.reshape(np.array(real_RLD_), (corenum*len(times),1))

    #np.savetxt("virt_data_"+soil+"_"+genotype+".txt",rm,fmt='%.2f')

    NRMSE = math.sqrt(sum(((np.subtract(RLD,real_RLD)))**2)/len(RLD))
    err = NRMSE
    print('error= ', err)
    return err

######################################################################################
#Parameters to be fitted and initial values
lmax= 10000 #cm
r = 3 #cm/d
ln2 = 4 #cm
ln3 = 1.2 #cm --> determines how many 2nd order laterals there should actually be
delayRC = 5 #d
delayB = 5 #d
tropismN=1 #[-]
tropismS=0.15 #rad
agegravi = 30 #d - after how many days does gravitropism start? 

p0=[r,ln2, ln3, delayRC,delayB,tropismN, tropismS, agegravi] #initial guess
interrow = 45 #cm
row = 20 #cm
soil = 'loam' #or sand 
genotype = 'WT'

corenum = 3 #number of soil cores taken 
times = [42, 63, 98, 154] #days after planting
r1 = 2.5  # core radius
exportVTP = False  # export single root system cores (for debugging)
rm = np.zeros((corenum,len(times))) #preallocation of result matrices

res = scipy.optimize.minimize(err, p0, method='Nelder-Mead',options={'xatol': 0.01,'fatol': 0.01}) #optimize params, the tolerance is set smaller than the default to speed up the optim process

print("initial guess= ", p0)
print("optimized params = ", res.x)

#Run final simulation and plot the final root system 
rs = pb.RootSystem()
rs.readParameters(path + name + ".xml")
p1 = rs.getRootRandomParameter(1)
#p1.lmax=res.x[0]
p1.r = res.x[0]
p1.ln = res.x[1]
rs.setOrganRandomParameter(p1)

p2 = rs.getRootRandomParameter(2)
p2.ln = res.x[2]
rs.setOrganRandomParameter(p2)

srp = rs.getRootSystemParameter()
srp.delayRC = res.x[3]
srp.delayB = res.x[4]
rs.setRootSystemParameter(srp)

rs.initialize()
mytropism2 = My_Age_Tropism(rs, res.x[5], res.x[6], res.x[7])
rs.setTropism(mytropism2,1)

rs.simulate(286, False)
ana = pb.SegmentAnalyser(rs)
ana.mapPeriodic(interrow, row)

rs.write("results/example_"+soil+"_"+genotype+".vtp")
vp.plot_roots(ana, "creationTime")
rs.writeParameters("rootsystem_params/Optimization_"+soil+"_"+genotype+".xml")

np.savetxt("virt_data_"+soil+"_"+genotype+".txt",rm,fmt='%.2f')




