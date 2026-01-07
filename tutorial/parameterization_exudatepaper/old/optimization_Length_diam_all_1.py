""" oprimization of xml file according to measured root length and diameter"""

import sys
sys.path.append("../../../../CPlantBox")
sys.path.append("../../../../CPlantBox/src")
import plantbox as pb
import visualisation.vtk_plot as vp
import vtk 
import math
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pandas as pd
import timeit

path = "rootsystem_params/"
name = "Zea_mays_5_Leitner_Streuber_2014_mod"
#name = "Zeamays_synMRI_modified"

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

def set_all_sd(rs):
    for p in rs.getRootRandomParameter():
        p.las = 0
        p.lbs = 0
        p.rs = 0
        p.lmaxs = 0
        p.thetas = 0

def err(fitparams):
    #r = fitparams[0]; ln1=fitparams[1]; ln1s=fitparams[2]; ln2 = fitparams[3];
    #ln2s = fitparams[4]; maxB=fitparams[5]; maxBs=fitparams[6]; delayB=fitparams[7];
    #delayRC=fitparams[8]; lmax1=fitparams[9]; lmax2=fitparams[10];

    #r = fitparams[0];  maxB=fitparams[1]; delayB=fitparams[2]; delayRC=fitparams[3];

    r = fitparams[0]; ln1=fitparams[1]; ln2 = fitparams[2]; maxB=fitparams[3]; delayB=fitparams[4];  lmax1=fitparams[5]; lmax2=fitparams[6]; lmax_crown = fitparams[7]; firstB = fitparams[8];firstSB = fitparams[9];

    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap root type
    p1.ln = ln1;
    #p1.lns = ln1s;
    p4 = rs.getRootRandomParameter(4)  # crown root type
    p4.r = r;
    p4.ln = ln1;
    #p4.lns = ln1s;
    p4.lmax = lmax_crown;

    p2 = rs.getRootRandomParameter(2)
    p2.ln = ln2
    #p2.lns = ln2s
    p2.lmax = lmax1

    p3 = rs.getRootRandomParameter(3)
    p3.lmax = lmax2
    
    srp = rs.getRootSystemParameter()
    srp.maxB = maxB
    #srp.maxBs = maxBs
    srp.delayB = delayB
    srp.firstB = firstB
    srp.firstSB = firstSB
    #srp.nC = nc
    #srp.delayRC = delayRC

    rs.setOrganRandomParameter(p1)
    rs.setOrganRandomParameter(p2)
    rs.setOrganRandomParameter(p3)
    rs.setOrganRandomParameter(p4)
    rs.setRootSystemParameter(srp)
    
    #make simulations
    set_all_sd(rs) #set the std zero
    #rs.setGeometry(bigcube) #set the geometry of the field cube - rather map periodic

    for k in range(0,sims):
        #rs.setSeed(0)
        rs.initializeDB(1,4)
        simtime = 1
        dummy = 0
        for j in range(0,times[-1]+1): 
            rs.simulate(simtime, True);
        
            if j in times and k == 0:
                rs.writeParameters("rootsystem_params/Optimization_all_test.xml")
                rs.write("rootsystem_params/Optimization_all_day"+str(j)+".vtp")
                
            if j in times: 
                length = np.array(rs.getParameter("length"))
                radius = np.array(rs.getParameter("radius"))
                rm_[dummy,k] =  np.sum(length)
                dummy+=1

    rm = np.round(np.mean(rm_,axis = 1),2)
    std_rm = np.round(np.std(rm_,axis = 1),2)
    virt_length = np.concatenate((rm, std_rm))
    d = {'DAS': times, 'mean_length': rm, 'std_length': std_rm}
    df_virt = pd.DataFrame(data=d)
    df_virt.to_csv("Length_diam/virtual_length_all.csv")


    #compare measured and simulated core RLDs 
    df_meas = pd.read_csv("Length_diam/measured_length_all.csv")
    #print(df_meas.columns.values.tolist())
    mean_l = df_meas["mean_length"].loc[:].values 
    std_l = df_meas["std_length"].loc[:].values
    real_length = np.concatenate((mean_l, std_l))

    print(rm)
    print(mean_l[:len(rm)])
    #print(std_l)
    #print(virt_length)
    #print(real_length)

    RMSE = math.sqrt(sum(((np.subtract(rm,mean_l[:len(rm)])))**2)/len(rm))
    #RMSE = math.sqrt(sum(((np.subtract(virt_length,real_length)))**2)/len(virt_length))
    err = RMSE
    print('error= ', err)
    return err

######################################################################################
#bigcube = pb.SDF_PlantBox(75, 20, 45)

#Parameters to be fitted and initial values
ln1 = 0.5 #0.3
lnc = 0.5 #0.2 #cm
ln2 = 1.2
lb1= 5
lbc = 5
la1= 10
lac = 10
firstB = 15
maxB = 40
delayB = 5
lmax2 = 6
lmax3 = 1.5

nc = 20 #number of crown roots per whirl 


p0=[ln1,lnc,ln2,lb1, lbc, la1, lac, firstB, maxB, delayB,lmax2,lmax3] #initial guess
#bnds = ((1, 5), (5, 20), (10, 40), (0, 40))
#p0=[r,maxB,delayB,delayRC] #initial guess
interrow = 45 #cm
row = 20 #cm
sims = 1

times = [42, 63, 98] #, 154] #days after planting
exportVTP = False  # export single root system cores (for debugging)
rm_ = np.zeros((len(times),sims)) #preallocation of result matrices
rm = np.zeros((len(times))) #preallocation of result matrices
std = np.zeros((len(times))) #preallocation of result matrices

starttime = timeit.default_timer()
#res = scipy.optimize.minimize(err, p0, method='L-BFGS-B', bounds=bnds, options={'maxiter':20})
#res = scipy.optimize.minimize(err, p0, method='L-BFGS-B')
res = scipy.optimize.minimize(err, p0, method='Nelder-Mead',options={'xatol': 0.01,'fatol': 0.01}) #optimize params, the tolerance is set smaller than the default to speed up the optim process
print("Optimization took:", timeit.default_timer() - starttime, " s")

print("initial guess= ", p0)
print("optimized params = ", res.x)
sys.exit()

#Run final simulation and plot the final root system 
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
rs.writeParameters("rootsystem_params/Optimization_all.xml")

sys.exit()
#rs.setGeometry(bigcube)
rs.initialize()
rs.simulate(21, False)
rs.write("results/p1_example_Maxime_"+mat+".vtp")
ana = pb.SegmentAnalyser(rs)
vp.plot_roots(ana, "creationTime")






