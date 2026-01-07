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
name = "Zea_mays_5_Leitner_Streuber_2014_mod1"

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

    a1 = fitparams[0]; a2=fitparams[1]; a3 = fitparams[2]; a4=fitparams[3]; a5=fitparams[4];  

    rs = pb.RootSystem()
    rs.readParameters(path + name + ".xml")

    #replace the relevant parameters with the fit parameters
    p1 = rs.getRootRandomParameter(1)  # tap root type
    p1.a = a1;

    p2 = rs.getRootRandomParameter(2)
    p2.a = a2

    p3 = rs.getRootRandomParameter(3)
    p3.a = a3

    p4 = rs.getRootRandomParameter(4)  # crown root type
    p4.a = a4;

    p5 = rs.getRootRandomParameter(5) #basal roots 
    p5.a = a5


    rs.setOrganRandomParameter(p1)
    rs.setOrganRandomParameter(p2)
    rs.setOrganRandomParameter(p3)
    rs.setOrganRandomParameter(p4)
    rs.setOrganRandomParameter(p5)
    
    #make simulations
    set_all_sd(rs) #set the std zero
    #rs.setGeometry(bigcube) #set the geometry of the field cube - rather map periodic

    for k in range(0,sims):
        #rs.setSeed(0)
        rs.initializeLB(5,4)
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
                meanrad = np.sum(length*radius)/np.sum(length)
                virt_diam[dummy] =  meanrad*20
                dummy+=1

    d = {'DAS': times, 'mean_diam': virt_diam}
    df_virt = pd.DataFrame(data=d)
    df_virt.to_csv("Length_diam/virtual_diam_all.csv")


    #compare measured and simulated core RLDs 
    df_meas = pd.read_csv("Length_diam/measured_diam_all.csv")
    #print(df_meas.columns.values.tolist())
    meas_diam = df_meas["meas_diam"].loc[:].values 

    RMSE = math.sqrt(sum(((np.subtract(virt_diam,meas_diam[:-1])))**2)/len(meas_diam[:-1]))

    err = RMSE
    print('error= ', err)
    return err

######################################################################################
#bigcube = pb.SDF_PlantBox(75, 20, 45)

#Parameters to be fitted and initial values
a1 = 0.12
a2 = 0.01
a3 = 0.001
a4 = 0.12
a5 = 0.1
 
p0=[a1, a2, a3, a4, a5] #initial guess
interrow = 45 #cm
row = 20 #cm
sims = 1

times = [42, 63, 98] #, 154] #days after planting
exportVTP = False  # export single root system cores (for debugging)
virt_diam = np.zeros((len(times))) #preallocation of result matrices

starttime = timeit.default_timer()
#res = scipy.optimize.minimize(err, p0, method='L-BFGS-B', bounds=bnds, options={'maxiter':20})
#res = scipy.optimize.minimize(err, p0, method='L-BFGS-B')
res = scipy.optimize.minimize(err, p0, method='Nelder-Mead') #,options={'xatol': 0.01,'fatol': 0.01}) #optimize params, the tolerance is set smaller than the default to speed up the optim process
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






