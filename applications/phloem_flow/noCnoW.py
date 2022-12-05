""" water movement within the root (static soil) """


import sys;
import os
CPBdir = "../.."

main_dir=os.environ['PWD']#dir of the file
directoryN = "/"+os.path.basename(__file__)[:-3]+"/"#"/a_threshold/"
results_dir = main_dir +"/results"+directoryN

sys.path.append("../.."); sys.path.append("../../src/python_modules")
#import matplotlib
#matplotlib.use('AGG') #otherwise "import matplotlib.pyplot" hangs

if not os.path.exists(results_dir):
    os.makedirs(results_dir)
else:
   test = os.listdir(results_dir)
   for item in test:
        try:
            os.remove(results_dir+item)
        except:
            pass
        
from phloem_flux import PhloemFluxPython  # Python hybrid solver
#from Leuning_speedup import Leuning #about 0.7 for both
#from photosynthesis_cython import Leuning
import plantbox as pb
#from plantbox import Photosynthesis as ph
#import vtk_plot as vp
import math
import os
import numpy as np
import vtk_plot as vp

#import matplotlib.pyplot as plt
from datetime import datetime, timedelta



isCluster = (os.environ['HOME'] == '/home/m.giraud')
if False:
    def print(*args, **kwargs):
        """ custom print() function.
            for cluster: can get output even if program stop
            unexpectedly (i.e., without creating the outputfile)
        """
        # Adding new arguments to the print function signature
        # is probably a bad idea.
        # Instead consider testing if custom argument keywords
        # are present in kwargs
        if 'sep' in kwargs:
            sep = kwargs['sep']
        else:
            sep = ' '
        home_dir = os.getcwd()
        dir_name =  "/results"+directoryN
        dir_name2 = home_dir + dir_name
        name2 = dir_name2 + 'prints.txt'
        with open(name2, 'a') as log:
            for arg in args: log.write(str(arg) + sep)
            log.write('\n')

#time lapse pea plants:
#https://www.youtube.com/watch?v=a2Mxpg5lNFM
#       qr, qs, alpha, n, ks (in [cm/d])
#vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
def theta2H(vg,theta):#(-) to cm
    thetar =vg[0]# 0.059
    thetas = vg[1]#0.445
    alpha = vg[2]#0.00644
    n = vg[3]#1.503
    nrev = 1/(1-1/n)
    H =-(((( (thetas - thetar)/(theta - thetar))**nrev) - 1)**(1/n))/alpha
    return(H)#cm

def sinusoidal(t):
    return (np.sin(np.pi*t*2)+1)/2

def qair2rh(qair, es_, press):
    e =qair * press / (0.378 * qair + 0.622)
    rh = e / es_
    rh=max(min(rh, 1.),0.)
    return rh


def weather(simDuration):
    vgSoil = [0.059, 0.45, 0.00644, 1.503, 1]
    
    Qmin = 0; Qmax = 960e-6 #458*2.1
    Tmin = 15.8; Tmax = 22
    specificHumidity = 0.0097
    Pair = 1010.00 #hPa
    thetaInit = 30/100

    coefhours = sinusoidal(simDuration)
    TairC_ = Tmin + (Tmax - Tmin) * coefhours
    Q_ = Qmin + (Qmax - Qmin) * coefhours
    cs = 850e-6 #co2 paartial pressure at leaf surface (mol mol-1)
    #RH = 0.5 # relative humidity
    es =  6.112 * np.exp((17.67 * TairC_)/(TairC_ + 243.5))
    RH = qair2rh(specificHumidity, es, Pair)
    
    pmean = theta2H(vgSoil, thetaInit)
    
    weatherVar = {'TairC' : TairC_,
                    'Qlight': Q_,
                    'cs':cs, 'RH':RH, 'p_mean':pmean, 'vg':vgSoil}
    print("Env variables at", round(simDuration//1),"d",round((simDuration%1)*24),"hrs :\n", weatherVar)
    return weatherVar

def div0(a, b, c):        
    return np.divide(a, b, out=np.full(len(a), c), where=b!=0)
    
def div0f(a, b, c):    
    if b != 0:
        return a/b
    else:
        return a/c
        
def write_file_array(name, data):
    name2 = 'results'+ directoryN+ name+ "_"+ strQ + "_"+strTh+"_"+strDecap+ '.txt'
    with open(name2, 'a') as log:
        log.write(','.join([num for num in map(str, data)])  +'\n')

def write_file_float(name, data):
    #name2 = 'results' + directoryN+  name+ "_"+ strQ + "_"+strTh+"_"+strDecap+ '.txt'
    name2 = 'results' + directoryN+  name+  '.txt'
    with open(name2, 'a') as log:
        log.write(repr( data)  +'\n')

""" Parameters """

weatherInit = weather(0)
simInit = 1
simDuration = simInit # [day] init simtime
simMax =50
depth = 60
dt = 1/24 #1h
verbose_simulate = False

# plant system 
pl = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
pl2 = pb.MappedPlant(seednum = 2) #pb.MappedRootSystem() #pb.MappedPlant()
path = CPBdir+"/modelparameter/plant/"
name = "UQ_1Leaf"

print(path + name + ".xml")

pl.readParameters(path+ name + ".xml")


#raise Exception
sdf = pb.SDF_PlantBox(np.Inf, np.Inf, depth )

pl.setGeometry(sdf) # creates soil space to stop roots from growing out of the soil


for p in pl.getOrganRandomParameter(pb.stem):
    if (p.subType > 0):
        print(p.subType, "radius", p.a, "lmax", p.lmax, p.ln, p.lb,  p.nob())

pl.initialize(verbose = False)#, stochastic = False)
pl.simulate(1, verbose_simulate)#, "outputpm15.txt")#simDuration

ana = pb.SegmentAnalyser(pl.mappedSegments())
#OrgIds =np.array(list(map(lambda x: np.array(x), pl.nodes)))
np.array([org.getId() for org in pl.getOrgans()])
#ana.addData("OrgId", Q_Exud_i_p)
ana.write("results"+directoryN+name + "_"+ str(0) +".vtp")

orgs = pl.getOrgans()
#print(p.segLength())

orgs = pl.getOrgans()
pl.getNumberOfOrgans()
ö = 0

#vp.plot_plant(pl,p_name = "organType")
#simDuration=50
while simDuration <= simMax:#7
    #a1 = np.array(list(map(lambda x: np.array(x), pl.segments)), dtype = np.int64)
    #a2 = np.array(list(map(lambda x: np.array(x), pl.plant().getSegments())), dtype = np.int64)
    print('simDuration:',simDuration )
    nodeId = np.array([i for i in range(len(pl.nodes))])
    #print(nodeId)
    ana = pb.SegmentAnalyser(pl.mappedSegments())
    ana.addData("nodeId",nodeId)
    #a2 = [[seg.x,seg.y] for seg in ana.segments]
    ana.write("results"+directoryN+"msx_"+ str(ö) +".vtp", ["organType","subType", "nodeId"])
    if(sum(np.isnan(nodeId))>0):
        raise Exception
    #raise Exception
    # #a3 = [[seg.x,seg.y] for seg in ana.segments]
    # # print("do plant")
    # ana = pb.SegmentAnalyser(pl.plant())
    # # print("write")
    # # a5 = [[seg.x,seg.y] for seg in ana.segments]
    # ana.write("results"+directoryN+"plx_"+ str(ö) +".vtp", ["organType","subType", "nodeId"])
    # a6 = [[seg.x,seg.y] for seg in ana.segments]
    
    ö +=1

    pl.simulate(dt, verbose_simulate)#, "outputpm15.txt") #time span in days /(60*60*24)
    #pl2.simulate(dt,  verbose_simulate)#, "outputpm15.txt")
    simDuration += dt
    #print(a1)
    #print(np.sort(a1, axis = 1))
    # print([[seg.x,seg.y] for seg in pl.segments]==a2)
    # print(a2==a3)
    # print(a1)
    

# for p in pl.getOrgans(pb.stem):
#     print(len(pl.getOrgans(pb.stem)), p.getId(), p.getLength())
# print("leaf")
print("simDuration", simDuration, "d", len(pl.getOrgans(pb.leaf)))
#vp.plot_plant(pl,p_name = "organType")
print("results"+directoryN)