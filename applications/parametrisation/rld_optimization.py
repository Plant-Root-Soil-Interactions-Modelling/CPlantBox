"""Optimization of CPlantBox Mdoel tap and basal root parameters"""
import sys; sys.path.append("../..")
import numpy as np
import matplotlib.pyplot as plt
import vtk_plot as vp
import plantbox as pb
import pickle as pkl
import math
from scipy.optimize import differential_evolution

def err (fitparams):    #parameters to be optimized
    lmaxp=fitparams[0]
    tropismNp=fitparams[1] 
    tropismSp=fitparams[2]
    rp=fitparams[3]

    simtime = 120
    M = 4  # number of plants in rows
    N = 2 # number of rows
    distp = 3  # distance between the root systems along row[cm]
    distr =12  # distance between the rows[cm]
    interrow=N*distr # inter-row spacing
    row=M*distp # row spacing

    r, depth, layers = 5, 100., 11 # Soil core analysis
    layerVolume = depth / layers * r * r * np.pi

    z_ = np.linspace(0, -1 * depth, layers)
    times = [120, 60, 30, 10]
    soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
    soilcolumn1 = pb.SDF_RotateTranslate(soilcolumn, 0, 0, pb.Vector3d(-6, 0, 0))
    soilcolumn2 = pb.SDF_RotateTranslate(soilcolumn, 0, 0, pb.Vector3d(6, 0, 0))

    soilcolumns=[soilcolumn1,soilcolumn, soilcolumn2]

    with open('rld.pkl','rb') as f:
        measured_RLD = pkl.load(f)
    real_RLD=np.reshape(measured_RLD,(measured_RLD.shape[0]*measured_RLD.shape[1]*measured_RLD.shape[2],1))
    rld=np.zeros([measured_RLD.shape[0],measured_RLD.shape[1],measured_RLD.shape[2]])
    # rld=np.zeros([len(soilcolumns),len(times),layers])

    path = "../../modelparameter/rootsystem/"
    name = "wheat"
    

    # fig, axes = plt.subplots(nrows = 1, ncols = len(soilcolumns), figsize = (16, 8))

    # Make a root length distribution along the soil cores  

    for k in range(len(soilcolumns)):
        # Initializes N*M root systems
        allRS = []
        for i in range(0, N):
            for j in range(0, M):
                rs = pb.RootSystem()
                rs.readParameters(path + name + ".xml")
                p1 = rs.getRootRandomParameter(1)  # tap and basal root type
                p1.lmax = lmaxp
                p1.tropismN=tropismNp
                p1.tropismS=tropismSp
                p1.r=rp
                for p in rs.getRootRandomParameter():
                    p.lns = 0
                    p.rs = 0
                    p.lmaxs = 0
                    p.thetas = 0
                    p.las = 0
                    p.lbs=0
                rs.setSeed(1)    
                rs.getRootSystemParameter().seedPos = pb.Vector3d(distr * i, distp * j, -3.)  # cm
                rs.initialize()
                allRS.append(rs)

        # Simulate
        for rs in allRS:
            rs.simulate(simtime)

        # Export results as single vtp files (as polylines)
        ana = pb.SegmentAnalyser()  # see example 3b
        for z, rs in enumerate(allRS):
            # vtpname = "results/rlds/plant_" + str(z) + ".vtp"
            # rs.write(vtpname)
            ana.addSegments(rs)  # collect all

        # Write all into single file (segments)
        # ana.write("results/rlds/all_plants.vtp")

        ana.mapPeriodic(interrow, row)
        # ana.write("results/rlds/plant_periodic.vtp")
        rl_ = []
        ana.crop(soilcolumns[k])
        ana.pack()
        # ana.write("results/rlds/core"+str(k+1)+".vtp")

        # axes[k].set_title('Soil core'+' ' +str(k+1))
        for j in range(len(times)):
            ana.filter("creationTime", 0, times[j])
            rl_.append(ana.distribution("length", 0., -depth, layers, True))
            # axes[k].plot(np.array(rl_[-1]) / layerVolume, z_)
        # axes[k].legend(["120 days", "60 days", "30 days", "15 days"])
        rld[k]=np.array(rl_)/layerVolume 

    # for a in axes:
    #     a.set_xlabel('RLD $(cm/cm^3)$')
    #     a.set_ylabel('Depth $(cm)$')
    #     a.set_xlim(0,np.max(rld))

    # fig.subplots_adjust()
    # plt.savefig("results/rlds/rld_plot.png")
    # plt.show()

    RLD=np.reshape(rld,(rld.shape[0]*rld.shape[1]*rld.shape[2],1))
    # print(rld)
    # with open('results/rlds/rld.pkl','wb') as f:
    #     pkl.dump(rld, f)
    # sys.exit()
    
    err = math.sqrt(sum(((np.subtract(RLD,real_RLD)))**2)/len(RLD)) #NRMSE
    return err


bnds = ([60,150],[0, 5],[0, 5],[0,5])     #[lmax,tropismN,tropismS,r]

res = differential_evolution(err, bnds)
print(res.x, res.fun)
print(res)

np.savetxt("opt_params.txt", np.around(res.x, decimals=9),fmt='%.9f',delimiter=',')
 
# lmaxp= 80 #cm
# tropismNp=1 #cm 
# tropismSp=0.5 #cm/d
# rp=2.2
# p0=[lmaxp,tropismNp,tropismSp,rp]
# print(err(p0))
