"""" CF (maize) Selhausen Rhizotron facility"""
import sys
sys.path.append("../..")
sys.path.append("../../src")

import plantbox as pb
import math
import numpy as np
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt

name = "Zea_mays_1_Leitner_2010"

simtime = 180

M = 8 # number of plants in rows
N = 6 # number of rows

distp = 12.5  # distance between the root systems along row[cm] (P-P)
distr = 75 # distance between the rows[cm] (R-R)
interrow=M*distp # intra-row spacing
row=N*distr # row spacing

depth = 130 #cm
layers = 13
soilvolume = (depth / layers) * interrow * row #volume of the soil layer

box = pb.SDF_PlantBox(5000, 5000, 5000)  # box
rhizotube = pb.SDF_PlantContainer(6.4/2, 6.4/2, 700, False)  # a single rhizotube
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(700 / 2, 0, 0))
rhizotubes_ = []
y_ = (-30, -18, -6, 6, 18, 30)
z_ = (-10, -20, -40, -60, -80, -120)
tube = []
for i in range(0, len(y_)):
    v = pb.Vector3d(0, y_[i], z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_) #six minirhizotubes at six different depths

rhizotube_mirror1=pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, distp*M, 0))
rhizotube_mirror2=pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, -distp*M, 0))
rhizotuben=pb.SDF_Union([rhizotube_mirror1,rhizotubes,rhizotube_mirror2]) # Mirrored rhizotubes along with existing

v_d  = 1.24 # viewing_depth in mm
rhizotubeo = pb.SDF_PlantContainer((6.4/2)+((v_d)/10), (6.4/2)+((v_d)/10), 700, False) # outer tube including v_d  in addition to existing rhizotube size 
rhizoXo = pb.SDF_RotateTranslate(rhizotubeo, 90, pb.SDF_Axis.yaxis, pb.Vector3d(700 / 2, 0, 0))

rhizotubes_o = []
tubeo = []
for i in range(0, len(y_)):
    vo = pb.Vector3d(0, y_[i], z_[i])
    tubeo.append(pb.SDF_RotateTranslate(rhizoXo, vo))
    rhizotubes_o.append(tubeo[i])

rhizotubeso = pb.SDF_Union(rhizotubes_o)

rhizoTube = pb.SDF_Difference(box, rhizotuben) # we set this geometry during root growth which means roots will grow in a big box geometry excluding the space where rhizotubes and mirrored rhizotubes exist

opening = pb.SDF_PlantBox(6, 4,10)  # image size inhouse facility box

#r and l are two different sides of tubes at 80 degrees. r0 and l0 are at depth 10 cm and r5 and l5 are at depth of 120 cm
r0 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[0], z_[0]))
l0 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[0], z_[0]))

r1 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[1], z_[1]))
l1 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[1], z_[1]))

r2 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[2], z_[2]))
l2 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[2], z_[2]))

r3 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[3], z_[3]))
l3 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[3], z_[3]))

r4 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[4], z_[4]))
l4 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[4], z_[4]))

r5 = pb.SDF_RotateTranslate(opening, 180-80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[5], z_[5]))
l5 = pb.SDF_RotateTranslate(opening, 180+80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[5], z_[5]))

# # array of 20 image positions on one side of rhizotube
# x_i = np.array([1039.5,1120.5,1255.5,1377,1498.5,
                           # 2038.5,2119.5,2254.5,2376,2497.5,
                           # 3037.5,3118.5,3253.5,3375,3496.5,
                           # 4036.5,4117.5,4252.5,4374,4495.5])
# pos = x_i/10
# lim_i = (pos[-1]-pos[0])
# lim_i_m = lim_i/2 #mid
# x_ = pos+lim_i_m-pos[0]-(700/2)

# array of additional image positions such that all images are adjacent to each other
x_i = np.array([1039.5,1120.5,1255.5,1377,1498.5,
                           2038.5,2119.5,2254.5,2376,2497.5,
                           3037.5,3118.5,3253.5,3375,3496.5,
                           4036.5,4117.5,4252.5,4374,4495.5])
pos = x_i/10
lim_i = (pos[-1]-pos[0])
lim_i_m = lim_i/2 #mid
x__ = pos+lim_i_m-pos[0]-(700/2)
x_ = np.arange(x__[0],x__[-1]+6,6)
                           

def run_benchmark():  # make a function to be called certain times
    # Initializes N*M root systems
    allRS = []
    for i in range(0, M):
        for j in range(0, N):
            rs = pb.RootSystem()
            rs.readParameters(name + ".xml")
            rs.getRootSystemParameter().seedPos = pb.Vector3d(distr * j-((N-1)*distr/2), distp * i-(distp*M/2), -3.) # wheat rows perpendicular to tubes
            rs.setGeometry(rhizoTube) #rhizotubes act as obstacles
            rs.initialize(False)
            # Simulate
            rs.simulate(simtime, False)
            allRS.append(rs)
            if i+ j == 0:
                allAna1 = pb.SegmentAnalyser(rs) 
            else:
                allAna1.addSegments(rs)  # collect all in a segAna object

    # allAna1.write("all_plants.vtp")
    allAna1.mapPeriodic(row, interrow)
    # allAna1.write("mp.vtp")

    rl_ = []
    rl_.append(allAna1.distribution("length", 0., -depth, layers, True))
    rl_ = np.array(rl_[-1]) / soilvolume  # convert to density

    allAna = pb.SegmentAnalyser(allAna1)
    allAna.crop(pb.SDF_Difference(rhizotubeso, rhizotubes))
    # allAna.write("hc.vtp")  #all roots potentially visible all around the tube

    sa= 6 * 4 # surface_area for surface density

    ls=[r0,l0, r1,l1, r2,l2, r3,l3, r4,l4, r5,l5] # list containing cropping geometries same as image size on both sides of tube and at different depths

    ls_not=['r0_','l0_', 'r1_','l1_', 'r2_','l2_', 'r3_','l3_', 'r4_','l4_', 'r5_','l5_'] #list of notations
    es=np.zeros((x_.size,len(ls)))
    

    for j in range(len(ls)):
        for i in range(es.shape[0]):
            v = pb.Vector3d(x_[i], 0, 0)
            g=pb.SDF_RotateTranslate(ls[j], v)
            ana = pb.SegmentAnalyser(allAna)
            ana.crop(g)
            ana.pack()
            rsd_each = (ana.getSummed("length"))/sa #root surface density
            es[i, j] = rsd_each

    t1 = es[:,:2] # pRLD from all images on two sides of tube 1 which is at 10 cm depth
    t2 = es[:,2:4]
    t3 = es[:,4:6]
    t4 = es[:,6:8]
    t5 = es[:,8:10]
    t6 = es[:,10:]
    pRLD = np.array([np.mean(t1), np.mean(t2), np.mean(t3), np.mean(t4), np.mean(t5), np.mean(t6)]) #mean pRLD at six different depths
    res_arr=np.zeros((pRLD.size, 6))
    
    #np.mean([rl_[1-1],rl_[2-1]]) gives vRLD @ 10 cm .... 
    vRLD = np.array([np.mean([rl_[1-1],rl_[2-1]]),
    np.mean([rl_[2-1],rl_[3-1]]),
    np.mean([rl_[4-1],rl_[5-1]]),
    np.mean([rl_[6-1],rl_[7-1]]),
    np.mean([rl_[8-1],rl_[9-1]]),
    np.mean([rl_[12-1],rl_[13-1]]),
     ]) #mean vRLD at six different depths
    
    res_arr[:,0] = z_ #depths
    res_arr[:,1] = pRLD
    res_arr[:,2] = vRLD
    res_arr[:,3] = vRLD / pRLD # CF (relevant).
    res_arr[:,4] = pRLD*sa # mean root length 
    res_arr[:,5] = (pRLD*sa) / vRLD # CF(RL_image/vRLD)
    return res_arr


if __name__ == '__main__': #for testing
    import time
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    start_time = time.time()
    outcomes = run_benchmark()
    df = pd.DataFrame(outcomes[:, 1:], columns=["pRLD", "vRLD", "CF(vRLD/pRLD)", "RL", "CF(RL/vRLD)"], index=outcomes[:, 0])
    print(df)   
    print("--- %s seconds, end benchmark ---" % (time.time() - start_time))