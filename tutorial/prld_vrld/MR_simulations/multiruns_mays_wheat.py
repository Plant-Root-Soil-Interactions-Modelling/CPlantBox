import argparse
from benchmark import run_benchmark
import os
import numpy as np
import plantbox as pb
import sys
import multiprocessing as mp
import time
import csv

def run_benchmark_(num):
    #print("run n*",num)
    seed = num+counts
    print('seed', seed) 
    fx = run_benchmark(seed, name, M, N, distr, distp, distl, simtime, w_field, l_field, hlayer, depth, planes, tube_diam, fieldbox_wo_rhizotubes, rhizotubes, rhizotubeso, sc1,sc2,sc3,sc4,sc_volume, imgw, imgl,size,x_standard, x_stand2, x_highr, x_cont, img_cont, y_, z_, tropismN, tropismS)
    return fx

def setup_scenario(plant, size, tube_diam, pdensity):
    tropismN_ = np.linspace(1,4,4) 
    tropismS_ = np.linspace(0.025, 0.2, 8) 

    #General parameters 
    if plant == 'maize': 
        name = "Zea_mays_3_Postma_2011"
        simtime = 120 #days
        if pdensity == 'base': 
            M = 8 # number of plants in one row
            N = 6 # number of rows
            distp = 12  # distance between the plants[cm] (P-P)
            distr = 75 # distance between the rows[cm] (R-R)
        elif pdensity == 'alt': 
            M = 6 # number of plants in one row
            N = 8 # number of rows
            distp = 20  # distance between the plants[cm] (P-P)
            distr = 50 # distance between the rows[cm] (R-R)
    elif plant == 'wheat': 
        name = "wheat_Morandage"
        simtime = 215 #days
        if pdensity == 'base': 
            M = 50 # number of plants in one row #2
            N = 10 # number of rows
            distp = 1.5 # distance between the plants[cm] (P-P)
            distr = 21 # distance between the rows[cm] (R-R)
        elif pdensity == 'alt': 
            M = 30 # number of plants in one row #2
            N = 16 # number of rows
            distp = 3 # distance between the plants[cm] (P-P)
            distr = 12 # distance between the rows[cm] (R-R)
            
    l_field=M*distp # length of field 
    w_field=N*distr # width of field 
    depth = 130
    layers = 13
    hlayer = depth / layers
    soilvolume = (depth / layers) * l_field * w_field #soil volume per layer
    tube_diam = float(tube_diam) #tube diameter
    tube_len = 700 #tube length 
    v_d  = 2 # viewing_depth in mm --> experiment
    ladd = 1.80 #length addition to minirhizotube diameter- should not overlap -->6.4+1.8*2= 10
    r_sc = 4.5 #radius of soil cores 
    sc_volume = r_sc**2*np.pi*hlayer
    distl = 5 #distance between perpendicular planes for coimputation of anisotropy factor 
    if size == 'small': 
        imgw = 2  # image width 
        imgl = 2  # image length 
    elif size == 'large': 
        imgw = 6
        imgl = 4 
    elif size == 'complete': 
        imgw = 2
        imgl = 20

    #Geometry parameters for rhizotron images and rhizotubes 
    y_ = [-25, -15, -5, 5, 15, 25] #y positions of rhizotubes 
    z_ = [-10, -20, -40, -60, -80, -120] #z positions of rhizotubes 

    #standard image positions 
    if plant == 'maize': 
        xstart = -(N-1)*distr/2
        x_standard = np.zeros(20) 
        dummy = 0
        for i in range(0,4): 
            for j in range(0,5): 
                x_standard[dummy] = xstart+i*93 + j*10
                dummy = dummy+1
    elif plant == 'wheat': 
        xstart = -(N-1)*distr/2
        x_standard = np.zeros(20) 
        dummy = 0
        for i in range(0,20): 
            x_standard[dummy] = xstart+i*10 #93 and 12 is spacing defined by Lena
            dummy = dummy+1
            
    #high resolution image positions 
    x_startend_ = np.zeros((2))
    arr = [0,N-1]
    for j in range(0,len(x_startend_)): 
        x_startend_[j] = distr * arr[j]-((N-1)*distr/2)
    x_stand2 = np.linspace(x_startend_[0]-distr/2, x_startend_[1]+distr/2, 20)    
    x_highr = np.linspace(x_startend_[0]-distr/2, x_startend_[1]+distr/2, 40) 

    #continuous image position 
    x_cont = x_startend_[0]-distr/2
    img_cont = (x_startend_[1]+distr/2)-(x_startend_[0]-distr/2)

    #rhizotubes, single batch 
    rhizotube_ = pb.SDF_PlantContainer(tube_diam/2, tube_diam/2, tube_len, False)  # a single rhizotube
    rhizotube = pb.SDF_RotateTranslate(rhizotube_, 90, pb.SDF_Axis.yaxis, pb.Vector3d(tube_len / 2, 0, 0))  # a single rhizotube with correct x position 
    rhizotubes_ = []
    for i in range(0, len(y_)):
        rhizotubes_.append(pb.SDF_RotateTranslate(rhizotube, pb.Vector3d(0, y_[i], z_[i])))
    rhizotubes = pb.SDF_Union(rhizotubes_) #single batch rhizotubes 
    # rs = pb.RootSystem()
    # rs.setGeometry(rhizotubes)
    # rs.write("test_results/rhizotubes.py")

    #rhizotubes increased by viewing depth, single batch 
    rhizotubeo_ = pb.SDF_PlantContainer(tube_diam/2+v_d/10, tube_diam/2+v_d/10, tube_len, False)
    rhizotubeo = pb.SDF_RotateTranslate(rhizotubeo_, 90, pb.SDF_Axis.yaxis, pb.Vector3d(tube_len / 2, 0, 0))
    rhizotubes_o = []
    for i in range(0, len(y_)):
        rhizotubes_o.append(pb.SDF_RotateTranslate(rhizotubeo, pb.Vector3d(0, y_[i], z_[i])))
    rhizotubeso = pb.SDF_Union(rhizotubes_o)

    #field box with rhizotubes as obstacles 
    fieldbox = pb.SDF_PlantBox(w_field*1.5, l_field*1.5, depth+70)  
    fieldbox_wo_rhizotubes = pb.SDF_Difference(fieldbox, rhizotubes)


    #soil cores  - two interrow, two interplant 
    sc_ = pb.SDF_PlantContainer(r_sc, r_sc, hlayer, False)  
    sc = []
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * np.floor(N/3)-((N-1)*distr/2)+distr/2, distp * np.floor(M/3)-(distp*M/2), 0)))  # interrow soil core 1
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * np.floor(N/3)-((N-1)*distr/2), distp * np.floor(M/3)-(distp*M/2)+distp/2, 0)))  # interplant soil core 1
    sc1 = pb.SDF_Union(sc)
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * 2*np.floor(N/3)-((N-1)*distr/2)+distr/2, distp * 2*np.floor(M/3)-(distp*M/2), 0)))  # interrow soil core 2
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * 2*np.floor(N/3)-((N-1)*distr/2), distp * 2*np.floor(M/3)-(distp*M/2)+distp/2, 0)))  # interplant soil core 2
    sc2 = pb.SDF_Union(sc)

    sc = []
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * np.floor(N/3)-((N-1)*distr/2)+distr/2, distp * np.floor(M/3)-(distp*M/2), 0)))  # interrow soil core 1
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * 2*np.floor(N/3)-((N-1)*distr/2)+distr/2, distp * 2*np.floor(M/3)-(distp*M/2), 0)))  # interrow soil core 2
    sc3 = pb.SDF_Union(sc)

    sc = []
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * np.floor(N/3)-((N-1)*distr/2), distp * np.floor(M/3)-(distp*M/2)+distp/2, 0)))  # interplant soil core 1
    sc.append(pb.SDF_RotateTranslate(sc_, pb.Vector3d(distr * 2*np.floor(N/3)-((N-1)*distr/2), distp * 2*np.floor(M/3)-(distp*M/2)+distp/2, 0)))  # interplant soil core 2
    sc4 = pb.SDF_Union(sc)
    # Computation of planes for anisotropy factor 
    #x
    x0x = pb.Vector3d(-1., 0., 0.)
    nxx = pb.Vector3d(-1., l_field, 0.)
    nyx = pb.Vector3d(-1, 0., depth)
    #y
    x0y = pb.Vector3d(0., -1, 0.)
    nxy = pb.Vector3d(w_field, -1, 0.)
    nyy = pb.Vector3d(0, -1, depth)
    #z
    x0z = pb.Vector3d(0., 0., -1.)
    nxz = pb.Vector3d(w_field, 0., -1.)
    nyz = pb.Vector3d(0., l_field, -1.)
    planes = [[x0x, nxx, nyx],[x0y, nxy, nyy],[x0z, nxz, nyz]]
    
    return name, M, N, distr, distp, distl, simtime, w_field, l_field, hlayer, depth, planes, tube_diam, pdensity, fieldbox_wo_rhizotubes, rhizotubes, rhizotubeso, sc1,sc2,sc3,sc4,sc_volume, imgw, imgl,size,x_standard, x_stand2, x_highr, x_cont, img_cont, y_, z_, tropismN_, tropismS_

def write_results(name, simtime, repetitions, imgw, imgl, tube_diam, pdensity, finalr): 
    #write out results
    header = ['depth', 'tropismN', 'tropismS', 'rRLD', 'An','CV',
        'prld_stand', 'prld_stand2', 'prld_highr','prld_cont',
        'vrld_sl', 'vrld_2sc', 'vrld_4sc','vrld_sc_ir', 'vrld_sc_ip',
        'prvd_stand','prvd_stand2','prvd_highr','prvd_cont',
        'vrvd_sl', 'vrvd_2sc', 'vrvd_4sc',
        'prcd_stand','prcd_stand2','prcd_highr','prcd_cont'] 
    with open('results/'+name+'_day'+str(simtime)+'_reps'+str(repetitions)+'_imgsize'+str(imgw)+'x'+str(imgl)+'tubediam'+str(tube_diam)+'_pdensity_'+pdensity+'.csv', 'w') as f:
        writer=csv.writer(f, delimiter=',', lineterminator='\n')
        writer.writerow(header)
        for row in finalr:
            writer.writerow(row)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Simulation options')
    parser.add_argument('plant', type = str, help = 'maize or wheat')
    parser.add_argument('size', type = str, help = 'small or large or complete')
    parser.add_argument('tubediam', type = str, help = '5.7, 6.4 or 7')
    parser.add_argument('pdensity', type = str, help = 'base or alt')

    args = parser.parse_args()
    name_id = args.plant + "_" + args.size + "_img_tubediam" + args.tubediam
    print()
    print(name_id, "\n")

    name, M, N, distr, distp, distl, simtime, w_field, l_field, hlayer, depth, planes, tube_diam, pdensity, fieldbox_wo_rhizotubes, rhizotubes, rhizotubeso, sc1,sc2,sc3,sc4,sc_volume, imgw, imgl,size,x_standard, x_stand2, x_highr, x_cont, img_cont, y_, z_, tropismN_, tropismS_ = setup_scenario(args.plant, args.size, args.tubediam, args.pdensity)

    repetitions = 100
    simnum = 0
    finalr = []
    finalr = np.asarray(finalr) 
    for j in range(0, len(tropismN_)): 
        for k in range(0, len(tropismS_)): 


            tropismN = tropismN_[j]
            tropismS = tropismS_[k]

            reps = repetitions
            counts = 0
            while reps>0: 
                start_time = time.time()
                processes = 25
                pool = mp.Pool(processes=processes)
                print('Simulation', simnum, ' with : ', name, 'tropismN = ', tropismN_[j], 'tropismS = ', tropismS_[k], 'repetitions left = ', reps) 
                result=np.array((pool.map(run_benchmark_, np.arange(processes))))
                result = result.reshape((result.shape[0]*result.shape[1]), result.shape[2])
                print("Simulations took ", (time.time() - start_time), 's') 
                if finalr.size == 0: 
                    finalr = result
                else: 
                    finalr = np.vstack((finalr,result))
                print(np.shape(finalr)) 
                counts = counts+processes
                reps = reps-processes
                
            simnum = simnum+1
    write_results(name, simtime, repetitions, imgw, imgl, tube_diam, pdensity, finalr)        


    # python3 multiruns_mays_wheat.py maize small 6.4 base