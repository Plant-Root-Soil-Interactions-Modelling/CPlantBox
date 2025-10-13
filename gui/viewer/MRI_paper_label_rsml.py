import sys; sys.path.append("../../src/python_modules/"); sys.path.append("../../")

import os
import pandas as pd
import numpy as np
import argparse
import plantbox as pb

import viewer_conductivities
from viewer_data import ViewerDataModel
import xylem_flux
from vtk_tools import * 
from collections import Counter 

"""
creates a csv file per rsml file containing krs values, and suf values per 1 mm layers

expected three arguments: file path, scenario index (1-4), 
optionally, --shoot, only use for measured monocots; --split, for multiple dicots, --z_shift to relocate seed to -3cm

python3 label_rsml.py /home/daniel/workspace/DUMUX/CPlantBox/gui/estimate/img/monocot/ 1 --shoot --z_shift
python3 label_rsml.py /home/daniel/workspace/DUMUX/CPlantBox/gui/estimate/img/dicot/ 1 --split --z_shift
"""
#			1			2			3			4			5		6
scenario_name = [ "Constant scenario 1", "Constant scenario 2", "Dynamic scenario 1", "Dynamic scenario 2", "MRI_scenario_constant", "MRI_scenario_dynamic"]

def label_file(file, z_shift, artificial_shoot, split, scenario_index, container_volume):
    """ opens rsml and computes krs and suf """

    data = ViewerDataModel()
    data.open_rsml(file, z_shift)  # same as in gui (to check)
    if artificial_shoot:
        data.add_artificial_shoot()

    r = data.xylem_flux
    j = scenario_index - 1
    if j == 0:
        viewer_conductivities.init_constant_scenario1(r)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(r)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(r)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(r)
    elif j == 4:
        viewer_conductivities.init_MRI_paper_scenario_constant(r) # include both scenarios of the MRI Paper
    elif j == 5:
        viewer_conductivities.init_MRI_paper_scenario_dynamic(r)
    

    print("Scenario: "+ str(scenario_name[j]) +  " ,see file viewer_conductivities.py")    

    kr, kx = r.get_conductivities(-1.)  # (default) -1. means: current time is maximum of node creation times
    data.analyser.addData("kr", kr)
    data.analyser.addData("kx", kx)
    print("kr min max mean [cm]: " + str(np.min(kr)) +" "+ str(np.max (kr)) +" "+ str(np.mean(kr)))
    print("kx min max mean [cm]: " + str(np.min(kx)) +" "+ str(np.max (kx)) +" "+ str(np.mean(kx)))
    #print("\nkr")
    #print(kr)
    #print("\n")
    #print("\nkx")
    #print(kr)
    #print("\n")


    #######   Root Metrics:   ########
    #Root system length TRL (cm): 
    ana = pb.SegmentAnalyser(data.analyser) 
    segment_length = [ana.getSegmentLength(i) for i in range(0, len(ana.segments))]
    TRL = np.sum(segment_length)
    mean_seg_length = np.mean(segment_length)
    data.analyser.addData("segment_length", segment_length)
    print("TRL [cm]: ", TRL)

    ###Root length density (cm/cm3):
    #volume_container MRI: 455 cm³, Radius = 2.8 cm, h = 18.5 cm 
    #volume_container CT:  885 cm³,  Radius = 3.5 cm, h = 23.0 cm 
    convH_MRI = 455
    convH_CT  = 885
    if container_volume == "MRI":
       convH = convH_MRI
    elif container_volume == "CT": 
       convH = convH_CT
    RLD = TRL / convH
    print("RLD [cm/cm³]: ", RLD, "   with a convex hull of ", convH, "cm³")

    ###Half mean distance HMD (cm):
    HMD = (np.pi*RLD)**-0.5
    print("HMD [cm]: ", HMD)

    ###Number of roots / tips: 
    print("Number of roots: ", len(data.polylines)) 
    
    ###Number of tips per order:
    tips=r.get_orders_tips() # node IDs of root tips 
    order = []
    for i in tips:
        orders = data.types[i] # get type of node ID
        order.append(orders)
    order_list = Counter(order)
    print("Roots per order:", sorted(order_list.items()))

    ###Number of tips >= order 5:
    large_order = []
    for i in order:
       if i >=5:
          large_order.append(i)
    if len(large_order) == 0:
       print("Roots >= Order 5: 0")
    else:
       print("Roots >= Order 5: ", len(large_order), "   highest order = " , np.max(large_order))

    ###Segment Radii (cm)
    segRadii = data.radii
    
    
    weight_mean_radii= segment_length / TRL
    #print(weight_mean_radii)
    #print(np.sum(weight_mean_radii))
    radius = np.array(ana.getParameter("radius"))
    weighted_mean_radii = weight_mean_radii * radius

    print("radius min max mean [cm]: " + str(np.min(segRadii)) +" "+ str(np.max (segRadii)) +" "+ str(np.mean(segRadii)) +" "+ "weighted mean: "+ str(np.sum(weighted_mean_radii))) # Tobi



    ###Number of segments:    
    print("Number of segments: ", len(r.get_segments()))
    
    ###Length of segments:  
    print("mean segment length [cm]: ", mean_seg_length)  
    
    ###Creation-Time [d]:
    # final_age = np.max(cts), see xylem_flux.py, this means max ct defines the age of the root system 

    cts2 = np.array(ana.data["creationTime"])
    final_age2 = np.max(cts2)
    print("Root System age: ", np.max(cts2))
    node_ages = final_age2 * np.ones(cts2.shape) - cts2  # from creation time to age
    print("segment creation time min max mean [d]" + str(np.min(cts2)) +" "+ str(np.max (cts2)) +" "+ str(np.mean(cts2)))
    print("segment ages          min max mean [d]" + str(np.min(node_ages)) +" "+ str(np.max (node_ages)) +" "+ str(np.mean(node_ages)))


        


    if split:

        # ViewerDataModel
        # data.base_nodes  # base nodes indices (of roots or multiple plants)
        # data.base_segs  # emerging segment indices from base nodes

        # XylemFluxPython
        # r.neumann_ind  # node indices for Neumann flux
        # r.dirichlet_ind  # node indices for Dirichlet flux
        # print(data.base_segs)
        # print(data.base_nodes)

        assert len(data.base_segs) == len(data.base_nodes), "base segments must emerge from different nodes"
        n = len(data.base_segs)
        krs = [None] * n
        depth = [None] * n

        y = []
        for i in data.base_nodes:
            y.append(data.analyser.nodes[i].y)
        ysorted = np.argsort(y)
        print("y situation:", y, ysorted)

        n2 = int(np.ceil(-data.analyser.getMinBounds().z))
        suf_ = np.zeros((n, int(n2) * 10))

        for i in ysorted:
            r.dirichlet_ind = [ data.base_nodes[i] ]  # krs is based on dirichlet boundary condition
            krs[i], _ = r.get_krs(data.max_ct, [data.base_segs[i]])
            r.neumann_ind = [data.base_nodes[i]]  # suf is based on neumann boundary condititon
            suf = r.get_suf(data.max_ct)
            data.analyser.addData("SUF", suf)
            suf_[i,:] = data.analyser.distribution("SUF", 0., float(-n2), int(n2) * 10, False)  # mm layers
            depth[i] = r.get_mean_suf_depth(data.max_ct)  # mean depth is based on suf, i.e. based on neumann bc
    else:

        krs0, _ = r.get_krs(data.max_ct, data.base_segs)
        krs = [krs0]
        suf = r.get_suf(data.max_ct)

        rx_krs = r.get_rx_krs(data.max_ct, data.base_segs) # get xylem pressure of krs calculation
        print("mean xylem pressure resulting from krs calculation: " + str(np.mean(rx_krs)))
        rx_p_s = r.get_rx_krs_p_s(data.max_ct, data.base_segs) # get soil water potential (hydrostatic equilibrium) applied in krs calculation
        print("while a mean soil water potential (hydrostatic equilibrium of " + str(np.mean(rx_p_s)) +"is applied.")
        print("Transpiration", r.collar_flux(data.max_ct, rx_krs, rx_p_s, cells = False), "cm3/day") 
        # @param sxx [cm] (here rx_p_s) soil matric potentials given per segment or per soil cell
        data.analyser.addData("age", node_ages)
        data.analyser.addData("rx_krs", rx_krs,)
        data.analyser.addData("SUF", suf)
        data.analyser.write(str(file)+".vtp", ["rx_krs","subType", "radius", "creationTime", "age", "kr", "kx", "segment_length", "SUF"])
        #print("Transpiration", r.collar_flux(data.max_ct, rx_krs, [p_s]), "cm3/day")
    


        n = int(np.ceil(-data.analyser.getMinBounds().z))
        suf_ = np.zeros((1, int(n) * 10))
        suf_[0,:] = data.analyser.distribution("SUF", 0., float(-n), int(n) * 10, False)  # mm layers
        depth0 = r.get_mean_suf_depth(data.max_ct)
        depth = [depth0]

    print("\n" + "mean SUF depth (zSUF): ", depth, "cm")
    
    return suf_, krs, depth



def write_csv(file, suf_, krs, si, depth, split, ind = 0):
    """ writes an xls sheet containing suf """
    if split:
        file_csv = file.rsplit('.', 1)[0] + '_suf' + str(si) + '_' + str(ind) + '.csv'
    else:
        file_csv = file.rsplit('.', 1)[0] + '_suf' + str(si) + '.csv'
    print("krs {:g}, suf from {:g} to {:g}, sums up to {:g}".format(krs, np.min(suf_), np.max(suf_), np.sum(suf_)))
    suf_ = np.insert(suf_, 0, depth)
    suf_ = np.insert(suf_, 0, krs)
    df = pd.DataFrame(suf_)
    df.to_csv(file_csv, index = False, header = False)
    print("done. \n\n")


if __name__ == '__main__':

    """ parse arguments """
    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', type = str, help = 'file path')
    parser.add_argument('scenario_index', type = int, help = 'scenario index (1-4)')
    parser.add_argument('--shoot', action = 'store_true', help = 'adds an artificial shoot')
    parser.add_argument('--split', action = 'store_true', help = 'splits output for multiple plants into different csv files')
    parser.add_argument('--z_shift', action = 'store_true', help = 'shifts seed to -3cm (based on the first seed, if multiple plants are present)')
    parser.add_argument('container_volume', type = str, help = 'MRI or CT')
    args = parser.parse_args()

    walk_dir = args.file_path
    print('walk_dir = ' + walk_dir)

    # If your current working directory may change during script execution, it's recommended to
    # immediately convert program arguments to an absolute path. Then the variable root below will
    # be an absolute path as well. Example:
    # walk_dir = os.path.abspath(walk_dir)
    print('walk_dir (absolute) = ' + os.path.abspath(walk_dir))
    artificial_shoot = args.shoot
    split = args.split
    z_shift = args.z_shift
    container_volume = args.container_volume

    scenario_index = args.scenario_index
    print("Scenario index {:g}, see file viewer_conductivities.py".format(scenario_index))

    if scenario_index < 1 or scenario_index > 6:
        raise

    print()

    """ walk the files """
    

    for root, subdirs, files in os.walk(walk_dir):
        sorted_files =  sorted(files) #list files to print console outputs in correct order
        sys.stdout = open("krs_zSUF_" + str(scenario_name[scenario_index-1]) + ".txt", "w")
        # print('--\nroot = ' + root)
        list_file_path = os.path.join(root, 'my-directory-list.txt')
        # print('list_file_path = ' + list_file_path)

        with open(list_file_path, 'wb') as list_file:

            for subdir in subdirs:
                print('\t- subdirectory ' + subdir)

            for filename in sorted_files: #listed files

                if filename.endswith('.rsml'):

                    file_path = os.path.join(root, filename)
                    print('file %s (full path: %s)\n' % (filename, file_path))

                    suf_, krs, depth = label_file(file_path, z_shift, artificial_shoot, split, scenario_index, container_volume)
                    if len(krs) > 1:
                        for i in range(0, len(krs)):
                            write_csv(file_path, suf_[i,:], krs[i], scenario_index, depth[i], split, i)
                    else:
                        write_csv(file_path, suf_[0,:], krs[0], scenario_index, depth[0], split)
    sys.stdout.close()

