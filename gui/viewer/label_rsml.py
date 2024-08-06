import sys; sys.path.append("../.."); sys.path.append("../../src/")

import os
import pandas as pd
import numpy as np
import argparse

import viewer_conductivities
from viewer_data import ViewerDataModel
import functional.xylem_flux as xylem_flux

"""
creates a csv file per rsml file containing krs values, and suf values per 1 mm layers

expected three arguments: file path, scenario index (1-4), 
optionally, --shoot, only use for measured monocots; --split, for multiple dicots, --z_shift to relocate seed to -3cm

e.g. 
python3 label_rsml.py ~/Downloads/second+round/monocot/maize 1
python3 label_rsml.py ~/Downloads/second+round/dicot/lupin 1

python3 label_rsml.py Benchmarking_ 1
python3 label_rsml.py ~/Downloads/second+round/dicot/lupin 1

python3 label_rsml.py ~/Downloads/second+round/dicot/lupin 1

python3 label_rsml.py /home/daniel/workspace/DUMUX/CPlantBox/gui/estimate/img/monocot/ 1 --shoot --z_shift
python3 label_rsml.py /home/daniel/workspace/DUMUX/CPlantBox/gui/estimate/img/dicot/ 1 --split --z_shift

python3 label_rsml.py Files_for_Daniel/lupin/ 3 --z_shift
python3 label_rsml.py Files_for_Daniel/maize/archisimple/ 3 --shoot --z_shift
python3 label_rsml.py Files_for_Daniel/maize/roottyp/ 3 --z_shift

results for Benjamin Delory from 6.12.2023:
python3 label_rsml.py Archive-2023-11-30/maize/ 3 --shoot --z_shift
python3 label_rsml.py Archive-2023-11-30/lupin/ 3 --z_shift
python3 label_rsml.py Archive-2023-11-30/reference/dicot/lupin 3 --split --z_shift
python3 label_rsml.py Archive-2023-11-30/reference/monocot/maize 3 --shoot --z_shift
"""


def label_file(file, z_shift, artificial_shoot, split, scenario_index):
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
        data.analyser.addData("SUF", suf)
        n = int(np.ceil(-data.analyser.getMinBounds().z))
        suf_ = np.zeros((1, int(n) * 10))
        suf_[0,:] = data.analyser.distribution("SUF", 0., float(-n), int(n) * 10, False)  # mm layers
        depth0 = r.get_mean_suf_depth(data.max_ct)
        depth = [depth0]

    print("depth: ", depth, "cm")
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
    parser.add_argument('--z_shift', action = 'store_true', help = 'shifts seed to -3 cm (based on the first seed, if multiple plants are present)')
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

    scenario_index = args.scenario_index
    print("Scenario index {:g}, see file viewer_conductivities.py".format(scenario_index))

    if scenario_index < 1 or scenario_index > 4:
        raise

    print()

    """ walk the files """
    for root, subdirs, files in os.walk(walk_dir):
        # print('--\nroot = ' + root)
        list_file_path = os.path.join(root, 'my-directory-list.txt')
        # print('list_file_path = ' + list_file_path)

        with open(list_file_path, 'wb') as list_file:

            for subdir in subdirs:
                print('\t- subdirectory ' + subdir)

            for filename in files:

                if filename.endswith('.rsml'):

                    file_path = os.path.join(root, filename)
                    print('file %s (full path: %s)\n' % (filename, file_path))

                    suf_, krs, depth = label_file(file_path, z_shift, artificial_shoot, split, scenario_index)
                    if len(krs) > 1:
                        for i in range(0, len(krs)):
                            write_csv(file_path, suf_[i,:], krs[i], scenario_index, depth[i], split, i)
                    else:
                        write_csv(file_path, suf_[0,:], krs[0], scenario_index, depth[0], split)

