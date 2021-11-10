import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import os
import pandas as pd
import numpy as np

import viewer_conductivities
import viewer_data
import xylem_flux

"""
creates a csv file per rsml file containing krs values, and suf values per 1 mm layers

expected three arguments: file path, artifical shoot (true/false), scenario index (1-4)

e.g. 
python3 label_rsml.py ~/workspace/DUMUX/CPlantBox/gui/maize true 1
python3 label_rsml.py ~/workspace/DUMUX/CPlantBox/gui/dicot false 1

"""


def label_file(file, artificial_shoot, scenario_index):
    """ opens rsml and computes krs and suf """
    data = viewer_data.DataModel()
    data.open_rsml(file)  # same as in gui (to check)
    if artificial_shoot:
        data.add_artificial_shoot()

    r = xylem_flux.XylemFluxPython(data.mapped_segments)
    j = scenario_index - 1
    if j == 0:
        viewer_conductivities.init_constant_scenario1(r)
    elif j == 1:
        viewer_conductivities.init_constant_scenario2(r)
    elif j == 2:
        viewer_conductivities.init_dynamic_scenario1(r)
    elif j == 3:
        viewer_conductivities.init_dynamic_scenario2(r)
    node_ind = data.get_base_node_indices()
    r.seg_ind = node_ind
    krs, _ = r.get_krs(data.max_ct)
    nop = len(node_ind)  # number of plants (we might want to multiply suf by it?)
    nop = 1
    suf = r.get_suf(data.max_ct) * nop
    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-data.analyser.getMinBounds().z))
    suf_ = data.analyser.distribution("SUF", 0., float(-n), int(n) * 10, False)  # mm layers
    return suf_, krs


def write_csv(file, suf_, krs, si):
    """ writes an xls sheet containing suf """
    file_csv = file.rsplit('.', 1)[0] + '_suf' + str(si) + '.csv'
    print("krs {:g}, suf from {:g} to {:g}, sums up to {:g}".format(krs, np.min(suf_), np.max(suf_), np.sum(suf_)))
    suf_ = np.insert(suf_, 0, krs)
    df = pd.DataFrame(suf_)
    df.to_csv(file_csv, index = False, header = False)
    print("done. \n")


if __name__ == '__main__':

    """ parse arguments """
    try:
        walk_dir = sys.argv[1]
        print('walk_dir = ' + walk_dir)
        # If your current working directory may change during script execution, it's recommended to
        # immediately convert program arguments to an absolute path. Then the variable root below will
        # be an absolute path as well. Example:
        # walk_dir = os.path.abspath(walk_dir)
        print('walk_dir (absolute) = ' + os.path.abspath(walk_dir))

        artificial_shoot = bool(sys.argv[2])
        print("Artifical shoot", str(artificial_shoot))

        scenario_index = int(sys.argv[3])
        print("Scenario index {:g}, see file viewer_conductivities.py".format(scenario_index))

        if scenario_index < 1 or scenario_index > 4:
            raise
    except:
        print("Expected three arguments: file path, artifical shoot (true/false), scenario index (1-4)")
        sys.exit(0)

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
                    print('\t- file %s (full path: %s)' % (filename, file_path))

                    suf_, krs = label_file(file_path, artificial_shoot, scenario_index)
                    write_csv(file_path, suf_, krs, scenario_index)

