import sys; sys.path.append("../src/python_modules/"); sys.path.append("../")

import os
import pandas as pd
import numpy as np

import viewer_data
import xylem_flux


def label_file(file, artificial_shoot,):
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
    r.seg_ind = node_ind

    krs, _ = r.get_krs(data.max_ct)
    nop = len(node_ind)  # number of plants (we might want to multiply suf by it?)
    suf = r.get_suf(max_ct) * nop

    data.analyser.addData("SUF", suf)
    n = int(np.ceil(-analyser.getMinBounds().z))
    d = analyser.distribution("SUF", 0., float(-n), int(n), True)
    # z_ = np.linspace(-0.5, -n + 0.5, n)
    # ax.plot(d, z_, "-*", label = "total")
    return suf, krs


def write_xls(file, suf, krs):
    """ writes an xls sheet containing suf """
    file_xls = file.rsplit('.', 1)[0] + '.xls'
    # df1 = pd.DataFrame(np.transpose(np.array(psi_x_)))
    # df1.to_excel(file1, index = False, header = False)


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
        print('--\nroot = ' + root)
        list_file_path = os.path.join(root, 'my-directory-list.txt')
        print('list_file_path = ' + list_file_path)

        with open(list_file_path, 'wb') as list_file:

            for subdir in subdirs:
                print('\t- subdirectory ' + subdir)

            for filename in files:

                if filename.endswith('.rsml'):

                    file_path = os.path.join(root, filename)
                    print('\t- file %s (full path: %s)' % (filename, file_path))

                    label_file(file_path, artificial_shoot, scenario_index)

