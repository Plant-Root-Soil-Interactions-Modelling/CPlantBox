import sys;  sys.path.append("../../src/python_modules/")
from vtk_tools import *
import os
import numpy as np
from io import StringIO  # StringIO behaves like a file object
import math


def read_rootsys(fname):

    with open(fname) as f:

        content = f.readlines()

        # read relevant data
        simtime = float(content[1])

        table1 = ""  # holds: segID#    x          y          z      prev or  br#  length   surface  mass
        table2 = ""  # holds: origination time

        i = 26  # start row
        while i < len(content):
            line = content[i]
            if len(line) < 40:  # lousy stopping criteria
                break
            table1 += (line + "\n")
            table2 += (content[i + 1] + "\n")
            i += 2

        id2, x, y, z, prev, order, bn, length, surface , mass = np.loadtxt(StringIO(table1), unpack=True)
        ctime = np.loadtxt(StringIO(table2))  #

        for i, id_ in enumerate(id2):
            if i != int(id_ - 1):
                print("ERROR: ids are not sequential, this is not implemented (yet)", i, id_ - 1)
                return

        # create polydata object
        nodes = np.zeros((len(x), 3))
        nodes[:, 0] = x
        nodes[:, 1] = y
        nodes[:, 2] = z
        points = vtk_points(nodes)

        bn_ = vtk.vtkIntArray()
        order_ = vtk.vtkIntArray()
        length_ = vtk.vtkFloatArray()
        surface_ = vtk.vtkFloatArray()
        radius_ = vtk.vtkFloatArray()
        mass_ = vtk.vtkFloatArray()
        time_ = vtk.vtkFloatArray()
        age_ = vtk.vtkFloatArray()

        order_.SetName("type")  # reader expects that name
        bn_.SetName("branch number")        
        length_.SetName("length")
        radius_.SetName("radius")
        mass_.SetName("mass")
        time_.SetName("node_creation_time")  # reader expects that name
        age_.SetName("age")
        surface_.SetName("surface")

        segs = []
        for i in range(0, len(id2)):
            if prev[i] != 0:
                i_ = int(id2[i] - 1)  # lets start at 0
                segs.append((int(prev[i] - 1), i_))
                bn_.InsertNextValue(int(bn[i_]))
                order_.InsertNextValue(int(order[i_]))  # used -7 for RootSys_verysimple
                length_.InsertNextValue(length[i_])
                surface_.InsertNextValue(surface[i_])
                radius_.InsertNextValue(surface[i_] / length[i_] / (2 * math.pi))
                mass_.InsertNextValue(mass[i_])
                time_.InsertNextValue(ctime[i_])
                age_.InsertNextValue(simtime - ctime[i_])
        
        segs = np.array(segs, dtype=int)  # copy list into numpy array
        cells = vtk_cells(segs)

        pd = vtk.vtkPolyData()
        pd.SetPoints(points)
        pd.SetLines(cells)
        pd.GetCellData().AddArray(order_)
        pd.GetCellData().AddArray(bn_)
        # pd.GetCellData().AddArray(surface_)
        # pd.GetCellData().AddArray(mass_)
        pd.GetCellData().AddArray(radius_)
        pd.GetCellData().AddArray(length_)
        pd.GetCellData().AddArray(time_)
        # pd.GetCellData().AddArray(age_)

    return pd


if __name__ == "__main__":

    name = "sys_files/RootSys_verysimple"
    pd = read_rootsys(name)
    
    write_rsml(name + ".rsml", pd, 1)  # index 0 == root order
    
    print("fin")

