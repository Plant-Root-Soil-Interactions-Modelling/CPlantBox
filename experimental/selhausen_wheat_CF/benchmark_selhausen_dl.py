"""" CF (wheat) Selhausen Rhizotron facility """
import sys; sys.path.append("../.."); sys.path.append("../../src")

import plantbox as pb
import numpy as np
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt
import time

name = "wheat"

simtime = 30  # 240

M = 35  # number of plants in rows
N = 40  # number of rows

distp = 3  # [cm] inter-plant distance, i.e. between the root systems along row (P-P)
distr = 12  # [cm] inter-row distance, i.e. distance between the rows (R-R)

width_ = M * distp  # field width
length_ = N * distr  # field length

depth = 120  # [cm]
layers = 12  #
soilvolume = (depth / layers) * width_ * length_  # [cm3]

# tube positions
y_ = (-30, -18, -6, 6, 18, 30)  # [cm]
z_ = (-10, -20, -40, -60, -80, -120)  # [cm]

tube_diameter = 6.4  # [cm]
tube_length = 700  # [cm]

v_d = 1.24  # [mm] camera viewing_depth

""" create rhizotubes """
box = pb.SDF_PlantContainer(100., 100., 700., True)  # box
rhizotube = pb.SDF_PlantContainer(tube_diameter / 2, tube_diameter / 2, tube_length, False)  # a single rhizotube
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(tube_length / 2, 0, 0))
rhizotubes_ = []
tube = []
for i in range(0, len(y_)):
    v = pb.Vector3d(0, y_[i], z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_)

rhizotube_mirror1 = pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, distp * M, 0))
rhizotube_mirror2 = pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, -distp * M, 0))
rhizotuben = pb.SDF_Union([rhizotube_mirror1, rhizotubes, rhizotube_mirror2])

# External geom (including viewing depth)
rhizotubeo = pb.SDF_PlantContainer((tube_diameter / 2) + (v_d / 10), (tube_diameter / 2) + (v_d / 10), 700, False)  # TODO in original file there was a factor *2 ?????
rhizoXo = pb.SDF_RotateTranslate(rhizotubeo, 90, pb.SDF_Axis.yaxis, pb.Vector3d(tube_length / 2, 0, 0))

rhizotubes_o = []
tubeo = []
for i in range(0, len(y_)):
    vo = pb.Vector3d(0, y_[i], z_[i])
    tubeo.append(pb.SDF_RotateTranslate(rhizoXo, vo))
    rhizotubes_o.append(tubeo[i])

rhizotubeso = pb.SDF_Union(rhizotubes_o)
rhizoTube = pb.SDF_Difference(box, rhizotuben)

opening = pb.SDF_PlantBox(6, 4, 10)  # image size inhouse facility box

ls = []
for i in range(0, 6):
    ls.append(pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[i], z_[i])))
    ls.append(pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[i], z_[i])))

lp = []  # list of image positions
for i in range(4):
    pim = np.arange(210 + (100 * i), 240 + (100 * i), 6)
    lp.append(pim)

x_ = np.asarray(lp).flatten() - (700 / 2)  # array([-140., -134., -128., -122., -116.,  -40.,  -34.,  -28.,  -22., -16.,   60.,   66.,   72.,   78.,   84.,  160.,  166.,  172., 178.,  184.])


def run_benchmark():  # make a function to be called certain times

    """ geometry """
    rs = pb.RootSystem()
    rs.setGeometry(rhizoTube)
    rs.write("rhizoTube.py")

    # Initializes N*M root systems
    allRS = []
    allAna1 = pb.SegmentAnalyser()
    for i in range(0, M):  # number of plants in rows
        for j in range(0, N):  # nubmer of rows

            rs = pb.RootSystem()
            rs.readParameters(name + ".xml")

            # print(rs.getRootRandomParameter(1).lmax)
            # rs.getRootRandomParameter(1).lmax = 180  # TODO is that right?
            # print(rs.getRootRandomParameter(1).lmax)
            # dd

            rs.getRootSystemParameter().seedPos = pb.Vector3d(distr * j - ((N - 1) * distr / 2), distp * i - (distp * M / 2), -3.)  # wheat rows perpendicular to tubes
            rs.setGeometry(rhizoTube)  # ## <----- ?????
            rs.initialize(False)
            # Simulate
            rs.simulate(simtime, False)
            allRS.append(rs)
            allAna1.addSegments(rs)

    allAna1.mapPeriodic(length_, width_)
    # allAna1.write("all_plants.vtp")

    # rl_ = []
    # rl_.append(allAna1.distribution("length", 0., -depth, layers, True))
    # rl_ = np.array(rl_[-1]) / soilvolume  # convert to density
    rl_ = np.array(allAna1.distribution("length", 0., -depth, layers, True)) / soilvolume  # convert to density

    allAna = pb.SegmentAnalyser(allAna1)
    allAna.crop(pb.SDF_Difference(rhizotubeso, rhizotubes))
    allAna.pack()
    # allAna.write("hc.vtp")
    # allAna.write("results/cg/hc.vtp")

    sa = 6 * 4  # surface_area for surface density

    ls_not = ['r0_', 'l0_', 'r1_', 'l1_', 'r2_', 'l2_', 'r3_', 'l3_', 'r4_', 'l4_', 'r5_', 'l5_']  # list of notations
    es = np.zeros((x_.size, len(ls)))  # (20 positions of images in one side in one tube, 12 sides of all tubes)

    for j in range(len(ls)):
        for i in range(es.shape[0]):
            v = pb.Vector3d(x_[i], 0, 0)
            g = pb.SDF_RotateTranslate(ls[j], v)
            ana = pb.SegmentAnalyser(allAna)
            ana.crop(g)
            ana.pack()
            rsd_each = (ana.getSummed("length")) / sa  # root surface density

            es[i, j] = rsd_each

    t1 = es[:,:2]
    t2 = es[:, 2:4]
    t3 = es[:, 4:6]
    t4 = es[:, 6:8]
    t5 = es[:, 8:10]
    t6 = es[:, 10:]
    rsd_a = np.array([np.mean(t1), np.mean(t2), np.mean(t3), np.mean(t4), np.mean(t5), np.mean(t6)])

    res_arr = np.zeros((rsd_a.size, 6))
    mrld = np.array([ rl_[1 - 1], rl_[2 - 1], rl_[4 - 1], rl_[6 - 1], rl_[8 - 1], rl_[12 - 1] ])  # average root length density below one plant in the crop stand
    res_arr[:, 0] = z_
    res_arr[:, 1] = rsd_a  # pRLD
    res_arr[:, 2] = mrld
    res_arr[:, 3] = mrld / rsd_a  # CF = vRLD/pRLD.
    res_arr[:, 4] = rsd_a * sa  # RL_image
    res_arr[:, 5] = (rsd_a * sa) / mrld  # CF(RL_image/vRLD)

    return res_arr


if __name__ == '__main__':  # for testing

    start_time = time.time()
    rlds = run_benchmark()
    print(rlds)
    print("\n--- %s seconds, end benchmark ---" % (time.time() - start_time))
