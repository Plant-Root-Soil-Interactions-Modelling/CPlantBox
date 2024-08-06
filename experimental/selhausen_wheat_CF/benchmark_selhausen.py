"""" CF (wheat) Selhausen Rhizotron facility """
import sys
sys.path.append("../..")
sys.path.append("../../src")
import plantbox as pb
import plantbox as pb
import math
import numpy as np
import pickle as pkl
import pandas as pd
import matplotlib.pyplot as plt

name = "wheat"

simtime = 240

M = 35  # number of plants in rows
N = 40  # number of rows

distp = 3  # distance between the root systems along row[cm] (P-P)
distr = 12  # distance between the rows[cm] (R-R)
interrow = M * distp  # intra-row spacing
row = N * distr  # row spacing

depth = 120  # TODO 130
layers = 12  # ## 13
soilvolume = (depth / layers) * interrow * row

box = pb.SDF_PlantBox(5000, 5000, 5000)  # box
rhizotube = pb.SDF_PlantContainer(6.4 / 2, 6.4 / 2, 700, False)  # a single rhizotube
rhizoX = pb.SDF_RotateTranslate(rhizotube, 90, pb.SDF_Axis.yaxis, pb.Vector3d(700 / 2, 0, 0))
rhizotubes_ = []
y_ = (-30, -18, -6, 6, 18, 30)
z_ = (-10, -20, -40, -60, -80, -120)
tube = []
for i in range(0, len(y_)):
    v = pb.Vector3d(0, y_[i], z_[i])
    tube.append(pb.SDF_RotateTranslate(rhizoX, v))
    rhizotubes_.append(tube[i])

rhizotubes = pb.SDF_Union(rhizotubes_)

rhizotube_mirror1 = pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, distp * M, 0))
rhizotube_mirror2 = pb.SDF_RotateTranslate(rhizotubes, pb.Vector3d(0, -distp * M, 0))
rhizotuben = pb.SDF_Union([rhizotube_mirror1, rhizotubes, rhizotube_mirror2])

v_d = 1.24  # viewing_depth in mm
# External geom
rhizotubeo = pb.SDF_PlantContainer((6.4 / 2) + ((v_d * 2) / 10), (6.4 / 2) + ((v_d * 2) / 10), 700, False)
rhizoXo = pb.SDF_RotateTranslate(rhizotubeo, 90, pb.SDF_Axis.yaxis, pb.Vector3d(700 / 2, 0, 0))

rhizotubes_o = []
tubeo = []
for i in range(0, len(y_)):
    vo = pb.Vector3d(0, y_[i], z_[i])
    tubeo.append(pb.SDF_RotateTranslate(rhizoXo, vo))
    rhizotubes_o.append(tubeo[i])

rhizotubeso = pb.SDF_Union(rhizotubes_o)
rhizoTube = pb.SDF_Difference(box, rhizotuben)

opening = pb.SDF_PlantBox(6, 4, 10)  # image size inhouse facility box

r0 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[0], z_[0]))
l0 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[0], z_[0]))
0,
r1 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[1], z_[1]))
l1 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[1], z_[1]))

r2 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[2], z_[2]))
l2 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[2], z_[2]))

r3 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[3], z_[3]))
l3 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[3], z_[3]))

r4 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[4], z_[4]))
l4 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[4], z_[4]))

r5 = pb.SDF_RotateTranslate(opening, 180 - 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[5], z_[5]))
l5 = pb.SDF_RotateTranslate(opening, 180 + 80, pb.SDF_Axis.xaxis, pb.Vector3d(0, y_[5], z_[5]))

lp = []  # list of image positions
for i in range(4):
    pim = np.arange(210 + (100 * i), 240 + (100 * i), 6)
    lp.append(pim)
x_ = np.asarray(lp).flatten() - (700 / 2)  # array([-140., -134., -128., -122., -116.,  -40.,  -34.,  -28.,  -22., -16.,   60.,   66.,   72.,   78.,   84.,  160.,  166.,  172., 178.,  184.])


def run_benchmark():  # make a function to be called certain times
    # Initializes N*M root systems
    allRS = []
    for i in range(0, M):
        for j in range(0, N):
            rs = pb.RootSystem()
            rs.readParameters(name + ".xml")
            rs.getRootRandomParameter(1).lmax = 180  ##################
            rs.getRootSystemParameter().seedPos = pb.Vector3d(distr * j - ((N - 1) * distr / 2), distp * i - (distp * M / 2), -3.)  # wheat rows perpendicular to tubes
            rs.setGeometry(rhizoTube)
            rs.initialize(False)
            # Simulate
            rs.simulate(simtime, False)
            allRS.append(rs)
            if i + j == 0:
                allAna1 = pb.SegmentAnalyser(rs)
            else:
                allAna1.addSegments(rs)

    # allAna1.write("all_plants.vtp")
    allAna1.mapPeriodic(row, interrow)
    # allAna1.write("mp.vtp")

    rl_ = []
    rl_.append(allAna1.distribution("length", 0., -depth, layers, True))
    rl_ = np.array(rl_[-1]) / soilvolume  # convert to density

    allAna = pb.SegmentAnalyser(allAna1)
    allAna.crop(pb.SDF_Difference(rhizotubeso, rhizotubes))
    # allAna.write("hc.vtp")
    # allAna.write("results/cg/hc.vtp")

    sa = 6 * 4  # surface_area for surface density

    ls = [r0, l0, r1, l1, r2, l2, r3, l3, r4, l4, r5, l5]

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
    rsd_a = np.array([np.mean(t1), np.mean(t2), np.mean(t3), np.mean(t4), np.mean(t5), np.mean(t6)])  # planar root length densities
    res_arr = np.zeros((rsd_a.size, 6))

    # rl_[0] : volumetric rld in [0, -10] cm [0, -20]
    # rl_[1] : volumetric rld in [-10, -20]
    # rl_[2] : volumetric rld in [-20, -30]
    # rl_[3] : volumetric rld in [-30, -40]
    vrld_t1 = mean(rl_[0], rl[1])  #  volumetric rld in [0, -20]  # tube located at -10
    vrld_t2 = mean(rl_[1], rl[2])  #  volumetric rld in [-10, -30] # tube located at -20
    vrld_t2 = mean(rl_[3], rl[4])  #  volumetric rld in [-30, -50] # tube located at -40

    mrld = np.array([ rl_[1 - 1], rl_[2 - 1], rl_[4 - 1], rl_[6 - 1], rl_[8 - 1], rl_[12 - 1] ])  # average root length density below one plant in the crop stand

    # OTHER approach
    # box [0, -15] # tube -10, second z coordinate in between tube locations
    # box [-15, -30] # tube -20
    # box [-30, -50] # tube -40
    # # and so on

    res_arr[:, 0] = z_
    res_arr[:, 1] = rsd_a  # pRLD
    res_arr[:, 2] = mrld
    res_arr[:, 3] = mrld / rsd_a  # CF = vRLD/pRLD.
    res_arr[:, 4] = rsd_a * sa  # RL_image
    res_arr[:, 5] = (rsd_a * sa) / mrld  # CF(RL_image/vRLD)
    return res_arr


if __name__ == '__main__':  # for testing
    import time
    start_time = time.time()
    rlds = run_benchmark()
    print(rlds)
    print("--- %s seconds, end benchmark ---" % (time.time() - start_time))
