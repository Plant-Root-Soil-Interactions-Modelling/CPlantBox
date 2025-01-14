"""root system surface density"""
import sys; sys.path.append("../..");
sys.path.append("../../src/")
import plantbox as pb
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functional.xylem_flux import XylemFluxPython

font = {'size'   : 16}
plt.rc('font', **font)

def maize_conductivities(r, skr = 1., skx = 1.):
    """ Hydraulic conductivities for maize modified from Couvreur et al. (2012), originally from Doussan et al. (1998) """

    kr00 = np.array([[0., 0.]])  # artificial shoot
    kx00 = np.array([[0., 1.e3]])  # artificial shoot

    kr0 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [8., 0.000181], [10, 0.0000648], [18, 0.0000648], [25, 0.0000173], [300, 0.0000173]])
    kr1 = np.array([[-1e4, 0.], [-0.1, 0.], [0., 0.000181], [10., 0.000181], [16, 0.0000173], [300, 0.0000173]])

    kx0 = np.array([[0., 0.000864], [5., 0.00173], [12., 0.0295], [15., 0.0295], [20., 0.432], [300., 0.432]])
    kx1 = np.array([[0., 0.000864], [5., 0.0000864], [10., 0.0000864], [12., 0.0006048], [20., 0.0006048], [23., 0.00173], [300., 0.00173]])

    kr01 = np.minimum(skr * kr0[:, 1], 1.)  # kr0[:, 1] are values
    kr11 = np.minimum(skr * kr1[:, 1], 1.)  # kr1[:, 1] are values
    r.setKrTables([kr00[:, 1], kr01, kr11, kr11, kr01, kr01],
                  [kr00[:, 0], kr0[:, 0], kr1[:, 0], kr1[:, 0], kr0[:, 0], kr0[:, 0]])
    kx01 = np.minimum(skx * kx0[:, 1], 1.)
    kx11 = np.minimum(skx * kx1[:, 1], 1.)
    r.setKxTables([kx00[:, 1], kx01, kx11, kx11, kx01, kx01],
                  [kx00[:, 0], kx0[:, 0], kx1[:, 0], kx1[:, 0], kx0[:, 0], kx0[:, 0]])

#simulation 
path = "rootsystem/"
names = ["RS_optimized_field_L_WT_new",  "RS_optimized_field_S_WT_new"]
title = ["Loam", "Sand"]
cols = ['r','b','g','c']
depth = 80
layers = 9
times = [42, 63, 98, 154]
BBCHs = [14,19,59,83]


fig1, ax1 = plt.subplots(1, 2)
fig2, ax2 = plt.subplots()
xx = []
yy = []
for i in range(0,len(names)):
    suf_ = []
    dummy = 0
    name = names[i]
    #rs = pb.RootSystem()
    rs = pb.MappedRootSystem()
    rs.readParameters(path + name + ".xml")

    rs.initializeLB(5,4)
    krs = []
    for j in range(0,times[-1]+1):
        if i ==0:
            if j<=98:
                rs.simulate(1, False)
        elif i == 1:
            rs.simulate(1, False)

        if j in times:
            r = XylemFluxPython(rs)  # hydraulic model
            maize_conductivities(r)
            suf = r.get_suf(j)
            krs_, _ = r.get_krs(j)
            krs.append(krs_)
            print("Krs: ", krs, "cm2/day")

            ana = pb.SegmentAnalyser(rs.mappedSegments())
            ana.addData("SUF", suf)
            n = int(np.ceil(-ana.getMinBounds().z))
            #z_ = np.linspace(-0.5, -n + 0.5, layers)
            z_ = np.linspace(0, -depth, layers)  # z - axis
            d = ana.distribution("SUF", 0., float(-n), int(layers), False)  # False!!!
            suf_.append(d)
            ax1[i].plot(d, z_, color = cols[dummy],label = "BBCH "+str(BBCHs[dummy]) )
            ax1[i].title.set_text(title[i])
            ax1[i].set_xlabel("SUF $(-)$")
            ax1[i].set_ylabel("Depth $(cm)$")
            ax1[i].set_ylim([-75,0])
            if i == 1: 
                ax1[i].legend()
            dummy = dummy+1


        if j == times[-1]:
            ax2.plot(times, krs, marker = '*', color = cols[i], label = title[i]) 
            ax2.set_xlabel("Time $(d)$")
            ax2.set_ylabel("Root system conductivity krs $(cm^2$ $d^{-1})$")
            ax2.legend()

    #krs for csv
    yy.extend(np.transpose(np.vstack((np.repeat(i,len(krs)),BBCHs, krs))))

    #suf for csv
    xx.extend(np.transpose(np.vstack((np.repeat(i,len(z_)),z_,suf_))))

np.savetxt('krs.csv', yy, delimiter = ",", header="treatment,BBCH,krs",fmt="%0.4f")
np.savetxt('SUF.csv', xx, delimiter = ",", header="treatment,depth,BBCH14,BBCH19,BBCH59,BBCH83",fmt="%0.4f")

plt.show()


