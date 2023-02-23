""" 
soil uptake fraction of a root system (soil is in hydrostatic equilibrium) 
"""
import sys; sys.path.append("../../python/modules/"); sys.path.append("../../../CPlantBox/"); 
sys.path.append("../../../CPlantBox/src/python_modules")

from xylem_flux import XylemFluxPython  # Python hybrid solver
import plantbox as pb
import vtk_plot as vp

import numpy as np
import matplotlib.pyplot as plt


# sets all standard deviation to a percantage of the mean value, i.e. value*s
def set_all_sd(rs, s):
    for p in rs.getRootRandomParameter():
        p.lmaxs = p.lmax * s
        p.lbs = p.lb * s
        p.las = p.la * s
        p.lns = p.ln * s
        p.rs = p.r * s
        p.a_s = p.a * s
        p.thetas = p.theta * s        
        p.dx = 0.1  # cm


""" Parameters """
name = "single_root"  # single_root, Zea_mays_1_Leitner_2010

simtime = 60  # days
kx_ = [4.e3, 4.e2, 4e1, 4, 4.e-1, 4.e-2, 4.e-3]  # 4000, 4, 0.004
kr = 0.002 

cols = ["-b*", "-g*", "-r*", "-m*", "-c*", "-y*", "b:*", "g:*" ]
labels = ["kx = {:g}".format(kx) for kx in kx_]

suf_ = []
for ii, kx in enumerate(kx_):
    
    print(ii, kx, cols[ii])
    suf_.append([])

    # artificial shoot
    kr0 = np.array([[-154, 0.], [0, 1.e-12], [1e20, 1.e-12]])
    kx0 = 1.  # p.array([[-154, 0.], [0, 1.], [1e20, 1.]])    
    # const
    kr1 = np.array([[-1e4, 1.e-12], [0., 1.e-12], [0, kr], [100, kr]])  # negative values must be samll for kr
    kr2 = kr1
    kr3 = kr2    
#     kr0 = kr1  # disable special artificial root conductivities
#     kx0 = kx
    
    """ root system """
    rs = pb.MappedRootSystem()    
    rs.openFile(name, "")
    # rs.setSeed(2)  # fix randomness
    set_all_sd(rs, 0)  # set all std to zero, and dx = 0.1 cm
    rs.getRootSystemParameter().seedPos.z = -3   
    rs.initialize()
    rs.simulate(simtime, False)  # simulate all, then use age dependent conductivities for predefined growth
    
#     # Plot, using vtk
#     vp.plot_roots(rs, "creationTime")    
    
    """ set up xylem parameters """
    r = XylemFluxPython(rs)
    kr4 = kr1  # basal
    kr5 = kr1  # shoot borne
    print()
    r.setKrTables([kr0[:, 1], kr1[:, 1], kr2[:, 1], kr3[:, 1], kr4[:, 1], kr5[:, 1]], [kr0[:, 0], kr1[:, 0], kr2[:, 0], kr3[:, 0], kr4[:, 0], kr5[:, 0]])
    r.setKx([kx0, kx, kx, kx, kx, kx])
    
#     """ for debugging """
#     r.test()
#     r.plot_conductivities()
#     shoot_segs = rs.getShootSegments()
#     print("Shoot segments", [str(s) for s in shoot_segs])
#     print("Shoot type", rs.subTypes[0])
    
    """ numerical solution of transpiration -1 cm3/day"""
    krs_, l_, jc_ = [], [], []
    for t in range(10, simtime + 1):
        
        # l_.append(rs.getSummed("length")) # does not work, since rs is precomputed
        ana = pb.SegmentAnalyser(rs)
        ana.filter("creationTime", -1, t + 1.e-4)
        l_.append(ana.getSummed("length"))  # does not work, all rs is precomputedana.getParameter("length")
             
        suf = r.get_suf(t)   
        
        krs, jc = r.get_krs(t)    
        krs_.append(krs)
        jc_.append(jc)
        
        if t > 10 and t % 10 == 0:  # for the suf plot
            suf_[-1].append(suf)
    
    time = np.linspace(10, simtime, simtime - 10 + 1)
    # plt.plot(time, krs_, cols[ii], label=labels[ii], alpha=0.7)
    plt.plot(time, jc_, cols[ii], label=labels[ii])

# """ Krs plot """
plt.xlabel("simulation time")
# plt.ylabel("root system conductance $krs$ $[cm^2 d^{-1}]$")
plt.ylabel("Transpiration at eswp -500 cm $[cm^3/day]$")
plt.legend()
plt.show()

""" Root system length """
plt.xlabel("simulation time")
plt.ylabel("root system length $[cm]$")
plt.plot(time, l_)
plt.show()

print("Final length", l_[-1])

""" SUF plot """
x_ = np.linspace(0, -100, 100)    
fig, (ax) = plt.subplots(len(kx_), len(suf_[0]))
for i, kx in enumerate(kx_):
    ax[i, 0].set_ylabel("kx = {:g} $[cm^3 / day]$".format(kx))
    for j, suf in enumerate(suf_[i]):
        cols = ["r", "g", "b", "m"]
        for k in range(0, 4):
            ana = pb.SegmentAnalyser(r.rs) 
            ana.addData("SUF", suf)    
            ana.filter("subType", k + 1)
            y_ = ana.distribution("SUF", 0, -100, 100, True)
            ax[i, j].plot(y_, x_, cols[k])
        if i == 0:
            ax[i, j].set_title(str(20 + 10 * j))
            ax[i, j].set_xlabel("SUF [-]")
            ax[i, j].legend([str(l) for l in range(1, 5)])
plt.show()

