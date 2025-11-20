"""root system length over time"""
import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import visualisation.vtk_plot as vp
import rsml.rsml_reader as rsml

import numpy as np
import matplotlib.pyplot as plt

path = ""
name = "WurzelAtlasParasTomate"
name_new = "WurzelAtlasTomate_parameters_out"
name_Johanna_WildType = "TomatoJohanna_WildType"
name_Johanna_RMC = "TomatoJohanna_RMC"

animation = False
pot = True
radius = 11.7
depth = 110


tomato_WAT = pb.MycorrhizalPlant()
tomato_J_WT = pb.MycorrhizalPlant()
tomato_J_RMC = pb.MycorrhizalPlant()

tomato_WAT.readParameters(path +name_new + ".xml", fromFile = True, verbose = True)
tomato_J_WT.readParameters(path +name_Johanna_WildType + ".xml", fromFile = True, verbose = True)
tomato_J_RMC.readParameters(path +name_Johanna_RMC + ".xml", fromFile = True, verbose = True)

if pot:
    # Pot dimensions 11.7 cm x 11.7 cm x 30 cm
    depth = 30
    pot = pb.SDF_PlantContainer(radius, radius, depth, False)
    tomato_WAT.setGeometry(pot)
    tomato_J_WT.setGeometry(pot)
    tomato_J_RMC.setGeometry(pot) 



tomato_WAT.initialize(True)
tomato_J_WT.initialize(True)
tomato_J_RMC.initialize(True)

### Set up simulation parameters
totalduration = 108
potduration = 73
pottime = totalduration-potduration
simtime = totalduration
fpd = 1 # one step per day
N = simtime * fpd
dt = simtime / N

for i in range(1,N+1):
    if i % 18 == 0:
        print('step',i, '/',N)
    tomato_WAT.simulate(dt, False)
    tomato_J_WT.simulate(dt, False)
    tomato_J_RMC.simulate(dt, False)
print('done')

def make_RLD_plot(rs, depth, layers,radius,name):
    z_ = np.linspace(0, -1 * depth, layers)
    #fig, axes = plt.subplots(nrows = 1, ncols = 1, figsize = (16, 8))
    fig, axes = plt.subplots(figsize = (6, 8))
    # Make a root length distribution
    ana = pb.SegmentAnalyser(rs)
    rad = ana.getParameter("radius")
    rad = sum(rad)/len(rad)
    rltot = ana.getSummed("length")
    print("Simulated average root radius and total root length after", simtime, "days:", "\n","  Average root radius (cm):", rad, "\n","  Total root length (cm):", rltot, "\n" )
    layerVolume = depth / layers * radius * radius * np.pi  # actually the only thing that changes
    ana.write("results/" + name + "_108days.vtp",["radius", "subType", "creationTime","organType"])
    rl0_ = ana.distribution("length", 0., -depth, layers, True)
    ana.filter("creationTime", 0, 80)
    rl1_ = ana.distribution("length", 0., -depth, layers, True)
    ana.filter("creationTime", 0, 50)
    rl2_ = ana.distribution("length", 0., -depth, layers, True)
    ana.filter("creationTime", 0, 35)
    rl3_ = ana.distribution("length", 0., -depth, layers, True)
    axes.set_title('All roots')
    axes.plot(np.array(rl3_) / layerVolume * 10, z_)
    axes.plot(np.array(rl2_) / layerVolume * 10, z_)
    axes.plot(np.array(rl1_) / layerVolume * 10, z_)
    axes.plot(np.array(rl0_) / layerVolume * 10, z_, color = 'goldenrod')
    axes.set_xlabel('$\dfrac{1}{10}$ RLD (cm/cm^3)')
    axes.set_ylabel('Depth (cm)')
    axes.legend(["35 days", "50 days", "80 days", "108 days"], loc = 'upper right')
    #axes.set_xlim(0,2)
    #axes.set_ylim(-31,0)
    fig.subplots_adjust()
    plt.savefig("results/" + name + "10th_RLD.pdf", dpi = 300)

    fig, ax = plt.subplots()
    ax.plot(np.array(rl3_) / layerVolume, z_)
    ax.plot(np.array(rl2_) / layerVolume, z_)
    ax.plot(np.array(rl1_) / layerVolume, z_)
    ax.plot(np.array(rl0_) / layerVolume, z_)
    ax.set_xlabel('RLD (cm/cm^3)')
    ax.set_ylabel('Depth (cm)')
    ax.legend(["35 days", "50 days", "80 days", "108 days"], loc = 'upper right')
    # ax.set_ylim(-31,0)
    fig.subplots_adjust()
    plt.savefig("results/" + name + "_RLD.pdf", dpi = 300)

    # plt.show()
make_RLD_plot(tomato_WAT, 30, 60, radius, name_new)
make_RLD_plot(tomato_J_WT, 110, 220, radius, name_Johanna_WildType)
make_RLD_plot(tomato_J_RMC, 110, 220, radius, name_Johanna_RMC)

# ana = pb.SegmentAnalyser(mycp)
# rad = ana.getParameter("radius")
# rad = sum(rad)/len(rad)
# rl = ana.getSummed("length")
# print("\n")
# print("Simulated Root Atlas Type:")
# print("Root length total (cm):", rl)
# # print("Difference to expected lengths (cm):", rl - WT[2])
# # print("Standard deviation of expected lengths (cm):", WT[3])
# print("Average root radius (cm):", rad)
# # print("Difference to expected radius (cm):", rad - WT[0])
# # print("Standard deviation of expected radius (cm):", WT[1])

# vp.plot_roots(mycp)
# ana.write(name + ".vtp", ["radius", "subType", "creationTime","organType"])
if pot:
    vp.write_container(pot, "results/pot.vtp")

