"""virtual sampling of root length densities over time"""

import matplotlib.pyplot as plt
import numpy as np

import plantbox as pb
from plantbox.visualisation import figure_style

path = "../../modelparameter/structural/rootsystem/"
filename = "wheat"

months = 8  # |\label{l2_3:timebegin}|
times = np.linspace(0, 30 * months, months + 1)  # |\label{l2_3:timeend}|

# 72 cm*45 cm size plot
n_plants = 16  # number of plants in rows #|\label{l2_3:plantsetbegin}|
n_colums = 7  # number of rows
distance_plants = 3  # distance between the root systems along row (cm)
distance_rows = 12  # distance between the rows (cm)
interrow = (n_plants - 1) * distance_plants  # intra-row spacing (cm)
row = (n_colums - 1) * distance_rows  # row spacing (cm) |\label{l2_3:plantsetend}|

r, depth, layers = 4.2 / 2, 160.0, 32  # Soil core analysis #|\label{l2_3:soilcorebegin}|
layerVolume = depth / layers * r * r * np.pi  # cm3
z_ = np.linspace(0, -depth, layers)  # slices of soil core

soilcolumn = pb.SDF_PlantContainer(r, r, depth, False)  # square = False

soilcor_x = interrow / 4
soilcor_y = distance_rows
x_ = [11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75, 11.25, 22.5, 33.75]
y_ = [66.0, 66.0, 66.0, 54.0, 54.0, 54.0, 42.0, 42.0, 42.0, 30.0, 30.0, 30.0, 18.0, 18.0, 18.0]
soilcolumns_ = [pb.Vector3d(x_ij, y_ij, 0) for x_ij, y_ij in zip(x_, y_)]
soilcolumns = [pb.SDF_RotateTranslate(soilcolumn, vi) for vi in soilcolumns_]  # |\label{l2_3:soilcoreend}|

soilSpace = pb.SDF_PlantContainer(500, 500, 500, True)

for i in range(0, n_plants):  # |\label{l2_3:simulationbegin}|
    for j in range(0, n_colums):
        plant = pb.Plant()
        plant.readParameters(path + filename + ".xml", fromFile = True, verbose = False)
        seed = plant.getOrganRandomParameter(pb.seed)[0]
        seed.seedPos = pb.Vector3d(distance_plants * i, distance_rows * j, -3.0)  # cm
        plant.setGeometry(soilSpace)
        plant.initialize(False)
        plant.simulate(30 * months, False)
        if i + j == 0:
            all_ana = pb.SegmentAnalyser(plant)
        else:
            all_ana.addSegments(plant)  # |\label{l2_3:simulationend}|
rld = np.zeros([len(soilcolumns) * len(times[1:]), layers])

for k, sc in enumerate(soilcolumns):  # |\label{l2_3:soilcolselectbegin}|
    ana = pb.SegmentAnalyser(all_ana)  # copy all
    ana.crop(sc)  # select soil column
    for j in range(len(times[1:])):
        ana.filter("creationTime", 0, np.flip(np.asarray(times))[j])
        ana.pack()
        distrib = ana.distribution("length", 0.0, -depth, layers, True)
        rld[(len(times[1:]) - 1 - j) * len(soilcolumns) + k] = np.array(distrib) / layerVolume  # |\label{l2_3:soilcolselectend}|

legend_lst = [str(int(t_i)) for t_i in times[1:]]  # |\label{l2_3:plotbegin}|
fig, axes = plt.subplots(nrows = 5, ncols = int(len(soilcolumns) / 5), sharex = True, sharey = True, figsize = (8, 16))

for k in range(len(soilcolumns)):
    axes.flat[k].set_title("Soil core" + " " + str(k + 1))
    for j in range(len(times[1:])):
        axes.flat[k].plot(np.array(rld[len(soilcolumns) * j + k]), z_)
        axes.flat[k].set_xlim(0, 5)

plt.setp(axes[-1,:], xlabel = "RLD $(cm/cm^3)$")
plt.setp(axes[:, 0], ylabel = "Depth $(cm)$")
plt.legend(np.asarray(legend_lst), loc = "lower center", bbox_to_anchor = (-0.8, -0.5), ncol = 8)
fig.subplots_adjust()
plt.savefig("results/rld_plot.png", dpi = 300, bbox_inches = "tight")
plt.show()  # |\label{l2_3:plotend}|
