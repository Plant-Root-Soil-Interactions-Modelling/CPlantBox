"""analysis of results using signed distance functions"""
import py_rootbox as rb
import numpy as np
import matplotlib.pyplot as plt

rs = rb.RootSystem()
name = "Brassica_oleracea_Vansteenkiste_2014"
rs.readParameters("modelparameter/" + name + ".xml")
rs.initialize()
rs.simulate(120)

# Soil core analysis
r, depth, layers = 10, 100., 100
soilcolumn = rb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
soilcolumn2 = rb.SDF_RotateTranslate(soilcolumn, 0, 0, rb.Vector3d(10, 0, 0))  # shift 10 cm

# pick one geometry for further analysis
geom = soilcolumn

z_ = np.linspace(0, -1 * depth, layers)
fig, axes = plt.subplots(nrows = 1, ncols = 4, figsize = (16, 8))
for a in axes:
    a.set_xlabel('RLD (cm/cm)')
    a.set_ylabel('Depth (cm)')

# Make a root length distribution
ana = rb.SegmentAnalyser(rs)
rl_ = ana.distribution("length", 0., depth, layers, True)
axes[0].set_title('All roots (120 days)')
axes[0].plot(rl_, z_)

# Make a root length distribution along the soil core
ana = rb.SegmentAnalyser(rs)
# ana.crop(geom)
ana.pack()
rl_ = ana.distribution("length", 0., depth, layers, True)
axes[1].set_title('Soil core (120 days)')
axes[1].plot(rl_, z_)

# How it looked after 30 days?
ana = rb.SegmentAnalyser(rs)
ana.filter("creationTime", 0, 30)
# ana.crop(geom)
ana.pack()
rl_ = ana.distribution("length", 0., depth, layers, True)
axes[2].set_title('Soil core (30 days)')
axes[2].plot(rl_, z_)

# Only laterals?
ana = rb.SegmentAnalyser(rs)
ana.filter("subType", 2)  # assuming laterals are of type 2
ana.crop(geom)
ana.pack()
rl_ = ana.distribution("length", 0., depth, layers, True)
axes[3].set_title('Soil core, lateral roots (120 days)')
axes[3].plot(rl_, z_)

fig.subplots_adjust()
plt.savefig("../results/example_3b.png")
plt.show()

print("done.")
