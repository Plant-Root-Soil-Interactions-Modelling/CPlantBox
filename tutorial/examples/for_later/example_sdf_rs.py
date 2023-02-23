import py_rootbox as rb

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.colors as colors

rs = rb.RootSystem()

# Open plant and root parameter from a file
name = "Anagallis_femina_Leitner_2010"
rs.openFile(name)

# Create and set geometry

# 1.creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50 cm, not square but circular
soilcore = rb.SDF_PlantContainer(5, 5, 40, False)

# 2. creates a square 27*27 cm containter with height 1.4 cm
rhizotron = rb.SDF_PlantBox(0.5, 27, 27)

# Pick 1, or 2
rs.setGeometry(rhizotron)  # soilcore, or rhizotron

# Initialize
rs.initialize()

# Simulate
rs.simulate(60)  # days

# Root system distance function
ana = rb.SegmentAnalyser(rs)
radii = ana.getParameter("radius")

print(radii[10])
rs_sdf = rb.SDF_RootSystem(ana.nodes, ana.segments, radii, 1)

# Distance Grid
nx = 1
ny = 1000
nz = 1000
X = np.linspace(-0.25 / 2, 0.25 / 2, nx)
Y = np.linspace(-27 / 2, 27 / 2, ny)
Z = np.linspace(0, -27, nz)
X_, Y_, Z_ = np.meshgrid(X, Y, Z, indexing = "ij")  # stupid matlab default
D = np.zeros(X_.shape)
print(D.shape)

for i in range(0, X_.shape[0]):
    for j in range(0, X_.shape[1]):
        for k in range(0, X_.shape[2]):
            D[i, j, k] = rs_sdf.getDist(rb.Vector3d(X_[i, j, k], Y_[i, j, k], Z_[i, j, k]))

D[D < -100] = -10

fig1 = plt.figure()
ax = plt.axes()

D_ = D[int(nx / 2), :, :]
levels = np.linspace(-5, 0.05, 100)
cs = ax.contourf(Y_[int(nx / 2), :, :], Z_[int(nx / 2), :, :], D_, levels = levels, cmap = 'jet')  # levels = levels, locator = ticker.LogLocator(),
ax.set_xlabel('x')
ax.set_ylabel('z')
plt.axis('equal')
cbar = fig1.colorbar(cs)
plt.show()

print("max ", max(D_.flatten()))

# Export final result (as vtp)
rs.write("../results/example_1b.vtp")

# Export container geometry as Paraview Python script
rs.write("../results/example_1b.py")
