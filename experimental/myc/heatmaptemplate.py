import sys; sys.path.append("../.."); sys.path.append("../../src/")

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.visualisation import figure_style
import numpy as np
import matplotlib.pyplot as plt
import time

nRings = 20
radius = radius/2
small_dish = pb.SDF_PlantContainer(radius*np.sqrt(1/nRings),radius*np.sqrt(1/nRings),height,False)
rings = []
rings.append(pb.SDF_Difference(small_dish, moved_helper_dish_hyphae))
for i in range(1, nRings):
    small_dish = pb.SDF_PlantContainer(radius*np.sqrt(i/nRings),radius*np.sqrt(i/nRings),height,False)
    small_hyphae_dish = pb.SDF_Difference(small_dish, moved_helper_dish_hyphae)
    old_dish = pb.SDF_Difference(pb.SDF_PlantContainer(radius*np.sqrt((i-1)/nRings),radius*np.sqrt((i-1) /nRings),height,False),moved_helper_dish_hyphae)
    small_hyphae_dish = pb.SDF_Difference(small_hyphae_dish,old_dish)
    rings.append(small_hyphae_dish)

times = np.linspace(0, max(mycp.getParameter("creationTime"))+0.1, 100)

def getParaDistperRing(parameter, times, plant, rings):
    paradenmat = np.zeros((len(rings),len(times[1:])))
    for k, ring in enumerate(rings):
        ringana = plant # need to copy the whole plant for segment analyzer since cripping to one ring removes all information outside
        ringana.crop(ring)
        for j in range(len(times[1:])):
            ringana.filter("creationTime", 0, np.flip(np.asarray(times))[j])
            ringana.pack()
            distrib = ringana.getSummed(parameter)
            paradenmat[k, len(times[1:])-1 -j] = np.array(distrib).sum() 
    return paradenmat

hld_matrix = getParaDistperRing("nodeTips", times, ana, rings)

print("shape:", hld_matrix.shape)
print("min:", np.min(hld_matrix))
print("max:", np.max(hld_matrix))
print("mean:", np.mean(hld_matrix))
plt.figure(figsize=(8, 5))

extent = [
    times[1], times[-1],   # x: Zeit
    1, len(rings)          # y: Ringe
]

plt.imshow(
    hld_matrix,
    aspect='auto',
    origin='lower',
    extent=extent,
    cmap='viridis'
)

plt.colorbar(label="Hyphal Length Density")

plt.xlabel("Time")
plt.ylabel("Ring (center → outside)")
plt.title("Radial movement of hyphal length density")

plt.show()