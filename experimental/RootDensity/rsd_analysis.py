"""3d surface densities"""

import matplotlib.pyplot as plt
import numpy as np
import scenario_setup as scenario
import vtk

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.Perirhizal import *

""" parameters """
soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)  # 0 = envirotype
xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
simtime = 25.5  # between 75-100 days
# cell_number = [76, 4, 200]
cell_number = [38, 2, 100]
# cell_number = [19, 1, 50]
# cell_number = [1, 1, 1]

# soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)  # 0 = envirotype
# xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
# simtime = 25  # between 75-100 days
# # cell_number = [76, 16, 200]
# # cell_number = [38, 8, 100]
# cell_number = [19, 4, 50]
# cell_number = [1, 1, 1]

width = max_b - min_b

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type=1)
r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth

r = r_.ms  # throw away (TODO ahve to change setup anyway...)

r.initialize(False)
r.simulate(simtime, True)

sn = np.prod(cell_number)
peri = PerirhizalPython(r)

""" add 3d soil surface density """
sd = peri.get_density("surface")
sd = np.minimum(sd, 1.1)  # limit for visualisation
grid = vp.uniform_grid(min_b, max_b, cell_number)
cell_sd = vtk.vtkDoubleArray()
cell_sd.SetName("surface_density")
cell_sd.SetNumberOfValues(sn)
for j in range(0, sn):
    cell_sd.SetValue(j, sd[j])
celldata = grid.GetCellData()
celldata.AddArray(cell_sd)

outer_radii = peri.get_outer_radii("voronoi")
# outer_radii = peri.get_outer_radii("length")

print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))

ana = pb.SegmentAnalyser(r)
outer_radii = np.minimum(outer_radii, 3.0)  # limit for visualisation
ana.addData("outer_r", outer_radii)
vp.plot_roots(ana, "outer_r")

ana.addData("radius", outer_radii)
vp.plot_mesh(grid, "surface_density")
vp.plot_mesh_cuts(grid, "surface_density")
vp.plot_roots(ana, "outer_r")
vp.plot_roots_and_mesh(ana, "outer_r", grid, "surface_density", True, width[0], width[1])

plt.hist(outer_radii, bins=100, rwidth=0.9, align="mid")
plt.show()

print("fin")
