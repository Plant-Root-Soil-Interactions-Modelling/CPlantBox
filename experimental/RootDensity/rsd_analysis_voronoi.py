"""small example"""

import matplotlib.pyplot as plt
import numpy as np
import scenario_setup as scenario
import vtk

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp
from plantbox.functional.Perirhizal import *

""" parameters """
# soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.soybean(0)  # 0 = envirotype
# xml_name = "data/Glycine_max_Moraes2020_opt2_modified.xml"  # root growth model parameter file
# simtime = 87.5  # 87.5  # between 75-100 days
# # cell_number = [76, 4, 200]
# cell_number = [38, 2, 100]
# # cell_number = [19, 1, 50]

soil_, table_name, min_b, max_b, cell_number, area, Kc = scenario.maize(0)  # 0 = envirotype
xml_name = "data/Zeamays_synMRI_modified.xml"  # root growth model parameter file
simtime = 25  # 95  # between 75-100 days
# cell_number = [76, 16, 200]
# cell_number = [38, 8, 100]
cell_number = [19, 4, 50]

width = max_b - min_b

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type=1)
r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# r_ = scenario.create_mapped_singleroot(min_b , max_b , cell_number, s)

r = r_.ms  # throw away (TODO have to change setup anyway...)

r.initialize(False)
r.simulate(simtime, False)

r.write("rootsystem_geometry_large.py")
ana = pb.SegmentAnalyser(r.mappedSegments())
ana.mapPeriodic(width[0], width[1])
ana.write("rootsystem_large.vtp")

sn = np.prod(cell_number)
peri = PerirhizalPython(r)

""" add 3d soil surface density """
sd = peri.get_density("surface")
sd = -np.minimum(sd, 1.1)  # limit for visualisation
grid = vp.uniform_grid(min_b, max_b, cell_number)
cell_sd = vtk.vtkDoubleArray()
cell_sd.SetName("surface_density")
cell_sd.SetNumberOfValues(sn)
for j in range(0, sn):
    cell_sd.SetValue(j, sd[j])
celldata = grid.GetCellData()
celldata.AddArray(cell_sd)

outer_radii = peri.get_outer_radii_voronoi()
# outer_radii = peri.get_outer_radii("surface")

print()
print("number of outer", outer_radii.shape)
print("open regions", np.sum(outer_radii == -1), np.sum(outer_radii == -1) / outer_radii.shape[0])
print("cell has point outside domain", np.sum(outer_radii == 0), np.sum(outer_radii == 0) / outer_radii.shape[0])
print()

print("outer_radii", np.min(outer_radii), np.max(outer_radii), "median", np.median(outer_radii), "mean", np.mean(outer_radii), np.std(outer_radii))
# vp.plot_mesh(grid, "surface_density")
# vp.plot_mesh_cuts(grid, "surface_density")
# vp.plot_roots(ana, "outer_r")

""" outer_r 3d visualisation """
ana = pb.SegmentAnalyser(r)
outer_radii = -np.minimum(outer_radii, 1.0)  # limit for visualisation
ana.addData("outer_r", outer_radii)
vp.plot_roots_and_mesh(ana, "subType", grid, "surface_density", True, width[0], width[1])  # outer_r

""" nice histogram """
fig, axes = plt.subplots(1, 1, figsize=(8, 8))
rr = peri.to_range_(outer_radii, 0.0, 2.0)
axes.hist(rr, bins=40, rwidth=0.9)
plt.show()

# Write the VTK unstructured grid to a ParaView VTU file
domain = pb.SDF_Cuboid(pb.Vector3d(min_b), pb.Vector3d(max_b))
grid = peri.get_voronoi_mesh(domain)
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("rootsystem_mesh_large.vtu")
writer.SetInputData(grid)
writer.Write()

print("fin")
