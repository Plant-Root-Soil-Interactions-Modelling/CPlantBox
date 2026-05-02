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
simtime = 4  # 95  # between 75-100 days
min_b = np.array([-5, -5, -20])
max_b = np.array([5, 5, 0])
cell_number = np.array([1, 1, 1])
width = max_b - min_b

s, soil = scenario.create_soil_model(soil_, min_b, max_b, cell_number, type=1)
r_ = scenario.create_mapped_rootsystem(min_b, max_b, cell_number, s, xml_name)  # pass parameter file for dynamic growth
# r_ = scenario.create_mapped_singleroot(min_b , max_b , cell_number, s)

r = r_.ms  # throw away (TODO ahve to change setup anyway...)

r.initialize(False)
r.simulate(simtime, False)

vp.plot_roots(r, "subType")  # , "outer_r") # OK

r.write("rootsystem_geometry.py")
ana = pb.SegmentAnalyser(r.mappedSegments())  # <------------DEAD
vp.plot_roots(ana, "subType")  # , "outer_r")
dd

ana.mapPeriodic(width[0], width[1])
ana.write("rootsystem.vtp")
vp.plot_roots(ana, "subType")  # , "outer_r")


sn = np.prod(cell_number)
peri = PerirhizalPython(r)

outer_radii = peri.get_outer_radii_voronoi()
# outer_radii = peri.get_outer_radii("surface")

print()
print("number of outer", outer_radii.shape)
print("open regions", np.sum(outer_radii == -1), np.sum(outer_radii == -1) / outer_radii.shape[0])
print("cell has point outside domain", np.sum(outer_radii == 0), np.sum(outer_radii == 0) / outer_radii.shape[0])
print()

print(
    "outer_radii:",
    "min",
    np.nanmin(outer_radii),
    "max",
    np.nanmax(outer_radii),
    "median",
    np.nanmedian(outer_radii),
    "mean",
    np.nanmean(outer_radii),
    "std",
    np.nanstd(outer_radii),
)

""" outer_r 3d visualisation """
ana = pb.SegmentAnalyser(r.mappedSegments())
ana.addData("outer_r", outer_radii)
grid = vp.uniform_grid(min_b, max_b, cell_number)
vp.plot_roots(ana, "subType")  # , "outer_r")

dd

""" nice histogram """
fig, axes = plt.subplots(1, 1, figsize=(8, 8))
outer_radii = peri.to_range_(outer_radii, 0.0, 2.0)
print(
    "outer_radii:",
    "min",
    np.min(outer_radii),
    "max",
    np.max(outer_radii),
    "median",
    np.median(outer_radii),
    "mean",
    np.mean(outer_radii),
    "std",
    np.std(outer_radii),
)
axes.hist(outer_radii, bins=40, rwidth=0.9)
plt.show()

# Write the VTK unstructured grid to a ParaView VTU file
domain = pb.SDF_Cuboid(pb.Vector3d(min_b), pb.Vector3d(max_b))
grid = peri.get_voronoi_mesh(domain)
writer = vtk.vtkXMLUnstructuredGridWriter()
writer.SetFileName("rootsystem_mesh.vtu")
writer.SetInputData(grid)
writer.Write()

print("fin")
