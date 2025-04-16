import sys; sys.path.append("../.."); sys.path.append("../../src/")

import vtk
import numpy as np

import plantbox as pb
import visualisation.vtk_plot as vp


def get_parameter_names():
    """ returns a list of plant parameter names with two values each, 
      first a short name, second exact filename """
    parameter_names = [
        ("Maize", "Zea_mays_1_Leitner_2010.xml"),
        ("Anagallis", "Anagallis_femina_Leitner_2010.xml"),
        ("FSPM", "fspm2023.xml") ]
    return parameter_names


def get_seed_slider_names():
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "First shoot borne root [day]": (1., 180.),
        "Shoot borne delay [day]": (1., 21.),
        "First basal root [day]": (1., 180),
        "Basal root delay [day]": (1., 180),
        "Maximal number of basal roots [1]": (0., 100.),
        "First tiller [day]": (1., 180),
        "Tiller delay [day]": (1., 21),
        "Maximal number of tillers [1]": (0., 30),
    }
    return parameter_sliders


def get_root_slider_names():
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 200),
        "Initial growth rate [cm/day]": (0.5, 10),
        "Initial angle [°]": (0., 90.),
        "Basal zone [cm]": (0.1, 20),
        "Interlateral distance [cm]": (0.1, 20),
        "Apical zone [cm]": (0.1, 20),
        "Radius [cm]": (1.e-3, 0.25),
        "Tropism strength [1]": (0, 6),
        "Tropism tortuosity [1]": (0., 1.),
    }
    return parameter_sliders


def vtk_polyline_to_dict(polydata):
    """ converts polylines to a dict """
    points = polydata.GetPoints()
    lines = polydata.GetLines()

    n_points = points.GetNumberOfPoints()
    n_lines = lines.GetNumberOfCells()

    # Points array
    pts = np.array([points.GetPoint(i) for i in range(n_points)], dtype = np.float32)
    pts = pts.flatten().tolist()

    # Lines connectivity array
    lines.InitTraversal()
    id_list = vtk.vtkIdList()
    conn = []

    for _ in range(n_lines):
        lines.GetNextCell(id_list)
        conn.append(id_list.GetNumberOfIds())
        for j in range(id_list.GetNumberOfIds()):
            conn.append(id_list.GetId(j))

    return {
        "points": pts,
        "lines": conn,
    }


def simulate_plant(plant_, time_slider_value):
    """ """
    # 1. open xml
    fname = get_parameter_names()[int(plant_)][1]
    # apply slider values
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    plant.initialize()
    plant.simulate(time_slider_value)
    return plant


def set_sliders_from_xml(plant_, data):
    """ for roots (todo) """
    fname = get_parameter_names()[int(plant_)][1]  # xml filename
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    rrp = plant.getOrganRandomParameter(pb.root)
    for i, p in enumerate(rrp[1:]):
        data[f"tab-{i}"] = [p.lmax, p.r, p.lb, p.ln, p.la, p.a, p.tropismN, p.sigma ]


def debug_params(plant_):
    fname = get_parameter_names()[int(plant_)][1]  # xml filename
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    seedrp = plant.getOrganRandomParameter(pb.seed)
    rrp = plant.getOrganRandomParameter(pb.root)
    stemrp = plant.getOrganRandomParameter(pb.stem)
    lrp = plant.getOrganRandomParameter(pb.leaf)
    print()
    print(fname)
    print("Seed", len(seedrp), seedrp[0].name)
    print("Root", len(rrp), rrp[1].name)
    print("Stem", len(stemrp), stemrp[1].name)
    print("Leaf", len(lrp), lrp[1].name)
    print()
    print(rrp[1])

