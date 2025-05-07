import sys; sys.path.append("../.."); sys.path.append("../../src/")

from vtk.util import numpy_support
import numpy as np

import plantbox as pb
import visualisation.vtk_plot as vp

from vtk_conversions import *


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
        "First shoot borne root [day]": (1., 30.),
        "Shoot borne delay [day]": (1., 21.),
        "First basal root [day]": (1., 30.),
        "Basal root delay [day]": (1., 21.),
        "Maximal number of basal roots [1]": (0., 50.),
        "First tiller [day]": (1., 30.),
        "Tiller delay [day]": (1., 21),
        "Maximal number of tillers [1]": (0., 25),
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


def simulate_plant(plant_, time_slider_value, seed_data, root_data, stem_data, leaf_data):
    """ simulates the plant xml parameter set with slider values """
    # 1. open base xml
    fname = get_parameter_names()[int(plant_)][1]
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    srp = plant.getOrganRandomParameter(pb.seed)
    rrp = plant.getOrganRandomParameter(pb.root)
    strp = plant.getOrganRandomParameter(pb.stem)
    lrp = plant.getOrganRandomParameter(pb.leaf)
    # 2. apply sliders to params
    apply_sliders(srp[0], seed_data["seed"], rrp, root_data, strp, stem_data, lrp, leaf_data)
    # 3. simulate
    N = 50
    t_ = np.linspace(0., time_slider_value, N + 1)
    v_, v1_, v2_, v3_ = np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    plant.initialize()
    for i, dt in enumerate(np.diff(t_)):
        plant.simulate(dt)
        ot = np.array(plant.getParameter("organType"))
        st = np.array(plant.getParameter("subType"))
        l = np.array(plant.getParameter("length"))
        v_[i] = np.sum(l)
        v1_[i] = np.sum(l[st == 1])
        v2_[i] = np.sum(l[st == 2])
        v3_[i] = np.sum(l[st == 3])

    # 4. make results store compatible (store pd stuff need for vtk.js, inlcuding different colours & 1D plots)
    pd = vp.segs_to_polydata(plant, 1., ["subType", "organType", "creationTime"])  # poly-data, "radius",
    tube = apply_tube_filter(pd)  # polydata + tube filter
    vtk_data = vtk_polydata_to_dashvtk_dict(tube)
    cellData = pd.GetCellData()
    cT = numpy_support.vtk_to_numpy(cellData.GetArray("creationTime"))
    vtk_data["creationTime"] = cT
    vtk_data["age"] = np.ones(cT.shape) * time_slider_value - cT
    vtk_data["subType"] = numpy_support.vtk_to_numpy(cellData.GetArray("subType"))
    vtk_data["organType"] = numpy_support.vtk_to_numpy(cellData.GetArray("organType"))
    vtk_data["time"] = t_[1:]

    # leaf geometry
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    leafes = plant.getOrgans(pb.leaf)
    for l in leafes:
        vp.create_leaf_(l, leaf_points, leaf_polys)
    pts_array = vtk.util.numpy_support.vtk_to_numpy(leaf_points.GetData()).astype(np.float32)
    polys_data = vtk.util.numpy_support.vtk_to_numpy(leaf_polys.GetData())
    vtk_data["leaf_points"] = pts_array.flatten().tolist()
    vtk_data["leaf_polys"] = polys_data.tolist()
    return vtk_data


def apply_sliders(srp, seed_data, rrp, root_data, strp, stem_data, lrp, leaf_data):
    tropism_names = { "Plagiotropism": 0, "Gravitropism":1, "Exotropism": 2 }
    # seed
    srp.firstSB = seed_data[0]
    srp.delaySB = seed_data[1]
    srp.firstB = seed_data[2]
    srp.delayB = seed_data[3]
    srp.maxB = seed_data[4]
    srp.firstTil = seed_data[5]
    srp.delayTil = seed_data[6]
    srp.maxTil = seed_data[7]
    # root
    print(root_data)
    for i, p in enumerate(rrp[1:]):
        print(p.name, p.subType)
        d = root_data[f"tab-{i+1}"]
        p.lmax = d[0]
        p.r = d[1]
        p.theta = d[2] / 180.*np.pi
        p.lb = d[3]
        p.ln = d[4]
        p.la = d[5]
        p.a = d[6]
        p.tropismN = d[7]
        p.tropismS = d[8]
        p.tropismT = tropism_names[d[9]]


def set_data(plant_, seed_data, root_data, root_typenames):
    print("set_data()")
    """ set root, seed, stem, and leaf data from xml """
    tropisms_names_ = { 0: "Plagiotropism", 1: "Gravitropism", 2: "Exotropism" }
    """ open xml """
    fname = get_parameter_names()[int(plant_)][1]  # xml filename
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    """ seed """
    srp = plant.getOrganRandomParameter(pb.seed)
    p = srp[0]
    seed_data["seed"] = [
        p.firstSB, p.delaySB,
        p.firstB, p.delayB, p.maxB,
        p.firstTil, p.delayTil, p.maxTil
    ]
    """ root """
    root_typenames.clear()
    rrp = plant.getOrganRandomParameter(pb.root)
    for i, p in enumerate(rrp[1:]):
        tropism_name = tropisms_names_[int(p.tropismT)]
        root_data[f"tab-{i+1}"] = [
            p.lmax, p.r, p.theta / np.pi * 180, p.lb, p.ln, p.la, p.a,
            p.tropismN, p.tropismS, tropism_name
        ]
        root_typenames[f"tab-{i+1}"] = p.name
    """ stem """

    """ leaf """


def param_to_dict(orp):

    print(orp.ogranType)
    if organType == pb.root:
        pass

    elif organType == pb.root:
        pass


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
    # print(rrp[1])

