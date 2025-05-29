import sys; sys.path.append("../.."); sys.path.append("../../src/")

from vtk.util import numpy_support
import numpy as np

import plantbox as pb
import visualisation.vtk_plot as vp

from vtk_conversions import *

tropism_names = { "Plagiotropism": 0, "Gravitropism":1, "Exotropism": 2, "Negative gravitropism": 4, "Variable gravitropism": 6}  # "Twist": 5,
tropism_names_ = { 0: "Plagiotropism", 1: "Gravitropism", 2: "Exotropism", 4: "Negative gravitropism", 6: "Variable gravitropism"}  # 5: "Twist",


def get_parameter_names():  # parameter xml file names
    """ returns a list of plant parameter names with two values each, first a short name, second exact filename """
    parameter_names = [
        ("Demo Leaf", "leaf_only.xml"),
        ("Demo Root", "root_only.xml"),
        ("Demo Stem", "stem_only.xml"),
#        ("Maize3", "P3.xml"),
        ("Maize", "P0.xml"),
        ("Wheat", "Triticum_aestivum_test_2021.xml"),  # Monas File
        ("FSPM", "fspm2023.xml") ]
    return parameter_names

# Felix Maximilian Bauer, Dirk Norbert Baker, Mona Giraud, Juan Carlos Baca Cabrera, Jan Vanderborght, Guillaume Lobet, Andrea Schnepf, Root system architecture reorganization under decreasing soil phosphorus lowers root system conductance of Zea mays, Annals of Botany, 2024;, mcae198, https://doi.org/10.1093/aob/mcae198


def get_seed_slider_names():  # see set_data, apply_sliders
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "First shoot borne root [day]": (1, 30),
        "Shoot borne delay [day]": (1, 21),
        "First basal root [day]": (1, 21),
        "Basal root delay [day]": (1, 21),
        "Maximal number of basal roots [1]": (0, 30),
        "First tiller [day]": (1, 21),
        "Tiller delay [day]": (1, 21),
        "Maximal number of tillers [1]": (0, 7),
    }
    return parameter_sliders


def get_root_slider_names():  # see set_data, apply_sliders
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 200),
        "Growth rate [cm/day]": (0.5, 10),
        "Initial angle [°]": (0., 90),
        "Basal zone [cm]": (0.1, 20),
        "Interlateral distance [cm]": (0.1, 20),
        "Apical zone [cm]": (0.1, 20),
        "Radius [cm]": (1.e-3, 0.25),
        "Tropism strength [1]": (0, 6),
        "Tropism tortuosity [1]": (0, 1),
    }
    return parameter_sliders


def get_stem_slider_names():  # see set_data, apply_sliders
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 200),
        "Growth rate [cm/day]": (0.5, 10),
        "Initial angle [°]": (0, 90),
        "Phytomer distance [cm]": (0.1, 20),
        "Radius [cm]": (1.e-3, 0.25),
        "Nodal growth start [day]": (0, 21),
        "Nodal growth time span [day]": (0, 21),
        "Fixed Rotation [°]": (0, 180),
        "Tropism strength [1]": (0, 6),
        "Tropism tortuosity [1]": (0, 1),
    }
    return parameter_sliders


def get_leaf_slider_names():  # see set_data, apply_sliders
    """ return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 50),
        "Growth rate [cm/day]": (0.5, 10),
        "Initial angle [°]": (0, 180),
        "Petiole length [cm]": (0.1, 10),  # lb
        "Fixed Rotation [°]": (0, 180),
        "Tropism strength [1]": (0, 6),
        "Tropism tortuosity [1]": (0, 1),
    }
    return parameter_sliders


def set_leaf_geometry(shapename, p):
    """ shape name "Defined", "Long", "Round", "Maple", "Flower" """
    if shapename == "Defined":
        return
    elif shapename == "Long":
        p.lb = 1  # length of leaf stem
        p.la, p.lmax = 3.5, 8.5
        p.areaMax = 10  # cm2, area reached when length = lmax
        phi = np.array([-90, -45, 0., 45, 90]) / 180. * np.pi
        l = np.array([3, 2.2, 1.7, 2, 3.5])  # distance from leaf center
    elif shapename == "Round":
        p.lb = 1  # length of leaf stem
        p.la, p.lmax = 5, 11
        p.areaMax = 3.145 * (p.la ** 2)
        phi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
        l_ = (p.lmax - p.lb) / 2  # == p.la
        l = np.array([l_ for x_i in range(len(phi))])  # ([2, 2, 2, 4,1,1, 4]) #distance from leaf center
    elif shapename == "Maple":
        p.lb = 1  # length of leaf stem
        p.areaMax = 50
        p.la, p.lmax = 5, 11
        phi = np.array([-90, -45, 0., 45, 67.5, 70, 90]) / 180. * np.pi
        l = np.array([5, 2, 2, 4, 1, 1, 5])  # distance from leaf center
    elif shapename == "Flower":
        p.lb = 1  # length of leaf stem
        p.areaMax = 50
        p.la, p.lb, p.lmax = 5, 1, 11
        phi = np.array([-90., -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180. * np.pi
        l = np.array([5., 1, 4.5, 1, 4, 1, 4.5, 1, 5])    
    p.createLeafRadialGeometry(phi, l, 100)


def set_leaf_sliders(slider_values):
    """ TODO cannot get the callbacks running in dash"""
    return slider_values


def fix_dx(rrp, strp, lrp):
    """ overrides the xml resolution settings to ensure smooth visualization """
    for r in rrp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.5
        r.dxMin = 1.e-6
    for r in strp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.5
        r.dxMin = 1.e-6
        r.nodalGrowth = 1  # <------- !
        # r.initBeta = 0.
        # r.betaDev = 0.
        print("delayNGStarts", r.delayNGStarts)
        print("delayNGEnds", r.delayNGEnds)
        #print("initBeta", r.initBeta)
        #print("betaDev", r.betaDev)
        #print("rotBeta", r.rotBeta)
        # r.rotBeta = 0.5
        #r.betaDev = 10
    for r in lrp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.25
        r.dxMin = 1.e-6
        r.initBeta = 0.
        r.betaDev = 0.
        r.rotBeta = 0.5


def simulate_plant(plant_, time_slider, seed_data, root_data, stem_data, leaf_data):
    """ simulates the plant xml parameter set with slider values """
    print("simulate_plant()")
    # 1. open base xml
    fname = get_parameter_names()[int(plant_)][1]
    plant = pb.Plant()
    plant.readParameters("params/" + fname)
    srp = plant.getOrganRandomParameter(pb.seed)
    rrp = plant.getOrganRandomParameter(pb.root)
    strp = plant.getOrganRandomParameter(pb.stem)
    lrp = plant.getOrganRandomParameter(pb.leaf)
    fix_dx(rrp, strp, lrp)
    number_r = len(rrp[1:])  # number of root types
    number_s = len(strp[1:])  # number of stem types
    # 2. apply sliders to params
    apply_sliders(srp[0], seed_data, rrp, root_data, strp, stem_data, lrp, leaf_data)
    srp[0].seedPos.x = 0.  # override position (always)
    srp[0].seedPos.y = 0.
    srp[0].seedPos.z = -3.
    srp[0].delayRC = 30.
    srp[0].delayDefinitionShoot = 2
    print("delaySB", srp[0].delaySB)
    print("firstSB", srp[0].firstSB)
    print("delayRC", srp[0].delayRC)
    print("nC", srp[0].nC)
    # 3. simulate
    N = 25
    t_ = np.linspace(0., time_slider, N + 1)
    root_length = np.zeros((number_r, N))
    stem_length = np.zeros((number_s, N))
    leaf_length = np.zeros((N,))
    plant.initialize()
    rld_, z_ = [], []
    for i, dt in enumerate(np.diff(t_)):
        plant.simulate(dt)
        ot = np.array(plant.getParameter("organType"))
        st = np.array(plant.getParameter("subType"))
        l = np.array(plant.getParameter("length"))
        for j in range(number_r):
            root_length[j, i] = np.sum(l[np.logical_and(st == j + 1, ot == pb.root)])
        for j in range(number_s):
            stem_length[j, i] = np.sum(l[np.logical_and(st == j + 1, ot == pb.stem)])
        leaf_length[i] = np.sum(l[ot == pb.leaf])
        if i in [4, 9, 14, 19, 24]:  # 4. make depth profiles
            ana = pb.SegmentAnalyser(plant)
            length = ana.getSummed("length")
            max_ = ana.getMaxBounds().z
            min_ = ana.getMinBounds().z
            rld_.append(np.array(ana.distribution("length", max_, min_, int(np.round(max_ - min_)), True)))
            z_.append(np.linspace(max_, min_, int(np.round(max_ - min_))))

    # 5. make results store compatible (store pd stuff need for vtk.js, inlcuding different colours & 1D plots)
    pd = vp.segs_to_polydata(plant, 1., ["subType", "organType", "radius", "creationTime"])  # poly-data, "radius",
    tube = apply_tube_filter(pd)  # polydata + tube filter
    vtk_data = vtk_polydata_to_dashvtk_dict(tube)
    cellData = pd.GetCellData()
    cT = numpy_support.vtk_to_numpy(cellData.GetArray("creationTime"))
    vtk_data["creationTime"] = cT  ################################################################### TODO somehow creationTime and Age are mixed up
    for i in range(0, len(rld_)):
        vtk_data[f"rld{i}"] = rld_[i] / length
        vtk_data[f"z{i}"] = z_[i]

    # vtk_data["age"] = np.ones(cT.shape) * time_slider_value - cT
    organType = numpy_support.vtk_to_numpy(cellData.GetArray("organType"))  #
    vtk_data["subType"] = numpy_support.vtk_to_numpy(cellData.GetArray("subType")) + 5 * (organType - np.ones(organType.shape) * 2)
    vtk_data["radius"] = numpy_support.vtk_to_numpy(cellData.GetArray("radius"))
    vtk_data["time"] = t_[1:]
    for j in range(number_r):
        vtk_data[f"root_length-{j+1}"] = list(root_length[j,:])
    for j in range(number_s):
        vtk_data[f"stem_length-{j+1}"] = list(stem_length[j,:])
    vtk_data["leaf_length"] = list(leaf_length[:])
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
    # general
    vtk_data["number_r"] = number_r
    vtk_data["number_s"] = number_s

    return vtk_data


def apply_sliders(srp, seed_data, rrp, root_data, strp, stem_data, lrp, leaf_data):
    """ slider values to cplantbox random parameters """
    # seed
    s = seed_data["seed"]
    srp.firstSB = s[0]
    srp.delaySB = s[1]
    srp.firstB = s[2]
    srp.delayB = s[3]
    srp.maxB = s[4]
    srp.firstTil = s[5]
    srp.delayTil = s[6]
    srp.maxTil = int(s[7] + 0.5)
    if not seed_data["shoot-checkbox"]:
        # print("no shootborne")
        srp.firstSB = 1.e6  # disable
        srp.delaySB = 1.e6
    if not seed_data["basal-checkbox"]:
        srp.firstB = 1.e6  # disable
    if not seed_data["tillers-checkbox"]:
        srp.firstTil = 1.e6  # disable
    # root
    for i, p in enumerate(rrp[1:]):
        # print(p.name, p.subType)
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
    # stem
    for i, p in enumerate(strp[1:]):
        # print(p.name, p.subType)
        d = stem_data[f"tab-{i+1}"]
        p.lmax = d[0]
        p.r = d[1]
        p.theta = d[2] / 180.*np.pi
        p.ln = d[3]
        p.a = d[4]
        p.delayNGStart = d[5]
        p.delayNGEnd = d[5] + d[6]
        p.rotBeta = d[7]/ 180.
        p.tropismN = d[8]
        p.tropismS = d[9]
        p.tropismT = tropism_names[d[10]]
    # Leaf
    if len(lrp) > 1:
        p = lrp[1]
        d = leaf_data["leaf"]
        set_leaf_geometry(d[0], p)
        p.lmax = d[1]
        p.r = d[2]
        p.theta = d[3] / 180 * np.pi
        p.lb = d[4]
        p.rotBeta = d[5] / 180.
        p.tropismN = d[6]
        p.tropismS = d[7]
        p.tropismT = tropism_names[d[8]]


def set_data(plant_, seed_data, root_data, stem_data, leaf_data, typename_data):
    """ sets all store data from the xml file """
    # print("set_data()")
    typename_data.clear()
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
    seed_data["shoot-checkbox"] = p.firstSB < 1.e3 and p.delaySB < 1.e3
    seed_data["basal-checkbox"] = p.firstB < 1.e3 and p.delayB < 1.e3 and p.maxB > 0
    seed_data["tillers-checkbox"] = p.firstTil < 1.e3 and p.delayTil < 1.e3 and p.maxTil > 0
    if not seed_data["basal-checkbox"]:  # defaults in case someone turns it on
        seed_data["seed"][2] = 7
        seed_data["seed"][3] = 7
        seed_data["seed"][4] = 5
    if not seed_data["tillers-checkbox"]:  # defaults in case someone turns it on
        seed_data["seed"][5] = 7
        seed_data["seed"][6] = 11
        seed_data["seed"][7] = 4
    seed_data["simulationTime"] = p.simtime  # where to put it
    """ root """
    rrp = plant.getOrganRandomParameter(pb.root)
    typename_data["number_roottypes"] = len(rrp[1:])
    for i, p in enumerate(rrp[1:]):
        tropism_name = tropism_names_[int(p.tropismT)]
        root_data[f"tab-{i+1}"] = [
            p.lmax, p.r, p.theta / np.pi * 180, p.lb, p.ln, p.la, p.a,
            p.tropismN, p.tropismS, tropism_name , len(p.successorST) > 0
        ]
        typename_data[f"root tab-{i+1}"] = p.name
    """ stem """
    strp = plant.getOrganRandomParameter(pb.stem)
    typename_data["number_stemtypes"] = len(strp[1:])
    for i, p in enumerate(strp[1:]):
        tropism_name = tropism_names_[int(p.tropismT)]
        stem_data[f"tab-{i+1}"] = [
            p.lmax, p.r, p.theta / np.pi * 180, p.ln, p.a,
            p.delayNGStart, p.delayNGEnd - p.delayNGStart, p.rotBeta * 180,
            p.tropismN, p.tropismS, tropism_name, len(p.successorST) > 0  # for p.rotBeta = 1 == 180 ° (???)
        ]
        typename_data[f"stem tab-{i+1}"] = p.name
    """ leaf """
    lrp = plant.getOrganRandomParameter(pb.leaf)
    if len(lrp) > 1:
        p = lrp[1]
        tropism_name = tropism_names_[int(p.tropismT)]
        leaf_data["leaf"] = [
            "Defined",
            p.lmax, p.r,
            p.theta / np.pi * 180, p.lb, p.rotBeta * 180,
            p.tropismN, p.tropismS, tropism_name
        ]
    else:
        leaf_data["leaf"] = None


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

