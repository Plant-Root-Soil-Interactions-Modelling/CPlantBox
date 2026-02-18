"""conversions between dash store data types  and gui, D. Leitner 2026"""

import os

import numpy as np
from dash import html
from vtk.util import numpy_support
from vtk_conversions import *

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp

tropism_names = {
    "Plagiotropism": 0,
    "Gravitropism": 1,
    "Exotropism": 2,
    "Negative gravitropism": 4,
    "Variable gravitropism": 6,
}  # "Twist": 5,

tropism_names_ = {
    0: "Plagiotropism",
    1: "Gravitropism",
    2: "Exotropism",
    4: "Negative gravitropism",
    6: "Variable gravitropism",
}  # 5: "Twist",


def get_parameter_names():  # parameter xml file names
    """returns a list of plant parameter names with two values each, first a short name, second exact filename"""
    parameter_names = [
        ("Maize", "P0.xml"),
        ("Wheat", "Triticum_aestivum_test_2021.xml"),  # Monas File
        ("FSPM", "fspm2023.xml"),
        ("Demo Leaf", "leaf_only.xml"),
        ("Demo Root", "root_only.xml"),
        ("Demo Stem", "stem_only.xml"),
        ("User Data", "xml-store"),
    ]
    return parameter_names


def get_seed_slider_names():  # see set_data, apply_sliders
    """return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "First shoot borne root [day]": (1, 30, 0.1),
        "Shoot borne delay [day]": (1, 21, 0.1),
        "First basal root [day]": (1, 21, 0.1),
        "Basal root delay [day]": (1, 21, 0.1),
        "Maximal number of basal roots [1]": (0, 30, 1),
        "First tiller [day]": (1, 21, 0.1),
        "Tiller delay [day]": (1, 21, 0.1),
        "Maximal number of tillers [1]": (0, 7, 1),
    }
    return parameter_sliders


def get_root_slider_names():  # see set_data, apply_sliders
    """return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 200, 0.1),
        "Growth rate [cm/day]": (0.5, 10, 0.01),
        "Initial angle [°]": (0.0, 90, 0.1),
        "Basal zone [cm]": (0.1, 20, 0.1),
        "Interlateral distance [cm]": (0.1, 20, 0.1),
        "Apical zone [cm]": (0.1, 20, 0.1),
        "Radius [cm]": (1.0e-3, 0.25, 1e-3),
        "Tropism strength [1]": (0.0, 6.0, 0.1),
        "Tropism tortuosity [1]": (0.0, 1.0, 0.01),
    }
    return parameter_sliders


def get_stem_slider_names():  # see set_data, apply_sliders
    """return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 200, 0.1),
        "Growth rate [cm/day]": (0.5, 10, 0.01),
        "Initial angle [°]": (0.0, 90, 0.1),
        "Phytomer distance [cm]": (0.1, 20, 0.1),
        "Radius [cm]": (1.0e-3, 0.25, 1e-3),
        "Nodal growth start [day]": (0, 21, 0.1),
        "Nodal growth time span [day]": (0, 21, 0.1),
        "Fixed Rotation [°]": (0, 180, 0.1),
        "Tropism strength [1]": (0, 6, 0.1),
        "Tropism tortuosity [1]": (0, 1, 0.01),
    }
    return parameter_sliders


def get_leaf_slider_names():  # see set_data, apply_sliders
    """return slider names as keys of dict and bounds as values"""
    parameter_sliders = {
        "Maximal length [cm]": (1, 50, 0.1),
        "Growth rate [cm/day]": (0.5, 10, 0.01),
        "Initial angle [°]": (0.0, 180, 0.1),
        "Petiole length [cm]": (0.1, 10, 0.1),  # lb
        "Fixed Rotation [°]": (0, 180, 0.1),
        "Tropism strength [1]": (0, 6, 0.1),
        "Tropism tortuosity [1]": (0, 1, 0.01),
    }
    return parameter_sliders


def set_leaf_geometry(shapename, p):
    """shape name "Defined", "Long", "Round", "Maple", "Flower" """
    if shapename == "Defined":
        return
    elif shapename == "Long":
        p.lb = 1  # length of leaf stem
        p.la, p.lmax = 3.5, 8.5
        p.areaMax = 10  # cm2, area reached when length = lmax
        phi = np.array([-90, -45, 0.0, 45, 90]) / 180.0 * np.pi
        l = np.array([3, 2.2, 1.7, 2, 3.5])  # distance from leaf center
    elif shapename == "Round":
        p.lb = 1  # length of leaf stem
        p.la, p.lmax = 5, 11
        p.areaMax = 3.145 * (p.la**2)
        phi = np.array([-90, -45, 0.0, 45, 67.5, 70, 90]) / 180.0 * np.pi
        l_ = (p.lmax - p.lb) / 2  # == p.la
        l = np.array([l_ for x_i in range(len(phi))])  # ([2, 2, 2, 4,1,1, 4]) #distance from leaf center
    elif shapename == "Maple":
        p.lb = 1  # length of leaf stem
        p.areaMax = 50
        p.la, p.lmax = 5, 11
        phi = np.array([-90, -45, 0.0, 45, 67.5, 70, 90]) / 180.0 * np.pi
        l = np.array([5, 2, 2, 4, 1, 1, 5])  # distance from leaf center
    elif shapename == "Flower":
        p.lb = 1  # length of leaf stem
        p.areaMax = 50
        p.la, p.lb, p.lmax = 5, 1, 11
        phi = np.array([-90.0, -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]) / 180.0 * np.pi
        l = np.array([5.0, 1, 4.5, 1, 4, 1, 4.5, 1, 5])
    p.createLeafRadialGeometry(phi, l, 100)


def set_leaf_sliders(slider_values):
    """TODO cannot get the callbacks running in dash"""
    return slider_values


def fix_dx(rrp, strp, lrp):
    """overrides the xml resolution settings to ensure smooth visualization"""
    for r in rrp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.5
        r.dxMin = 1.0e-6
    for r in strp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.5
        r.dxMin = 1.0e-6
        r.nodalGrowth = 1  # <------- !
        # r.initBeta = 0.
        # r.betaDev = 0.
        # print("delayNGStarts", r.delayNGStarts)
        # print("delayNGEnds", r.delayNGEnds)
        # print("initBeta", r.initBeta)
        # print("betaDev", r.betaDev)
        # print("rotBeta", r.rotBeta)
        # r.rotBeta = 0.5
        # r.betaDev = 10
    for r in lrp:
        # print(r.subType, ":", r.dx, r.dxMin)
        r.dx = 0.25
        r.dxMin = 1.0e-6
        r.initBeta = 0.0
        r.betaDev = 0.0
        r.rotBeta = 0.5


def apply_sliders(srp, seed_data, rrp, root_data, strp, stem_data, lrp, leaf_data):
    """slider values to cplantbox random parameters"""
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
        srp.firstSB = 1.0e6  # disable
        srp.delaySB = 1.0e6
    if not seed_data["basal-checkbox"]:
        srp.firstB = 1.0e6  # disable
    if not seed_data["tillers-checkbox"]:
        srp.firstTil = 1.0e6  # disable
    # root
    for i, p in enumerate(rrp[1:]):
        # print(p.name, p.subType)
        d = root_data[f"tab-{i+1}"]
        p.lmax = d[0]
        p.r = d[1]
        p.theta = d[2] / 180.0 * np.pi
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
        p.theta = d[2] / 180.0 * np.pi
        p.ln = d[3]
        p.a = d[4]
        p.delayNGStart = d[5]
        p.delayNGEnd = d[5] + d[6]
        p.rotBeta = d[7] / 180.0
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
        p.rotBeta = d[5] / 180.0
        p.tropismN = d[6]
        p.tropismS = d[7]
        p.tropismT = tropism_names[d[8]]


def set_data(plant_, seed_data, root_data, stem_data, leaf_data, typename_data, xml_data):
    """sets all store data from the xml file"""
    # print("set_data()")
    typename_data.clear()
    """ open xml """
    fname = get_parameter_names()[int(plant_)][1]
    plant = pb.Plant()
    if fname == "xml-store":
        print("xlm data", type(xml_data["xml"]), len(xml_data["xml"]))
        plant.readParameters(xml_data["xml"], "plant", False, True)  # read from string, with verbose output (for debugging)
    else:
        my_dir = os.path.dirname(os.path.abspath(__file__))  # still works, if started from ohter folder
        plant.readParameters(my_dir + "/params/" + fname)
    """ seed """
    srp = plant.getOrganRandomParameter(pb.seed)
    p = srp[0]
    seed_data["seed"] = [
        p.firstSB,
        p.delaySB,
        p.firstB,
        p.delayB,
        p.maxB,
        p.firstTil,
        p.delayTil,
        p.maxTil,
    ]
    seed_data["shoot-checkbox"] = p.firstSB < 1.0e3 and p.delaySB < 1.0e3
    seed_data["basal-checkbox"] = p.firstB < 1.0e3 and p.delayB < 1.0e3 and p.maxB > 0
    seed_data["tillers-checkbox"] = p.firstTil < 1.0e3 and p.delayTil < 1.0e3 and p.maxTil > 0
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
            p.lmax,
            p.r,
            p.theta / np.pi * 180,
            p.lb,
            p.ln,
            p.la,
            p.a,
            p.tropismN,
            p.tropismS,
            tropism_name,
            len(p.successorST) > 0,
        ]
        typename_data[f"root tab-{i+1}"] = p.name
    """ stem """
    strp = plant.getOrganRandomParameter(pb.stem)
    typename_data["number_stemtypes"] = len(strp[1:])
    for i, p in enumerate(strp[1:]):
        tropism_name = tropism_names_[int(p.tropismT)]
        stem_data[f"tab-{i+1}"] = [
            p.lmax,
            p.r,
            p.theta / np.pi * 180,
            p.ln,
            p.a,
            p.delayNGStart,
            p.delayNGEnd - p.delayNGStart,
            p.rotBeta * 180,
            p.tropismN,
            p.tropismS,
            tropism_name,
            len(p.successorST) > 0,  # for p.rotBeta = 1 == 180 ° (???)
        ]
        typename_data[f"stem tab-{i+1}"] = p.name
    """ leaf """
    lrp = plant.getOrganRandomParameter(pb.leaf)
    if len(lrp) > 1:
        p = lrp[1]
        tropism_name = tropism_names_[int(p.tropismT)]
        leaf_data["leaf"] = [
            "Defined",
            p.lmax,
            p.r,
            p.theta / np.pi * 180,
            p.lb,
            p.rotBeta * 180,
            p.tropismN,
            p.tropismS,
            tropism_name,
        ]
    else:
        leaf_data["leaf"] = None


def into_panel(content, range_=None):
    """puts content into a panel"""
    if range_ is None:
        range_ = range(len(content))
    c_ = []
    for i in range_:
        c_.append(content[i])
    panel = html.Div(c_, className="panel")
    return panel
