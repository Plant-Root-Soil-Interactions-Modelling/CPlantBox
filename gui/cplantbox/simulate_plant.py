"""simulates the plant xml parameter set with slider values, D. Leitner 2026"""

import os

import numpy as np
from conversions import *  # auxiliary stuff
from vtk.util import numpy_support
from vtk_conversions import *

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


def get_plant(plant_, seed_data, root_data, stem_data, leaf_data, xml_data):
    """returns the plant object with slider values applied"""
    print("get_plant()")
    # 1. open base xml
    fname = get_parameter_names()[int(plant_)][1]
    plant = pb.Plant()
    if fname == "xml-store":
        plant.readParameters(xml_data["xml"], "plant", False, True)  # read from string, with verbose output (for debugging)
    else:
        my_dir = os.path.dirname(os.path.abspath(__file__))  # still works, if started from ohter folder
        plant.readParameters(my_dir + "/params/" + fname)
    srp = plant.getOrganRandomParameter(pb.seed)
    rrp = plant.getOrganRandomParameter(pb.root)
    strp = plant.getOrganRandomParameter(pb.stem)
    lrp = plant.getOrganRandomParameter(pb.leaf)
    fix_dx(rrp, strp, lrp)
    number_r = len(rrp[1:])  # number of root types
    number_s = len(strp[1:])  # number of stem types
    # 2. apply sliders to params
    apply_sliders(srp[0], seed_data, rrp, root_data, strp, stem_data, lrp, leaf_data)
    srp[0].seedPos.x = 0.0  # override position (always)
    srp[0].seedPos.y = 0.0
    srp[0].seedPos.z = -3.0
    srp[0].delayRC = 30.0
    srp[0].delayDefinitionShoot = 2
    print("get_plant() - delaySB", srp[0].delaySB)
    print("get_plant() - firstSB", srp[0].firstSB)
    print("get_plant() - delayRC", srp[0].delayRC)
    print("get_plant() - nC", srp[0].nC)
    return plant, number_r, number_s


def get_xml(plant_, seed_data, root_data, stem_data, leaf_data, xml_data):
    """returns the plant xml parameter set"""
    print("get_xml()")
    plant, _, _ = get_plant(plant_, seed_data, root_data, stem_data, leaf_data, xml_data)
    return plant.writeParameters("", "plant", False, True)  # write into string, with comments


def simulate_plant(plant_, time_slider, seed_data, root_data, stem_data, leaf_data, random_seed, xml_data):
    """simulates the plant xml parameter set with slider values"""
    print("simulate_plant()")

    # 1. open base xml and 2. apply sliders to params
    plant, number_r, number_s = get_plant(plant_, seed_data, root_data, stem_data, leaf_data, xml_data)

    # 3. simulate
    N = time_slider  # makes dt = 1
    t_ = np.linspace(0.0, time_slider, N + 1)
    root_length = np.zeros((number_r, N))
    stem_length = np.zeros((number_s, N))
    leaf_length = np.zeros((N,))
    plant.setSeed(random_seed)
    plant.initialize()
    rld_, z_, time_ = [], [], []
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

        # 4. make depth profiles
        if i + 1 in N // 5 * np.array(range(1, 6)):  # every 1/6 - 5/6 of the simulation time
            ana = pb.SegmentAnalyser(plant)
            length = ana.getSummed("length")
            max_ = ana.getMaxBounds().z
            min_ = ana.getMinBounds().z
            time_.append(i + 1)
            rld_.append(np.array(ana.distribution("length", max_, min_, int(np.round(max_ - min_)), True)))
            z_.append(np.linspace(max_, min_, int(np.round(max_ - min_))))

    # 5. make results store compatible (store pd stuff need for vtk.js)
    pd = vp.segs_to_polydata(plant, 1.0, ["subType", "organType", "radius", "creationTime"])  # poly-data, "radius",
    tube = apply_tube_filter(pd)  # polydata + tube filter
    vtk_data = vtk_polydata_to_dashvtk_dict(tube)  # addd "points" and "polys"
    cellData = pd.GetCellData()
    cT = numpy_support.vtk_to_numpy(cellData.GetArray("creationTime"))
    vtk_data["creationTime"] = encode_array(cT)  ################## TODO somehow creationTime and Age are mixed up
    # vtk_data["age"] = np.ones(cT.shape) * time_slider_value - cT
    organType = numpy_support.vtk_to_numpy(cellData.GetArray("organType"))  #
    vtk_data["subType"] = encode_array(numpy_support.vtk_to_numpy(cellData.GetArray("subType")) + 5 * (organType - np.ones(organType.shape) * 2))
    vtk_data["radius"] = encode_array(numpy_support.vtk_to_numpy(cellData.GetArray("radius")))

    # leaf geometry
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    leafes = plant.getOrgans(pb.leaf)
    for l in leafes:
        vp.create_leaf_(l, leaf_points, leaf_polys)
    pts_array = vtk.util.numpy_support.vtk_to_numpy(leaf_points.GetData()).astype(np.float32)
    polys_data = vtk.util.numpy_support.vtk_to_numpy(leaf_polys.GetData())
    vtk_data["leaf_points"] = encode_array(pts_array)
    vtk_data["leaf_polys"] = encode_array(polys_data)

    # 6. add dynamic results to result_data
    result_data = {}
    result_data["time"] = encode_array(t_[1:])
    for i in range(0, len(rld_)):
        result_data[f"rld{i}"] = encode_array(rld_[i] / length)
        result_data[f"z{i}"] = encode_array(z_[i])
        result_data[f"time{i}"] = time_[i]
    for j in range(number_r):
        result_data[f"root_length-{j+1}"] = encode_array(root_length[j, :])
    for j in range(number_s):
        result_data[f"stem_length-{j+1}"] = encode_array(stem_length[j, :])
    result_data["leaf_length"] = encode_array(leaf_length[:])
    # general
    result_data["number_r"] = number_r
    result_data["number_s"] = number_s

    return vtk_data, result_data
