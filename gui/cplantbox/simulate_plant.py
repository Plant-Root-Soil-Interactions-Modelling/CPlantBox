"""simulates the plant xml parameter set with slider values, D. Leitner 2026"""

import os

import numpy as np
from conversions import *  # auxiliary stuff
from vtk.util import numpy_support
from vtk_conversions import *

import plantbox as pb
import plantbox.visualisation.vtk_plot as vp


def _fit_values_to_length(values, target_len):
    """Resize scalar values to the target cell count with repeat+trim fallback."""
    values = np.asarray(values)
    if target_len <= 0:
        return values[:0]
    if values.size == target_len:
        return values
    if values.size == 0:
        return np.zeros(target_len)
    rep = int(np.ceil(target_len / values.size))
    return np.repeat(values, rep)[:target_len]


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
    srp[0].delayRC = 30.0  # root crown implementation is DEPRECATED (fix values)
    srp[0].simtime = min(45, srp[0].simtime)  # limit with slider max value
    print("simtime", srp[0].simtime)
    # srp[0].delayDefinitionShoot = 2
    # print("get_plant() - delaySB", srp[0].delaySB)
    # print("get_plant() - firstSB", srp[0].firstSB)
    # print("get_plant() - delayRC", srp[0].delayRC)
    # print("get_plant() - nC", srp[0].nC)
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
    pd_cell_data = pd.GetCellData()
    tube_cell_data = tube.GetCellData()
    n_tube_cells = tube.GetNumberOfCells()

    cT_arr = tube_cell_data.GetArray("creationTime")
    if cT_arr is not None:
        cT = numpy_support.vtk_to_numpy(cT_arr)
    else:
        cT = _fit_values_to_length(numpy_support.vtk_to_numpy(pd_cell_data.GetArray("creationTime")), n_tube_cells)
    vtk_data["creationTime"] = encode_array(cT)

    tube_subtype_arr = tube_cell_data.GetArray("subType")
    tube_organtype_arr = tube_cell_data.GetArray("organType")
    if tube_subtype_arr is not None and tube_organtype_arr is not None:
        sub_type = numpy_support.vtk_to_numpy(tube_subtype_arr)
        organ_type = numpy_support.vtk_to_numpy(tube_organtype_arr)
    else:
        sub_type = _fit_values_to_length(numpy_support.vtk_to_numpy(pd_cell_data.GetArray("subType")), n_tube_cells)
        organ_type = _fit_values_to_length(numpy_support.vtk_to_numpy(pd_cell_data.GetArray("organType")), n_tube_cells)
    vtk_data["subType"] = encode_array(sub_type + 5 * (organ_type - np.ones(organ_type.shape) * 2))

    tube_radius_arr = tube_cell_data.GetArray("radius")
    if tube_radius_arr is not None:
        radius = numpy_support.vtk_to_numpy(tube_radius_arr)
    else:
        radius = _fit_values_to_length(numpy_support.vtk_to_numpy(pd_cell_data.GetArray("radius")), n_tube_cells)
    vtk_data["radius"] = encode_array(radius)

    # leaf geometry
    leaf_points = vtk.vtkPoints()
    leaf_polys = vtk.vtkCellArray()  # describing the leaf surface area
    leafes = plant.getOrgans(pb.leaf)
    for l in leafes:
        vp.create_leaf_(l, leaf_points, leaf_polys)
    pts_array = numpy_support.vtk_to_numpy(leaf_points.GetData()).astype(np.float32)
    polys_data = numpy_support.vtk_to_numpy(leaf_polys.GetData())
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

    # print("simulate_plant():", vtk_data.keys(), result_data.keys())

    return vtk_data, result_data


def generate_mp4(plant_, time_slider, seed_data, root_data, stem_data, leaf_data, random_seed, xml_data, result_tab="VTK3D"):
    """Generate an MP4 animation (40 frames = 4 s at 10 fps) using offscreen VTK rendering.

    Simulates the plant twice: once to determine final scene bounds for a
    stable fixed camera, then again to capture 40 evenly-spaced frames.
    Returns a base64-encoded MP4 byte string suitable for dcc.Download with
    base64=True.
    """
    p_name = "age" if result_tab == "VTK3DAge" else "subType"
    import base64
    import shutil
    import subprocess
    import tempfile

    NUM_FRAMES = 40  # 4 s × 10 fps

    # ---- 1. First-pass simulation to determine final scene bounds ----
    plant_b, _, _ = get_plant(plant_, seed_data, root_data, stem_data, leaf_data, xml_data)
    plant_b.setSeed(random_seed)
    plant_b.initialize()
    N = max(1, int(time_slider))
    for dt_b in np.diff(np.linspace(0.0, time_slider, N + 1)):
        plant_b.simulate(dt_b)
    # Compute bounding-sphere the same way prepare_vtk_render_data does it in main.py
    pd_b = vp.segs_to_polydata(plant_b, 1.0, ["subType", "radius"])
    if pd_b.GetNumberOfPoints() > 0:
        pts_b = numpy_support.vtk_to_numpy(pd_b.GetPoints().GetData()).reshape(-1, 3)
        center = pts_b.mean(axis=0)
        radius = np.linalg.norm(pts_b - center, axis=1).max()
        radius = max(float(radius), 1.0)
    else:
        center = np.array([0.0, 0.0, -5.0])
        radius = 10.0
    distance = 1.5 * radius

    # ---- 2. Set up offscreen renderer with fixed oblique camera ----
    colors = vtk.vtkNamedColors()
    ren = vtk.vtkRenderer()
    ren.SetBackground(colors.GetColor3d("White"))
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(900, 750)
    renWin.SetOffScreenRendering(1)
    renWin.AddRenderer(ren)

    # Add the final-state actors temporarily so ResetCamera() can compute
    # the correct ParallelScale (same approach as render_window in vtk_plot.py)
    ana_final = pb.SegmentAnalyser(plant_b)
    ana_final.addAge(float(time_slider))
    tmp_actors, _ = vp.plot_plant(ana_final, p_name, render=False)
    for a in tmp_actors:
        ren.AddActor(a)
    ren.ResetCamera()
    for a in tmp_actors:
        ren.RemoveActor(a)

    camera = ren.GetActiveCamera()
    camera.ParallelProjectionOn()
    camera.SetFocalPoint(center.tolist())
    cam_pos = center + np.array([-distance, -distance, distance])
    camera.SetPosition(cam_pos.tolist())
    camera.SetViewUp(0, 0, 1)
    camera.OrthogonalizeViewUp()
    camera.SetClippingRange(0.1, distance * 10)

    # ---- 3. Second-pass simulation + frame capture ----
    plant_a, _, _ = get_plant(plant_, seed_data, root_data, stem_data, leaf_data, xml_data)
    plant_a.setSeed(random_seed)
    plant_a.initialize()

    tmpdir = tempfile.mkdtemp(prefix="cplantbox_mp4_")
    t_steps = np.linspace(0.0, time_slider, NUM_FRAMES + 1)
    cur_actors, cur_bar = [], None
    try:
        for i, (t0, t1) in enumerate(zip(t_steps, t_steps[1:])):
            plant_a.simulate(t1 - t0)
            ana = pb.SegmentAnalyser(plant_a)
            ana.addAge(t1)
            actors, scalar_bar = vp.plot_plant(ana, p_name, render=False)

            for a in cur_actors:
                ren.RemoveActor(a)
            if cur_bar is not None:
                ren.RemoveActor2D(cur_bar)
            for a in actors:
                ren.AddActor(a)
            if scalar_bar:
                ren.AddActor2D(scalar_bar)
            cur_actors, cur_bar = actors, scalar_bar

            vp.write_jpg(renWin, os.path.join(tmpdir, f"frame_{i}"), magnification=1)

        # ---- 4. Assemble MP4 with ffmpeg and return base64-encoded bytes ----
        mp4_path = os.path.join(tmpdir, "animation.mp4")
        subprocess.run(
            [
                "ffmpeg", "-y", "-r", "10",
                "-i", os.path.join(tmpdir, "frame_%d.jpg"),
                "-c:v", "libx264", "-pix_fmt", "yuv420p",
                mp4_path,
            ],
            check=True,
            capture_output=True,
        )
        with open(mp4_path, "rb") as f:
            return base64.b64encode(f.read()).decode("ascii")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
