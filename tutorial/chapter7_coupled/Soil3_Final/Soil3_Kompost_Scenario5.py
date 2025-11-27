"""
two independent plants, each with its own hydraulic + perirhizal model,
coupled to one DuMuX soil, using MultiPerirhizalManager to avoid
double-demand in cells where both plants have segments.

Now with 18 soil domains (6 in z, 3 in x) like the old rye setup:
z:  [0, -15, -30, -45, -60, -75, -90]
x:  [-50, -12], [-12, 12], [12, 50]
left/right = same VG, middle = other VG
"""

import sys; sys.path.append("../../.."); sys.path.append("../../../src/")
sys.path.append("../../../../dumux-rosi/build-cmake/cpp/python_binding/")
sys.path.append("../../../../dumux-rosi/python/modules/")

import plantbox as pb
from plantbox import New_MultiPerirhizalManager
import visualisation.vtk_plot as vp

from functional.PlantHydraulicParameters import PlantHydraulicParameters
from functional.PlantHydraulicModel import HydraulicModel_Doussan
from functional.PerirhizalHeterogeneous import PerirhizalHeterogeneous
import functional.van_genuchten as vg

from rosi_richards import RichardsSP
from richards import RichardsWrapper

import time
import timeit
import numpy as np
import matplotlib.pyplot as plt
import math
import json 
import pandas as pd
import os
import inspect 



# ---- Output-Pfade automatisch aus Skriptnamen ableiten ----
try:
    this_file = __file__
except NameError:
    # fallback, falls __file__ nicht gesetzt (z.B. in manchen Umgebungen)
    this_file = inspect.getsourcefile(lambda: 0)

script_name = os.path.splitext(os.path.basename(this_file))[0]  # z.B. "Soil3_Kompost_Baseline"
result_dir = os.path.join("results", script_name)
os.makedirs(result_dir, exist_ok=True)

# ein gemeinsamer Prefix für alle Dateien dieses Skripts
prefix = os.path.join(result_dir, script_name)
print("Saving results to:", result_dir)



# ------------------------------------------------------------
# Load treatment-specific k(z) profiles
#   Kompost (furrow / middle domains)
#   Kontrolle (off-furrow / left+right domains)
# ------------------------------------------------------------
with open("PR_calibration_kompost.json") as f:
    pr_calib_kom = json.load(f)

with open("PR_calibration_kontrolle.json") as f:
    pr_calib_kon = json.load(f)

k_depths_cm_kom = np.array(pr_calib_kom["depth_cm"], dtype=float)
k_profile_kom   = np.array(pr_calib_kom["k_profile"], dtype=float)

k_depths_cm_kon = np.array(pr_calib_kon["depth_cm"], dtype=float)
k_profile_kon   = np.array(pr_calib_kon["k_profile"], dtype=float)




# -------------------------------------------------------------------------
# helpers for mapping soil → root segments
# -------------------------------------------------------------------------
def build_cell_to_domain(soil, domain_ids):
    """
    Collect all DuMuX cells that belong to a set of soil 'domains'
    and return a dict: {cell_id: domain_id}.
    """
    cell2dom = {}
    for dom in domain_ids:
        cells = soil.getLayerCellIndx(dom)
        for c in cells:
            cell2dom[int(c)] = dom
    return cell2dom


def map_segments_to_domains(ms, cell2dom, default_domain=0):
    """
    Use CPB's ms.cell2seg to assign each segment the domain of its soil cell.
    """
    cell2seg = ms.cell2seg
    nseg = len(ms.radii)
    seg2dom = np.full(nseg, int(default_domain), dtype=int)
    for cell_id, segs in cell2seg.items():
        if cell_id in cell2dom:
            d = cell2dom[cell_id]
            for s in segs:
                seg2dom[s] = d
    return seg2dom


def map_segments_to_cell_value(ms, cell2value, default_value=1.0):
    """
    Generic version for things like VR_eff: per-cell value → per-segment value.
    """
    cell2seg = ms.cell2seg
    nseg = len(ms.radii)
    seg_values = np.full(nseg, float(default_value), dtype=float)
    for cell_id, segs in cell2seg.items():
        if cell_id in cell2value:
            v = cell2value[cell_id]
            for s in segs:
                seg_values[s] = v
    return seg_values


def sinusoidal(t):
    return np.sin(2. * np.pi * np.array(t) - 0.5 * np.pi) + 1.

def growth_factor(rs_age, full=60.0):
    g = rs_age / full
    return 1.0 if g >= 1.0 else g



# ------------------------------------------------------------
# domain / params  (rye-style box)
# ------------------------------------------------------------

# old baseline
min_b = [-48., -2., -90.]
max_b = [ 48.,  2.,   0.]
#cell_number = [50, 4, 90]

# grid resolution
#cell_number = [24, 1, 30] # fast testing  
cell_number = [24, 3, 30] # simulation target resolution

path = "../../../modelparameter/structural/rootsystem/"
#name = "wheat_1997_for_australia"
name = "hydro_wheat_1997_for_australia"

# 90 would be a sensible value for trans, roughly
trans = 30            # cm3/day per plant
wilting_point = -15000.   # cm
rs_age = 0.5            # days dont set to 0, then takeoff is too rough

# two "types" of soils (off-furrow vs at-furrow)
#soil3_sand  = [0.056, 0.372, 0.0162, 4.48, 218.0]  # off-furrow
#soil3_sand2 = [0.078, 0.43, 0.036, 1.56, 24.96]    # furrow/middle

initial = -70.
sim_time = 60     #Real observation of TRL was done on 02.06 (), with roughly 370cm on furrow, and 202 cm off-furrow, harvest was 15.07
dt = 150. / (24 * 3600)

# two plants, left/right # periodic tile is 1.0 m × 0.04 m with 8 plants → 8 / 0.04 = 200 plants/m²

plant_positions = [
    pb.Vector3d(-43.75, 0.0, -3.0),
    pb.Vector3d(-31.25, 0.0, -3.0),
    pb.Vector3d(-18.75, 0.0, -3.0),
    pb.Vector3d( -6.25, 0.0, -3.0),
    pb.Vector3d(  6.25, 0.0, -3.0),
    pb.Vector3d( 18.75, 0.0, -3.0),
    pb.Vector3d( 31.25, 0.0, -3.0),
    pb.Vector3d( 43.75, 0.0, -3.0),
]


length = min_b[0] * -1 + max_b[0]
width  = min_b[1] * -1 + max_b[1]
depth  = min_b[2] * -1 + max_b[2]
nx, ny, nz = cell_number
vr_lookup = pb.EquidistantGrid3D(length, width, depth, nx, ny, nz)


# ------------------------------------------------------------
# soil with 18 domains (6 z-slabs × 3 x-bands)
# ------------------------------------------------------------
s = RichardsWrapper(RichardsSP())
s.initialize()
s.createGrid(min_b, max_b, cell_number, periodic=True)
s.setHomogeneousIC(initial, False) # WICHTIG: Wenn nicht mit spin-up, besser hydrostatic equilibrium setzen!!


#s.setTopBC("noFlux")
s.setBotBC("constantPressure", -80.)

# READ PRECIPIPATION DATA
# read Excel
weather = pd.read_excel("Weather_data_CF3.xlsx", sheet_name="w_weather_data")
# define simulation time in days, starting at 0 for the first date / day
t0 = weather["date"].iloc[0]
times = (weather["date"] - t0).dt.days.astype(float).tolist()  # [0, 1, 2, 3, ...]

# convert precipitation from mm/day -> cm/day
E = 0.4   # cm/day dummy evaporation
fluxes = (weather["precipitation (mm/day)"] * 0.1 -0.4).tolist()

# build climate table for atmospheric BC
climate = [times, fluxes]

# set top BC (value_top is ignored for "atmospheric")
s.setTopBC("atmospheric", 0.0, climate)

# ------------------------------------------------------------
# DuMuX soil: 18 domains (6 depths × 3 columns)
# source params from Tobias' batch script (top→bottom there)
# ------------------------------------------------------------
# AF = auf Furche (middle) RICHTIGE MESSWERTE!!! NUR LAYER 0-30 in 0-15, 15-30 gespiegelt, da Messung nur auf 0-30 gemacht worden ist
af_params_top_to_bot = [
    [0.056, 0.372, 0.0162, 4.48, 218.0],   # layer0_af (TOP)
    [0.056, 0.372, 0.0162, 4.48, 218.0],   # layer1_af
    [0.047, 0.386, 0.0191, 3.82, 529.0],   # layer2_af
    [0.01307, 0.39421, 0.00624, 1.43671, 19.11943],   # layer3_af SCENARIO5
    [0.033, 0.340, 0.0127, 2.046, 96.0],   # layer4_af
    [0.019, 0.253, 0.0120, 1.922, 36.67],  # layer5_af (BOTTOM)
]

# NF = neben Furche (left + right)
nf_params_top_to_bot = [
    [0.048, 0.351, 0.0190, 4.887, 421.67],  # layer0_nf (TOP)
    [0.048, 0.351, 0.0190, 4.887, 421.67],  # layer1_nf
    [0.062, 0.371, 0.0160, 3.573, 373.3],   # layer2_nf
    [0.046, 0.368, 0.0152, 2.736, 105.67],  # layer3_nf
    [0.021, 0.304, 0.01590, 1.89, 127.67],  # layer4_nf
    [0.028, 0.243, 0.0124, 1.999, 7.89],    # layer5_nf (BOTTOM)
]

# our domain order now is BOTTOM→TOP,
# so depth 0 (bottom) must use index 5 from the lists above
vg_params = []
for depth in range(6):               # 0=bottom ... 5=top
    src_idx = 5 - depth              # flip because tables are top→bottom
    # left
    vg_params.append(nf_params_top_to_bot[src_idx])
    # middle
    vg_params.append(af_params_top_to_bot[src_idx])
    # right
    vg_params.append(nf_params_top_to_bot[src_idx])

# feed them to DuMuX in exactly this order
s.setVGParameters(vg_params)


############### VR_eff loooop

# Penetration_resistance und VR_eff --------------------------------
# Messdaten Bulk Density Roggen 2021
af_bulk_top_to_bot = [
    1.64,
    1.64,
    1.45,
    1.087,
    1.60,
    1.63,
]

nf_bulk_top_to_bot = [
    1.65,
    1.65,
    1.64,
    1.63,
    1.58,
    1.68,
]

# bulk density pro Domain – gleiche Reihenfolge wie vg_params
bulk_params = []
for depth in range(6):
    src = 5 - depth
    # left
    bulk_params.append(nf_bulk_top_to_bot[src])
    # middle
    bulk_params.append(af_bulk_top_to_bot[src])
    # right
    bulk_params.append(nf_bulk_top_to_bot[src])

# BD Normalisierung für Shehans Paper 
bd_min = 0.84
bd_max = max(bulk_params)
bd_span = bd_max - bd_min  # später auf 0 prüfen

z_breaks_domain = [0., -15., -30., -45., -60., -75., -90.]


def domain_center_depth_cm(depth_index):
    """
    depth_index: 0 = bottom ... 5 = top (same as in vg_params / bulk_params)
    returns positive depth in cm (0 at surface, increasing downward)
    """
    # bottom→top indexing vs top→bottom breakpoints:
    top_z = z_breaks_domain[5 - depth_index]   # e.g. -75 for depth_index=0
    bot_z = z_breaks_domain[6 - depth_index]   # e.g. -90
    top_d = abs(top_z)
    bot_d = abs(bot_z)
    return 0.5 * (top_d + bot_d)               # e.g. 82.5 cm

k_per_domain = []
for depth_index in range(6):   # 0 = bottom ... 5 = top
    z_center = domain_center_depth_cm(depth_index)  # positive depth in cm

    # interpolate Kompost and Kontrolle profiles separately
    k_here_kon = float(np.interp(
        z_center, k_depths_cm_kon, k_profile_kon,
        left=k_profile_kon[0], right=k_profile_kon[-1]
    ))
    k_here_kom = float(np.interp(
        z_center, k_depths_cm_kom, k_profile_kom,
        left=k_profile_kom[0], right=k_profile_kom[-1]
    ))

    # domain order: bottom→top, within each layer: left, middle, right
    # left  & right  = off-strip / Kontrolle
    # middle         = compost strip / Kompost
    k_per_domain.append(k_here_kon)  # left (off / Kontrolle)
    k_per_domain.append(k_here_kom)  # middle (Kompost)
    k_per_domain.append(k_here_kon)  # right (off / Kontrolle)

        
        

def theta_from_vg(h_cm, vg):
    """van Genuchten-θ(h) für eine Schicht."""
    theta_r, theta_s, alpha, n, Ks = vg
    if h_cm >= 0.0:
        return theta_s
    m = 1.0 - 1.0 / n
    Se = 1.0 / (1.0 + (alpha * abs(h_cm))**n)**m
    return theta_r + (theta_s - theta_r) * Se

# water content at wilting point per layer 
theta_p_per_domain = []
for dom_vg in vg_params:
    # wilting point bei -15000 cm
    theta_wp = theta_from_vg(-15000.0, dom_vg)
    theta_p_per_domain.append(theta_wp)


# Konstanten aus dem Paper
psi_w = -15000.0   # cm
PR_half = 2.0      # MPa


def compute_vr_eff_for_all_cells(s, cell2dom, debug=False):
    hs = s.getSolutionHead()      # matrix of pressure heads
    thetas = s.getWaterContent()  # matrix of volumetric water content
    vr_list = []

    PR_values = []  # for optional debug summary

    for cid in range(len(hs)):
        dom = cell2dom[cid]

        # --- moisture saturation term ---
        theta_v = thetas[cid]
        theta_p = theta_p_per_domain[dom]   # θ_p from VG(-15000)
        theta_s = vg_params[dom][1]         # θ_s

        denom = (theta_s - theta_p)
        if denom <= 0.0:
            sp = 0.0
        else:
            sp = (theta_v - theta_p) / denom
            sp = max(0.0, min(1.0, sp))

        # --- bulk density normalization ---
        bd_raw = bulk_params[dom]
        if bd_span > 0.0:
            bd = (bd_raw - bd_min) / bd_span
        else:
            bd = 0.0

        # --- PR (Vaz) * calibrated correction ---
        PR_raw = math.exp(1.5 + 2.18 * bd - 4.0 * sp)
        k_corr = k_per_domain[dom]
        PR = k_corr * PR_raw

        PR_values.append(PR)

        # --- VR_eff ---
        psi0 = hs[cid]
        vr_eff = -(psi0 / psi_w) + math.exp(-0.6931 * (PR / PR_half))
        vr_eff = max(0.0, min(1.0, vr_eff))

        vr_list.append(vr_eff)

    # --- Optional debug summary: PR statistics ---
    if debug and len(PR_values) > 0:
        PR_arr = np.array(PR_values)
        print(f"[DEBUG PR ALL] min={PR_arr.min():.2f} MPa, "
              f"mean={PR_arr.mean():.2f} MPa, "
              f"max={PR_arr.max():.2f} MPa")

        # -----------------------------
        # Per-layer × (compost/control)
        # -----------------------------
        # get cell depths in cm (positive downward)
        cc = s.getCellCenters()           # shape (ncells, 3)
        z = cc[:, 2]                      # z-coordinate
        depth_cm_all = np.abs(z)  # assumes domain in meters with z<0 downwards

        # layer bounds in cm
        z_layers = [0., 15., 30., 45., 60., 75., 90.]
        nL = len(z_layers) - 1

        # buckets: [layer][values...]
        compost_layers = [[] for _ in range(nL)]
        control_layers = [[] for _ in range(nL)]

        for cid, PR_val in enumerate(PR_values):
            dom = cell2dom[cid]
            dcm = depth_cm_all[cid]

            # find layer index
            layer_idx = None
            for li in range(nL):
                if z_layers[li] <= dcm < z_layers[li + 1]:
                    layer_idx = li
                    break
            if layer_idx is None:
                continue  # outside 0–90 cm, ignore

            if dom in furrow_dom_ids:       # compost strip
                compost_layers[layer_idx].append(PR_val)
            elif dom in off_dom_ids:        # control / off-furrow
                control_layers[layer_idx].append(PR_val)

        # print stats per layer
        for li in range(nL):
            z_top = z_layers[li]
            z_bot = z_layers[li + 1]

            if compost_layers[li]:
                arr_c = np.array(compost_layers[li])
                print(f"[DEBUG PR COMPOST {z_top:.0f}-{z_bot:.0f}cm] "
                      f"min={arr_c.min():.2f}, "
                      f"mean={arr_c.mean():.2f}, "
                      f"max={arr_c.max():.2f}")

            if control_layers[li]:
                arr_o = np.array(control_layers[li])
                print(f"[DEBUG PR CONTROL {z_top:.0f}-{z_bot:.0f}cm] "
                      f"min={arr_o.min():.2f}, "
                      f"mean={arr_o.mean():.2f}, "
                      f"max={arr_o.max():.2f}")



    return vr_list






#-----------------------------------------------------------
# z from top to bottom (for the breakpoints)
z_breaks = [0., -15., -30., -45., -60., -75., -90.]
# x bands like in rye
#x_breaks = [-50., -12., 12., 50.]
x_breaks = [-48.0, -12.0, 12.0, 48.0]


# Regularization tests

s.setParameter("Soil.SourceSlope", "3000")
s.setParameter("Problem.EnableLogging", "true")
s.setParameter("Newton.MaxSteps", "12")
s.setParameter("Newton.MaxRelativeShift", "1e-6")
s.setParameter("Newton.ResidualReduction", "1e-4")  # erstmal nicht 1e-2


s.initializeProblem()
#s.printParameters()

yL, yR = min_b[1], max_b[1]
dom_id = 0
for iz in range(5, -1, -1):        # bottom → top
    z_top = z_breaks[iz]
    z_bot = z_breaks[iz + 1]
    for ix in range(3):            # left → middle → right
        x_left  = x_breaks[ix]
        x_right = x_breaks[ix + 1]
        p1 = [x_left,  yL, z_bot]
        p2 = [x_right, yR, z_top]
        s.addVanGenuchtenDomain(p1, p2, dom_id)
        dom_id += 1

print("domain indexing scheme: bottom→top, and within each layer: left(0) → middle(1) → right(2)")
for d in range(18):
    cells = [int(c) for c in s.getLayerCellIndx(d)]
    if cells:
        print(f"domain {d}: {len(cells)} cells, id range [{min(cells)} … {max(cells)}]")
    else:
        print(f"domain {d}: 0 cells (!) check box definition")


s.setCriticalPressure(wilting_point)
#s.setParameter("SpatialParams.Swr", "0.15")
# build general cell→domain map for 18 domains
cell2dom = build_cell_to_domain(s, list(range(18)))


# ------------------------------------------------------------
# build two independent plants
# ------------------------------------------------------------
hydros = []
peris = []
mapped_plants = []      # keep original pb.MappedPlant for re-mapping
inner_r_all = []
rho_all = []


# Growth logging
vr_log = {
    "t": [],
    "furrow_mean": [],
    "off_mean": [],
    "dom_means": {}   # key = dom_id, value = list over time
}

furrow_dom_ids = [1, 4, 7, 10, 13, 16]   # alle mittleren
off_dom_ids    = [0,2, 3,5, 6,8, 9,11, 12,14, 15,17]



def log_vr_eff(t_cur, vr_eff_cells, cell2dom, furrow_dom_ids, off_dom_ids, vr_log):
    fur_vals = []
    off_vals = []

    # pro Dom sammeln
    dom_acc = {}   # dom -> [values]

    for cid, vr in enumerate(vr_eff_cells):
        dom = cell2dom[cid]
        dom_acc.setdefault(dom, []).append(vr)

        if dom in furrow_dom_ids:
            fur_vals.append(vr)
        elif dom in off_dom_ids:
            off_vals.append(vr)

    # Mittel über furrow/off
    fur_mean = sum(fur_vals) / len(fur_vals) if fur_vals else 0.0
    off_mean = sum(off_vals) / len(off_vals) if off_vals else 0.0

    vr_log["t"].append(t_cur)
    vr_log["furrow_mean"].append(fur_mean)
    vr_log["off_mean"].append(off_mean)

    # pro Dom Mittel ablegen
    for dom, vals in dom_acc.items():
        m = sum(vals) / len(vals)
        vr_log["dom_means"].setdefault(dom, []).append(m)

    # falls es Domains gibt, die in diesem step keine Zellen hatten (sollte nicht passieren),
    # kannst du hier noch ein padding machen





picker = lambda x, y, z: s.pick([x, y, z])

for i, pos in enumerate(plant_positions):
    rs = pb.MappedPlant()
    rs.readParameters(path + name + ".xml")
    """
    try:
            # 1. Get the LIST of all root parameter sets
            root_parameter_sets = rs.getOrganRandomParameter(pb.OrganTypes.root)
            
            # 2. Loop through the LIST and modify each parameter set
            for params in root_parameter_sets:
                params.dxMin = 0.05  # Set to a safe value (e.g., 0.1 cm)
                
            print(f"Plant {i}: Set dxMin=0.1 for {len(root_parameter_sets)} root types.")
            
    except Exception as e:
        print(f"Warning: Plant {i}: Could not set dxMin. Error: {e}")
    """    
    # make them different if you like
    rs.setSeed(1)

    # move seed
    organ_params = rs.getOrganRandomParameter(pb.OrganTypes.seed)
    organ_params[0].seedPos = pos

    # confine vertical domain (keep from your earlier script)
    base = pb.SDF_PlantBox(np.inf, np.inf, max_b[2] - min_b[2] - 0.6)
    sdf  = pb.SDF_RotateTranslate(base, pb.Vector3d(0., 0., -0.5))
    rs.setGeometry(sdf)

    # initial mapping to soil grid
    rs.setSoilGrid(picker)
    rs.setRectangularGrid(pb.Vector3d(min_b), pb.Vector3d(max_b),
                          pb.Vector3d(cell_number), False, True)

    rs.initialize(True)
    rs.simulate(rs_age, True)

    # hydraulics
    params = PlantHydraulicParameters()
    params.read_parameters("../../../modelparameter/functional/plant_hydraulics/couvreur2012")
    hm = HydraulicModel_Doussan(rs, params)
    hm.wilting_point = wilting_point

    # perirhizal with **18** soil domains
    peri = PerirhizalHeterogeneous(hm.ms)

    # paths to your 12 tables (adjust folder/suffix to your actual layout)
    af_tables_top_to_bot = [
        "Soil3_look_ups/layer0_af_mpi",
        "Soil3_look_ups/layer1_af_mpi",
        "Soil3_look_ups/layer2_af_mpi",
        "Soil3_scenario_look_ups/layer3_af_scenario5_mpi",
        "Soil3_look_ups/layer4_af_mpi",
        "Soil3_look_ups/layer5_af_mpi",
    ]
    nf_tables_top_to_bot = [
        "Soil3_look_ups/layer0_nf_mpi",
        "Soil3_look_ups/layer1_nf_mpi",
        "Soil3_look_ups/layer2_nf_mpi",
        "Soil3_look_ups/layer3_nf_mpi",
        "Soil3_look_ups/layer4_nf_mpi",
        "Soil3_look_ups/layer5_nf_mpi",
    ]

    # register 18 domains: bottom→top, left/mid/right
    for depth in range(6):                  # 0=bottom … 5=top
        src_idx = 5 - depth                 # flip top→bottom → bottom→top

        # left = NF of that depth
        peri.add_lookup_table(nf_tables_top_to_bot[src_idx])
        # middle = AF of that depth
        peri.add_lookup_table(af_tables_top_to_bot[src_idx])
        # right = NF again
        peri.add_lookup_table(nf_tables_top_to_bot[src_idx])

    
    
    """
    # register 18 domain parameter sets in the same left/mid/right pattern
    for dom in range(18):
        if dom % 3 == 1:  # middle columns
            peri.add_vg_parameters(vg.Parameters(soil3_sand2))
        else:             # left/right
            peri.add_vg_parameters(vg.Parameters(soil3_sand))
    """

    # initial seg → domain
    seg2dom = map_segments_to_domains(hm.ms, cell2dom, default_domain=0)
    peri.set_segment_soil_index(seg2dom)
    print("seg → soil:", np.bincount(seg2dom))

    # outer/inner radii
    outer_r = peri.get_outer_radii("length")
    inner_r = hm.ms.radii
    rho_ = np.divide(outer_r, np.array(inner_r))

    hydros.append(hm)
    peris.append(peri)
    mapped_plants.append(rs)
    inner_r_all.append(inner_r)
    rho_all.append(rho_)


# collect both plants into one analyser (for VTP)
plants = [hm.ms.plant() for hm in hydros]
ana = pb.SegmentAnalyser(plants[0])
for p in plants[1:]:
    ana.addSegments(p)

dx = max_b[0] - min_b[0]
dy = max_b[1] - min_b[1]
ana.mapPeriodic(dx, dy)

output_every = 10
output_idx = 0

# ------------------------------------------------------------
# shared perirhizal manager
# ------------------------------------------------------------
mgr = New_MultiPerirhizalManager()
for pidx, hm in enumerate(hydros):
    mgr.addPlant(peris[pidx], hm.ms)

mgr.recomputeSharedOuterRadii(2)

rho_all = []
for pidx, hm in enumerate(hydros):
    outer_r = np.array(mgr.getSharedOuter(pidx))
    inner_r = np.array(hm.ms.radii)
    rho_all.append(outer_r / inner_r)

# ------------------------------------------------------------
# time loop
# ------------------------------------------------------------
start_time = timeit.default_timer()

x_ = []
y_all = [[] for _ in hydros]
z_ = []

hs = s.getSolutionHead_()
hs_seg = [hm.ms.getHs(hs) for hm in hydros]
hsr_list = [seg.copy() for seg in hs_seg]
hl_log_per_plant_layer = []      # Stores positive values (Lift)
uptake_log_per_plant_layer = []  # Stores negative values (Source Extraction)

def print_water_profile(s, cell_number, min_b, max_b, cell2dom, vg_params, x_breaks):
    # areas
    width_cm  = max_b[0] - min_b[0]    # 96
    depth_cm  = max_b[1] - min_b[1]    # 4
    area_total_cm2 = width_cm * depth_cm
    area_total_m2  = area_total_cm2 / 1e4

    furrow_width_cm = x_breaks[2] - x_breaks[1]   # 24
    area_furrow_cm2 = furrow_width_cm * depth_cm  # 96
    area_furrow_m2  = area_furrow_cm2 / 1e4       # 0.0096
    area_off_m2     = area_total_m2 - area_furrow_m2

    # total from DuMuX
    init_water_cm3 = s.getWaterVolume()
    init_water_l   = init_water_cm3 / 1000.0
    init_water_lm2 = init_water_l / area_total_m2

    print(
        f"[water stats] total = {init_water_cm3:.2f} cm³ "
        f"({init_water_l:.3f} L, {init_water_lm2:.3f} L/m²)"
    )

    # depth-discretization
    hs = s.getSolutionHead_()
    nx, ny, nz = cell_number
    dx = (max_b[0] - min_b[0]) / nx
    dy = (max_b[1] - min_b[1]) / ny
    dz = abs((max_b[2] - min_b[2]) / nz)
    cell_vol_cm3 = dx * dy * dz

    bin_tops = list(range(0, 90, 10))
    nbins = len(bin_tops)
    water_tot_bins = [0.0] * nbins
    water_fur_bins = [0.0] * nbins
    water_off_bins = [0.0] * nbins

    def theta_vg(h_cm, params):
        qr, qs, alpha, n, ks = params
        if h_cm >= 0.0:
            Se = 1.0
        else:
            m = 1.0 - 1.0 / n
            Se = 1.0 / (1.0 + (alpha * abs(h_cm))**n)**m
        return qr + (qs - qr) * Se

    for cid, h in enumerate(hs):
        k = cid // (nx * ny)
        z_center = min_b[2] + (k + 0.5) * dz
        depth_cm = abs(z_center)
        bin_idx = int(depth_cm // 10)
        if bin_idx >= nbins:
            bin_idx = nbins - 1

        dom_id = cell2dom[cid]
        theta = theta_vg(h, vg_params[dom_id])
        w_cm3 = theta * cell_vol_cm3

        is_furrow = (dom_id % 3 == 1)
        water_tot_bins[bin_idx] += w_cm3
        if is_furrow:
            water_fur_bins[bin_idx] += w_cm3
        else:
            water_off_bins[bin_idx] += w_cm3

    print("[water stats] per 10 cm layer (total / furrow / off):")
    for i in range(nbins):
        z1 = i * 10
        z2 = z1 + 10
        w_tot_cm3 = water_tot_bins[i]
        w_fur_cm3 = water_fur_bins[i]
        w_off_cm3 = water_off_bins[i]

        w_tot_l = w_tot_cm3 / 1000.0
        w_fur_l = w_fur_cm3 / 1000.0
        w_off_l = w_off_cm3 / 1000.0

        w_tot_lm2 = w_tot_l / area_total_m2
        w_fur_lm2 = w_fur_l / area_furrow_m2 if area_furrow_m2 > 0 else 0.0
        w_off_lm2 = w_off_l / (area_total_m2 - area_furrow_m2) if area_total_m2 > area_furrow_m2 else 0.0

        print(
            f"  {z1:2d}-{z2:2d} cm: total {w_tot_cm3:8.2f} cm³ ({w_tot_l:5.3f} L, {w_tot_lm2:5.2f} L/m²) | "
            f"furrow {w_fur_cm3:8.2f} cm³ ({w_fur_l:5.3f} L, {w_fur_lm2:5.2f} L/m²) | "
            f"off {w_off_cm3:8.2f} cm³ ({w_off_l:5.3f} L, {w_off_lm2:5.2f} L/m²)"
        )

print_water_profile(s, cell_number, min_b, max_b, cell2dom, vg_params, x_breaks)



    
N = round(sim_time / dt)
t = 0.0

growth_every=15


# ---- depth bookkeeping for "uptake per z-layer" ----
nz = cell_number[2]
dz = (max_b[2] - min_b[2]) / nz
ztop = max_b[2]
depth_centers_cm = [
    ztop - (min_b[2] + (k + 0.5) * dz) for k in range(nz)
]

per_plant_depth_uptake_ts = []   # keep instantaneous profiles (optional)
cum_per_plant_depth_uptake = [np.zeros(nz) for _ in range(len(hydros))]



trl_log = [[] for _ in hydros]   # one list per plant
trl_time = []


for i in range(N):

    step_start = time.time()
    water_before = s.getWaterVolume()
    max_c = 0
    did_growth = False
    rad_time = 0.0
    growth_maint_time = 0.0
    plants_time = 0.0   # (or set before use)
    soil_time = 0.0

    groups = [
        (0, 7),  # boundary edge
        (1, 2),  # off-furrow left
        (3, 4),  # on furrow
        (5, 6),  # off-furrow right
    ]
    trl_groups = [0.0 for _ in groups]

    hs = s.getSolutionHead()
    hs_seg = [hm.ms.getHs(hs) for hm in hydros]

    # -------------------------------------------------
    # growth step
    # -------------------------------------------------
    if i > 0 and i % growth_every == 0:
        did_growth = True
        grow_dt = growth_every * dt

        # 1) VR_eff dynamisch aus dem Boden holen
        vr_eff_cells = compute_vr_eff_for_all_cells(s, cell2dom, debug=True)
        
        #log VR for plot
        sim_time_days = i * dt
        log_vr_eff(sim_time_days, vr_eff_cells, cell2dom,
                   furrow_dom_ids, off_dom_ids, vr_log)

        # 2) in das gemeinsame Grid schreiben
        vr_lookup.data = vr_eff_cells

        # 3) dieses eine Grid an alle Pflanzen hängen (robust)
        for rs in mapped_plants:
            root_params = rs.getOrganRandomParameter(pb.OrganTypes.root)
            if not root_params:
                continue
            for rp in root_params:
                if rp is None:
                    continue
                if hasattr(rp, "f_se"):
                    rp.f_se = vr_lookup

        # 3b) (optional) Mittelwert pro Pflanze ausgeben, nur zum Monitoring
        cell2vr = {cid: v for cid, v in enumerate(vr_eff_cells)}
        for pidx, hm in enumerate(hydros):
            seg_vr = map_segments_to_cell_value(hm.ms, cell2vr, default_value=1.0)
            mean_vr = float(np.mean(seg_vr))
            print(f"[VR_eff] plant {pidx}: mean VR_eff (segments) = {mean_vr:.3f} (nseg={len(seg_vr)})")

        # 4) wachsen lassen
        for hm in hydros:
            hm.ms.simulate(grow_dt, False)

        rad_start = time.time()
        mgr.recomputeSharedOuterRadii(2)
        rad_time = time.time() - rad_start
        print(f"[growth step {i}] recomputeSharedOuterRadii: {rad_time:.4f}s")

        # new soil head + segment heads
        hs = s.getSolutionHead()
        hs_seg = [hm.ms.getHs(hs) for hm in hydros]

        # time the heavy per-plant post-growth block
        growth_maint_start = time.time()

        trl_time.append(sim_time_days) # save time for TRL plot
        for pidx, hm in enumerate(hydros):
            new_inner = np.array(hm.ms.radii)
            new_outer = np.array(mgr.getSharedOuter(pidx))

            rho = new_outer / new_inner

            inner_r_all[pidx] = new_inner
            rho_all[pidx] = rho

            # resize stored hsr if plant grew
            if len(hsr_list[pidx]) != len(hs_seg[pidx]):
                old = hsr_list[pidx]
                new = np.empty(len(hs_seg[pidx]))
                new[:len(old)] = old
                new[len(old):] = hs_seg[pidx][len(old):]
                hsr_list[pidx] = new

            # TRL ...
            trl = float(np.sum(hm.ms.segLength()))
            trl_log[pidx].append(trl) #save TRL per plant for plot

            
            for g_idx, g in enumerate(groups):
                if pidx in g:
                    trl_groups[g_idx] += trl
                    break

            # re-map this plant to soil grid
            rs = mapped_plants[pidx]
            rs.setSoilGrid(picker)
            rs.setRectangularGrid(pb.Vector3d(min_b),
                                  pb.Vector3d(max_b),
                                  pb.Vector3d(cell_number),
                                  False,
                                  True)

            # rebuild seg → soil-domain for 18 domains
            seg2dom = map_segments_to_domains(hm.ms, cell2dom, default_domain=0)
            peris[pidx].set_segment_soil_index(seg2dom)

        growth_maint_time = time.time() - growth_maint_start

        print(f"[growth step {i}] TRL groups (cm):")
        print(f"  G1(0,7): {trl_groups[0]:.2f} (Target: 606 cm)")
        print(f"  G2(1,2): {trl_groups[1]:.2f} (Target: 606 cm)")
        print(f"  G3(3,4): {trl_groups[2]:.2f} (Target: 1100 cm)")
        print(f"  G4(5,6): {trl_groups[3]:.2f} (Target: 606 cm)")
        print(f"[growth step {i}] remap+seg2dom time: {growth_maint_time:.4f}s")


    else:
        # no growth, just get new soil head
        hs = s.getSolutionHead()
        hs_seg = [hm.ms.getHs(hs) for hm in hydros]
    

    # -------------------------------------------------
    # plant hydraulics
    # -------------------------------------------------
    plant_propose_start = time.time()
    hx_all = []

    # scaling attempts transpiration
    day_shape = sinusoidal(t)               # 0..2
    size_shape = growth_factor(rs_age + t)  # 0..1
    Tpot_step = trans * day_shape * size_shape  # cm3/day per plant


    for pidx, hm in enumerate(hydros):
        peri = peris[pidx]
        inner_r = inner_r_all[pidx]
        rho_ = rho_all[pidx]
        hs_seg_p = hs_seg[pidx]
        hsr = hsr_list[pidx]

        # first plant solve with scaled Tpot
        hx = hm.solve(rs_age + t, -Tpot_step, hsr, cells=False)
        hx_old = hx.copy()

        # segment conductivities
        kr_ = hm.params.getKr(rs_age + t)
        inner_kr_ = np.multiply(inner_r, kr_)

        # iterate plant ↔ perirhizal until change small
        err = 1.e6
        c = 0
        while err > 100. and c < 100:
            # update interface heads from perirhizal model
            hsr = peri.soil_root_interface_potentials(
                hx[1:],    # xylem head per segment
                hs_seg_p,  # soil head at segment
                inner_kr_,
                rho_,
            )
            # solve plant again with updated interface heads
            hx = hm.solve_again(rs_age + t, -Tpot_step, hsr, cells=False)
            err = np.linalg.norm(hx - hx_old)
            hx_old = hx.copy()
            c += 1

        # remember worst iteration count
        max_c = max(max_c, c)

        # compute ranges for this plant
        hsr_min = float(np.min(hsr))
        hsr_max = float(np.max(hsr))
        hx_min = float(np.min(hx))
        hx_max = float(np.max(hx))

        # give merged manager our final interface heads
        mgr.setProposal(pidx, hsr)

        # store for next time step
        hsr_list[pidx] = hsr
        hx_all.append(hx)

    plant_propose_time = time.time() - plant_propose_start



    # -------------------------------------------------
    # merge
    # -------------------------------------------------
    #merge_start = time.time()
    #mgr.mergeSharedCells()
    #merge_time = time.time() - merge_start

    # -------------------------------------------------
    # flux + soil
    # -------------------------------------------------
    plant_flux_start = time.time()
    total_source = {}
    per_plant_srcs = []

    for pidx, hm in enumerate(hydros):
        hsr_merged = hsr_list[pidx]

        # re-solve xylem with the merged interface
        hx = hm.solve_again(rs_age + t, -Tpot_step, hsr_merged, cells=False)
        hx_all[pidx] = hx  # keep it in sync for VTK / later use

        # now flux with the SAME pair (hx, hsr_merged)
        fluxes = hm.radial_fluxes(rs_age + t, hx, hsr_merged, cells=False)
        src = hm.sumSegFluxes(fluxes)
        per_plant_srcs.append(src)

        # das MUSS bleiben: alle Pflanzen-Sources aufsummieren
        for cid, val in src.items():
            total_source[cid] = total_source.get(cid, 0.0) + val

        # transpiration loggen wie vorher
        y_all[pidx].append(
            hm.get_transpiration(rs_age + t, hx.copy(), hsr_merged.copy())
        )
        hsr_list[pidx] = hsr_merged

    plant_flux_time = time.time() - plant_flux_start
    plants_time = plant_propose_time + plant_flux_time

    ### NEW: aggregate REAL per-plant uptake per actual z-layer
    nx, ny, nz = cell_number
    plant_depth_uptake = [np.zeros(nz) for _ in range(len(hydros))]

    for pidx, src in enumerate(per_plant_srcs):
        arr = plant_depth_uptake[pidx]
        for cid, q in src.items():
            k = cid // (nx * ny)      # which z-layer this cell is in
            arr[k] += -q              # make uptake positive (cm³/day)

    # dt is in days → arr is cm³/day → arr * dt = cm³ in this step
    for pidx, arr in enumerate(plant_depth_uptake):
        cum_per_plant_depth_uptake[pidx] += arr * dt

    # optional: still keep the per-timestep profile
    per_plant_depth_uptake_ts.append(plant_depth_uptake)

    
    # ---- per-plant uptake by depth bands (using global z_breaks) ----
    # depth bands (already have z_breaks = [0., -15., ...])
    depth_breaks = [abs(z) for z in z_breaks]
    nb = len(depth_breaks) - 1

    nx, ny, nz = cell_number
    ztop = max_b[2]
    dz = (max_b[2] - min_b[2]) / nz

    # per plant srcs already collected above
    # per_plant_srcs[pidx] = {cell_id: q, ...}

    # 4 groups: (0,7), (1,2), (3,4), (5,6)
    groups = [
        (0, 7),
        (1, 2),
        (3, 4),
        (5, 6),
    ]

    uptake = [[0.0]*len(groups) for _ in range(nb)]

    for pidx, src in enumerate(per_plant_srcs):
        for cid, q in src.items():
            k = cid // (nx*ny)
            z_center = min_b[2] + (k + 0.5)*dz
            depth = ztop - z_center  # cm below top

            # which depth band?
            for b in range(nb):
                if depth_breaks[b] <= depth < depth_breaks[b+1]:
                    # which group is this plant in?
                    for g_idx, g in enumerate(groups):
                        if pidx in g:
                            uptake[b][g_idx] += -q   # make uptake positive
                            break
                    break

    # -------------------------------------------------
    # Hydraulic Flux Analysis: Lift vs. Uptake
    # -------------------------------------------------
    nx, ny, nz = cell_number
    
    # Create TWO matrices for this step
    step_lift   = np.zeros((len(hydros), nz)) 
    step_uptake = np.zeros((len(hydros), nz))
    
    for pidx, src in enumerate(per_plant_srcs):
        for cid, q in src.items():
            k = cid // (nx * ny)          # Z-layer
            vol = q * dt                  # Volume in cm3
            
            if q > 0.0:
                # POSITIVE: Add to Lift Matrix
                step_lift[pidx][k] += vol
            else:
                # NEGATIVE: Add to Uptake Matrix (keep sign negative)
                step_uptake[pidx][k] += vol
    
    hl_log_per_plant_layer.append(step_lift)
    uptake_log_per_plant_layer.append(step_uptake)

    # Monitor strictly the LIFT (Positive addition to soil)
    net_lift_rate = np.sum(step_lift) / dt  # Back to rate for printing
    if net_lift_rate > 0.01:
        print(f"[Hydraulic Lift] Time {t:.3f}: Net release {net_lift_rate:.4f} cm3/day into soil!")
    # -------------------------------------------------


    # ---- soil solve ----
    soil_start = time.time()
    s.setSource(total_source)
    s.solve(dt, saveInnerFluxes_ = True)
    soil_time = time.time() - soil_start


    if i % 10 == 0:
        print(f"\n=== water stats at step {i}, t = {t:.3f} d ===")
        print_water_profile(s, cell_number, min_b, max_b, cell2dom, vg_params, x_breaks)

    # fresh soil (like baseline)
    hs = s.getSolutionHead()
    soil_min = float(np.min(hs))
    soil_max = float(np.max(hs))

    # aggregate over ALL plants to mimic single-plant behaviour
    all_iface_mins = [float(np.min(hsr_list[p])) for p in range(len(hsr_list))]
    all_iface_maxs = [float(np.max(hsr_list[p])) for p in range(len(hsr_list))]
    iface_min = min(all_iface_mins)
    iface_max = max(all_iface_maxs)

    all_xyl_mins = [float(np.min(hx_all[p])) for p in range(len(hx_all))]
    all_xyl_maxs = [float(np.max(hx_all[p])) for p in range(len(hx_all))]
    xyl_min = min(all_xyl_mins)
    xyl_max = max(all_xyl_maxs)

    # progress bar
    n = round(float(i) / float(N) * 100.)

    print("[{}{}], {} iterations, soil hs [{:.1f}, {:.1f}], "
          "interface [{:.1f}, {:.1f}] cm, root [{:.1f}, {:.1f}] cm, {:.6g} days"
          .format(
              "*" * n,
              " " * (100 - n),
              max_c,
              soil_min, soil_max,
              iface_min, iface_max,
              xyl_min, xyl_max,
              s.simTime,
          ))




    soil_water = (s.getWaterVolume() - water_before) / dt
    z_.append(soil_water)

    x_.append(t)
    t += dt

    loop_time = time.time() - step_start
    other_time = loop_time - plants_time - soil_time - growth_maint_time - rad_time


    if did_growth:
        print(
            f"[TIMING] step {i:6d}: "
            f"loop={loop_time:.4f}s, "
            f"soil={soil_time:.4f}s, "
            f"plants={plants_time:.4f}s, "
            f"growth_geom={growth_maint_time:.4f}s, "
            f"outerR={rad_time:.4f}s, "
            f"other={other_time:.4f}s"
        )


    n = round(float(i) / float(N) * 100.)
    print(f"[{'*' * n}{' ' * (100 - n)}] step {i+1}/{N}, "
          f"max per-plant iters = {max_c}, soil time = {s.simTime:.4f} d")

    if i % output_every == 0:
        vp.write_soil(
            f"{prefix}_soil_{output_idx:06d}",
            s,
            min_b,
            max_b,
            cell_number,
        )

        ana = pb.SegmentAnalyser()
        xylem_seg = []

        for pidx, hm in enumerate(hydros):
            ms = hm.ms                 # <-- use the multisegment system
            ana.addSegments(ms)        # geometry and segment order from ms
            hx = hx_all[pidx]          # node xylem pressures for THIS ms

            # --- SOLUTION ---
            # Assign the *average* pressure of the segment's
            # proximal (x) and distal (y) nodes.
            for seg in ms.segments:    # or ms.getSegments()
                nx = int(seg.x)  # Proximal node index
                ny = int(seg.y)  # Distal node index
                val = 0.5 * (float(hx[nx]) + float(hx[ny]))
                xylem_seg.append(val)

        ana.addData("xylem_pressure", xylem_seg)
        ana.mapPeriodic(dx, dy)
        ana.write(
            f"{prefix}_roots_{output_idx:06d}.vtp",
            ["radius", "subType", "creationTime", "age",
             "organType", "length", "xylem_pressure"],
        )

        output_idx += 1



print("done in", timeit.default_timer() - start_time, "s")

# ------------------------------------------------------------
# plotting
# ------------------------------------------------------------
# ------------------------------------------------------------
# plotting (grouped)
# ------------------------------------------------------------
tt = np.array(x_)

groups = [
    (0, 7),  # boundary edge
    (1, 2),  # off-furrow left
    (3, 4),  # on furrow
    (5, 6),  # off-furrow right
]
group_labels = [
    "G1 (0,7) boundary",
    "G2 (1,2) off-furrow L",
    "G3 (3,4) on furrow",
    "G4 (5,6) off-furrow R",
]

# prescribed per plant
pot_per_plant = (
    trans * sinusoidal(tt) * np.vectorize(growth_factor)(rs_age + tt)
)
plt.figure()
plt.plot(tt, pot_per_plant, "k", label="potential / plant")

# actual: gruppieren und MITTELN
for (g, label) in zip(groups, group_labels):
    g_series = np.zeros_like(tt, dtype=float)
    for pidx in g:
        g_series += np.array(y_all[pidx])
    g_series /= len(g)          # <-- hier das Mitteln
    plt.plot(tt, -g_series, label=label)

plt.legend()
plt.xlabel("time [d]")
plt.ylabel("transpiration [cm³/d]")
plt.tight_layout()
# save plot
plt.savefig(f"{prefix}_transpiration_groups.png", dpi=300)

plt.show()



# y_all[pidx] = [Tact_step0, Tact_step1, ...]  in cm³/day
plant_total_T = []
for pidx, series in enumerate(y_all):
    total = sum(val * dt for val in series)   # cm³ over whole sim
    plant_total_T.append(total)

group_defs = [
    ("on furrow", (3, 4)),
    ("left beside furrow", (1, 2)),
    ("right beside furrow", (5, 6)),
    ("furthest from furrow", (0, 7)),
]

for k, d in enumerate(depth_centers_cm):
    print(f"layer {k:2d}: {d:.2f} cm")


plt.figure()
for name, plants in group_defs:
    # sum soil-side depth profiles of this group
    profile = np.zeros_like(cum_per_plant_depth_uptake[0])
    soil_total = 0.0
    plant_total = 0.0
    for pidx in plants:
        profile += cum_per_plant_depth_uptake[pidx]
        soil_total += cum_per_plant_depth_uptake[pidx].sum()  # soil view
        plant_total += plant_total_T[pidx]                     # plant view

    label = f"{name} (soil={soil_total:.1f}, T={plant_total:.1f} cm³)"
    plt.plot(profile, depth_centers_cm, marker="o", label=label)

plt.gca().invert_yaxis()
plt.xlabel("accumulated uptake [cm³]")
plt.ylabel("depth [cm]")
plt.title("Accumulated uptake per depth layer by position")
plt.legend()
plt.tight_layout()
# --- save it ---
plt.savefig(f"{prefix}_uptake_depth_groups.png", dpi=300)
# ---------------

plt.show()

# ---- VR_eff pro Tiefenband plotten (oben → unten) ----
t = vr_log["t"]

# z_breaks kommen bei dir von oben nach unten: [0., -15., -30., ...]  :contentReference[oaicite:0]{index=0}
depth_breaks = [abs(z) for z in z_breaks]   # [0, 15, 30, 45, 60, 75, 90]
# daraus ergeben sich 6 Bänder:
# 0–15, 15–30, 30–45, 45–60, 60–75, 75–90

# Domains sind aber bottom → top angelegt:
# unten: 0,1,2 ... oben: 15,16,17  :contentReference[oaicite:1]{index=1}
# Mitte (furrow) je Schicht: 1,4,7,10,13,16  → das drehen wir jetzt um:
furrow_dom_ids_by_depth = [16, 13, 10, 7, 4, 1]   # top → bottom

# off-furrow je Schicht auch drehen:
off_dom_ids_by_depth = [
    [15, 17],   # 0–15 cm
    [12, 14],   # 15–30 cm
    [9, 11],    # 30–45 cm
    [6, 8],     # 45–60 cm
    [3, 5],     # 60–75 cm
    [0, 2],     # 75–90 cm
]

plt.figure()

# Furrow-Linien
for band_idx, dom_id in enumerate(furrow_dom_ids_by_depth):
    if dom_id in vr_log["dom_means"]:
        z_top = depth_breaks[band_idx]
        z_bot = depth_breaks[band_idx + 1]
        label = f"furrow {z_top:.0f}–{z_bot:.0f} cm"
        plt.plot(t, vr_log["dom_means"][dom_id], label=label)

# Off-furrow gemittelt
for band_idx, dom_pair in enumerate(off_dom_ids_by_depth):
    z_top = depth_breaks[band_idx]
    z_bot = depth_breaks[band_idx + 1]
    label = f"off {z_top:.0f}–{z_bot:.0f} cm"

    series = []
    for k in range(len(t)):
        vals = []
        for dom_id in dom_pair:
            if dom_id in vr_log["dom_means"]:
                vals.append(vr_log["dom_means"][dom_id][k])
        series.append(sum(vals) / len(vals) if vals else 0.0)

    plt.plot(t, series, linestyle="--", label=label)

plt.xlabel("time (d)")
plt.ylabel("VR_eff")
plt.legend()
plt.tight_layout()
plt.savefig(f"{prefix}_VR_eff_depths.png", dpi=300)
plt.show()



#TRL plot
plt.figure(figsize=(12, 8))

for pidx in range(len(hydros)):   # loop over plants directly
    plt.plot(trl_time, trl_log[pidx], label=f"Plant {pidx}")

plt.xlabel("Time (days)")
plt.ylabel("Total Root Length (cm)")
plt.title("TRL per Plant Over Time")
plt.legend()
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig(f"{prefix}_trl_over_time.png", dpi=300)
plt.show()


# ============================================================
# 1) TRL over time – groups only
# ============================================================

plt.figure(figsize=(10, 6))

# group-averaged TRL curves
for g, label in zip(groups, group_labels):
    g_trl = np.zeros(len(trl_time), dtype=float)
    for pidx in g:
        g_trl += np.array(trl_log[pidx])

    plt.plot(trl_time, g_trl, linewidth=2.5, label=label)

plt.xlabel("time [d]")
plt.ylabel("total root length [cm]")
plt.title("Total root length (TRL) per group")
plt.legend()
plt.tight_layout()
plt.savefig(f"{prefix}_trl_groups.png", dpi=300)
plt.show()


# ============================================================
# Root Length (absolute) & Root Length Density (RLD)
# ============================================================

nx, ny, nz = cell_number

dx = (max_b[0] - min_b[0]) / nx   # cm
dy = (max_b[1] - min_b[1]) / ny   # cm
dz = abs((max_b[2] - min_b[2]) / nz)

cell_vol_cm3 = dx * dy * dz
layer_vol_cm3 = nx * ny * cell_vol_cm3

# depth centers
ztop = max_b[2]
depth_centers = np.array([
    (ztop - (min_b[2] + (k + 0.5) * dz)) for k in range(nz)
])

# ---------- Root length per plant ----------
length_profiles = []
rld_profiles = []

for pidx, hm in enumerate(hydros):
    ms = hm.ms
    seg_lengths = np.array(ms.segLength())

    profile_len = np.zeros(nz, dtype=float)  # absolute length
    profile_rld = np.zeros(nz, dtype=float)  # length density

    cell2seg = ms.cell2seg
    for cid, segs in cell2seg.items():
        k = cid // (nx * ny)
        for s in segs:
            L = seg_lengths[s]
            profile_len[k] += L
            profile_rld[k] += L

    length_profiles.append(profile_len)                # cm per layer
    rld_profiles.append(profile_rld / layer_vol_cm3)   # cm/cm³

length_profiles = np.array(length_profiles)  # (plants, nz)
rld_profiles    = np.array(rld_profiles)     # (plants, nz)

# ---------- Per-group averages ----------
length_group_profiles = []
rld_group_profiles = []

for g in groups:
    lp = np.zeros(nz)
    rp = np.zeros(nz)

    for pidx in g:
        lp += length_profiles[pidx]
        rp += rld_profiles[pidx]

    length_group_profiles.append(lp)
    rld_group_profiles.append(rp)

length_group_profiles = np.array(length_group_profiles)
rld_group_profiles    = np.array(rld_group_profiles)






# Kontrolle: Summe über z vs. TRL aus trl_log
print("=== Check TRL vs Profil-Summen ===")
for g_idx, (g, label) in enumerate(zip(groups, group_labels)):
    # Summe der Längen im Profil
    prof_sum = np.sum(length_group_profiles[g_idx])

    # Summe der TRL der Pflanzen in dieser Gruppe (letzter Zeitschritt)
    trl_sum = 0.0
    for pidx in g:
        trl_sum += trl_log[pidx][-1]

    print(f"{label}: profil sum = {prof_sum:.2f} cm, TRL sum = {trl_sum:.2f} cm")


# ============================================================
# 1) Plot absolute root length profile (cm per depth layer)
# ============================================================

plt.figure(figsize=(8, 10))

for lp, label in zip(length_group_profiles, group_labels):
    plt.plot(lp, depth_centers, linewidth=2.5, marker="o" ,label=label)

plt.gca().invert_yaxis()
plt.xlabel("Root length per depth layer (cm)")
plt.ylabel("Depth (cm)")
plt.title("Root length profiles per group (final state)")
plt.legend()
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig(f"{prefix}_root_length_profiles_groups_final.png", dpi=300)
plt.show()


# ============================================================
# 2) Plot RLD profile (cm/cm³ per depth layer)
# ============================================================

plt.figure(figsize=(8, 10))

for rp, label in zip(rld_group_profiles, group_labels):
    plt.plot(rp, depth_centers, linewidth=2.5, label=label)

plt.gca().invert_yaxis()
plt.xlabel("Root length density (cm cm$^{-3}$)")
plt.ylabel("Depth (cm)")
plt.title("Root length density profiles per group (final state)")
plt.legend()
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig(f"{prefix}_rld_profiles_groups_final.png", dpi=300)
plt.show()


# ------------------------------------------------------------
# Vergleich: Modellprofil (G3 auf Furche) vs. Feldprofil FU+Ko_af
# ------------------------------------------------------------

# Messdaten (2 Pflanzen pro 20 cm each. Die Messungen beschreiben Zaehlungen in 4 grids (jeweils 5x5cm) mit 5 cm Hoehe. Wir nehmen wie ueblich eine Tiefe von 1cm an und Skalieren die Tiefe hier mit 3, um die Messwerte mit der Simulation vergleichbar zu machen.
depth_meas = np.array([5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90], dtype=float)
FU_Ko_zf = np.array([37.0,53.3,41.0,16.7,28.7,14.3,0.7,0.0,0.7,0.0,0.0,3.0,0.7,0.0,1.5,2.7,1.7,0.3])
FU_Ko_af = np.array([25.0,45.3,26.7,20.7,26.3,31.0,36.7,52.0,55.7,19.3,2.7,3.0,4.7,5.3,3.3,4.3,1.7,5.3])
FU_Ko_zf_corrected = FU_Ko_zf * 3 # to account for 3cm simdomain, ca 606  cm fuer zwei Pflanzen zur Bluete
FU_Ko_af_corrected = FU_Ko_af * 3 # to account for 3cm simdomain, ca 1100 cm fuer zwei Pflanzen zur Bluete


# G3 = (3,4) on furrow → Index 2 in groups/group_labels
g3_idx = 2
g2_idx = 1
sim_depth = depth_centers            # cm
sim_len1   = length_group_profiles[g3_idx]  # cm pro Layer, 2 Pflanzen
sim_len2   = length_group_profiles[g2_idx]  # cm pro Layer, 2 Pflanzen
plt.figure(figsize=(6, 10))

# Modellprofil
plt.plot(sim_len1, sim_depth, label="Modell G3 Furrow(2 Pflanzen)", linewidth=2.5)
plt.plot(sim_len2, sim_depth, label="Modell G2 off_furrow(2 Pflanzen)", linewidth=2.5)
# Feldprofil: auf Furche (FU+Ko_af)
plt.plot(FU_Ko_af_corrected, depth_meas, "o--", label="Feldmessung Furche (2 Pflanzen)")
plt.plot(FU_Ko_zf_corrected, depth_meas, "o--", label="Feldmessung off_furrow (2 Pflanzen)")
plt.gca().invert_yaxis()
plt.xlabel("Root length per depth layer [cm]")
plt.ylabel("Depth [cm]")
plt.title("Root length profile on furrow: Model vs Field")
plt.legend()
plt.grid(True, linestyle=":", alpha=0.6)
plt.tight_layout()
plt.savefig(f"{prefix}_compare_profile_on_furrow.png", dpi=300)
plt.show()


import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

# ============================================================
# SAVE DATA
# ============================================================
hl_data_array = np.array(hl_log_per_plant_layer)
uptake_data_array = np.array(uptake_log_per_plant_layer)

try:
    np.savez_compressed(
        f"{prefix}_timeseries_data.npz",
        sim_time_days=np.array(x_),
        transpiration_per_plant_cm3_day=np.array(y_all),
        trl_time_days=np.array(trl_time),
        trl_per_plant_cm=np.array(trl_log),
        cum_uptake_per_plant_per_layer_cm3=np.array(cum_per_plant_depth_uptake),
        layer_depth_centers_cm=np.array(depth_centers_cm),
        # Detailed 3D arrays
        hydraulic_lift_per_plant_layer_cm3=hl_data_array,
        uptake_per_plant_layer_cm3=uptake_data_array 
    )
    print("... All time series saved (including detailed Lift & Uptake).")
except Exception as e:
    print(f"Fehler beim Speichern: {e}")


# ============================================================
# PREPARE DATA & SCALING (Shared for both plots)
# ============================================================

# Combine for Net Flux: Blue=Lift(Pos), Red=Uptake(Neg)
net_flux_array = hl_data_array + uptake_data_array 

# Convert time list to array for plotting
time_arr = np.array(x_)

# 1. Determine Dynamic Scaling (TwoSlopeNorm)
# This forces 0 to be White, while scaling Red and Blue independently
# to their respective maximums.
data_min = np.min(net_flux_array)
data_max = np.max(net_flux_array)

# Safety checks to avoid crashes if simulation is empty or one-sided
if data_max <= 0: data_max = 1e-5      # Ensure we have a positive range
if data_min >= 0: data_min = -1e-5     # Ensure we have a negative range

print(f"Flux Scaling: Uptake(Red) down to {data_min:.4e}, Lift(Blue) up to {data_max:.4e}")

# Create the normalization object
norm = mcolors.TwoSlopeNorm(vmin=data_min, vcenter=0., vmax=data_max)

# Create meshgrid (Shared)
time_grid, depth_grid = np.meshgrid(time_arr, np.array(depth_centers_cm))


# ============================================================
# PLOT 1: Diverging Flux Heatmaps (FULL TIMELINE)
# ============================================================

fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(18, 6), sharey=True, sharex=True)

for idx, (grp, label) in enumerate(zip(groups, group_labels)):
    ax = axes[idx]
    
    # Select & Mean per group
    group_data = net_flux_array[:, grp, :]
    mean_flux = np.mean(group_data, axis=1).T # Transpose to (Layers, Time)
    
    # Plot using 'norm' instead of vmin/vmax
    c = ax.pcolormesh(time_grid, depth_grid, mean_flux, 
                      cmap='RdBu', norm=norm, shading='auto')
    
    ax.set_title(label, fontsize=10, fontweight='bold')
    ax.set_xlabel("Time [d]")
    if idx == 0: ax.set_ylabel("Depth [cm]")
    if not ax.yaxis_inverted(): ax.invert_yaxis()

# Colorbar with Manual Ticks
fig.subplots_adjust(right=0.92)
cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
cb = fig.colorbar(c, cax=cbar_ax, label="Net Flux Volume ($cm^3$ / step)\nRed=Uptake, Blue=Lift")

# --- FORCE TICKS AT EXTREMES ---
d_min = norm.vmin
d_max = norm.vmax
ticks = [d_min, d_min/2, 0, d_max/2, d_max]
cb.set_ticks(ticks)
cb.set_ticklabels([f"{t:.4f}" for t in ticks]) 
# -------------------------------

plt.suptitle(f"Root Water Flux Profiles (Full Simulation)", fontsize=14)
plt.savefig(f"{prefix}_flux_diverging_heatmap_FULL.png", dpi=300)
print("Plot saved as 'flux_diverging_heatmap_FULL.png'")
plt.show()


# ============================================================
# PLOT 2: Diverging Flux Heatmaps (LAST DAY ONLY)
# ============================================================

fig2, axes2 = plt.subplots(nrows=1, ncols=4, figsize=(18, 6), sharey=True, sharex=True)

# Define window: Last 1.0 day
t_end = time_arr[-1]
t_start = max(0, t_end - 1.0)

print(f"Generating Last-Day Zoom ({t_start:.2f} - {t_end:.2f} d)...")

for idx, (grp, label) in enumerate(zip(groups, group_labels)):
    ax = axes2[idx]
    
    # Select & Mean per group
    group_data = net_flux_array[:, grp, :]
    mean_flux = np.mean(group_data, axis=1).T 
    
    # Plot using the SAME 'norm' so colors are comparable to the full plot
    c = ax.pcolormesh(time_grid, depth_grid, mean_flux, 
                      cmap='RdBu', norm=norm, shading='auto')
    
    ax.set_title(label, fontsize=10, fontweight='bold')
    ax.set_xlabel("Time [d]")
    if idx == 0: ax.set_ylabel("Depth [cm]")
    if not ax.yaxis_inverted(): ax.invert_yaxis()
    
    # Limit X Axis
    ax.set_xlim(t_start, t_end)

# Colorbar with Manual Ticks (Same logic as above)
fig2.subplots_adjust(right=0.92)
cbar_ax2 = fig2.add_axes([0.93, 0.15, 0.015, 0.7])
cb2 = fig2.colorbar(c, cax=cbar_ax2, label="Net Flux Volume ($cm^3$ / step)\nRed=Uptake, Blue=Lift")

# --- FORCE TICKS AT EXTREMES ---
cb2.set_ticks(ticks)
cb2.set_ticklabels([f"{t:.4f}" for t in ticks]) 
# -------------------------------

plt.suptitle(f"Root Water Flux Profiles (Last 24h Zoom)", fontsize=14)
plt.savefig(f"{prefix}_flux_diverging_heatmap_LAST_DAY.png", dpi=300)
print("Plot saved as 'flux_diverging_heatmap_LAST_DAY.png'")
plt.show()

# ============================================================
# FINAL SUMMARY METRICS (PER GROUP)
# ============================================================
print("\n" + "="*50)
print(" HYDRAULIC LIFT SUMMARY (CUMULATIVE VOLUMES)")
print("="*50)

# Header for a clean table
print(f"{'Group':<25} | {'Uptake (cm3)':<12} | {'Lift (cm3)':<12} | {'Lift/Uptake (%)':<15}")
print("-" * 72)

# Iterate through groups
for grp, label in zip(groups, group_labels):
    # 1. Select data for plants in this group
    # Shape becomes (Time, GroupSize, Layers)
    grp_hl_data = hl_data_array[:, grp, :]
    grp_up_data = uptake_data_array[:, grp, :]
    
    # 2. Sum over all dimensions (Time, Plants in group, Layers)
    total_lift = np.sum(grp_hl_data)
    total_uptake = np.sum(np.abs(grp_up_data)) # abs because uptake is negative
    
    # 3. Calculate Ratio
    ratio = 0.0
    if total_uptake > 0:
        ratio = (total_lift / total_uptake) * 100
    
    # 4. Print Row
    # Shorten label if needed to fit column
    short_label = (label[:22] + '..') if len(label) > 22 else label
    print(f"{short_label:<25} | {total_uptake:12.2f} | {total_lift:12.2f} | {ratio:14.2f} %")

print("-" * 72)
print(f"{'TOTAL SYSTEM':<25} | {np.sum(np.abs(uptake_data_array)):12.2f} | {np.sum(hl_data_array):12.2f} |")
print("="*50 + "\n")