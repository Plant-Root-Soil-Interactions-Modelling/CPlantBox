"""
Central configuration of the sensitivity analysis (Petri dish setup, TwoHyphaePlusBAS.xml).

Defined here:
  - which parameters (factors) are varied and over which range,
  - how the `dev` values are handled,
  - seeds, replicates, simulation time after t*, evaluation times and outputs.

Ranges with `source=None` are placeholders (±25 % around the nominal value) and are reported
as open by `check_config.py`.
"""
import math
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
MYC_DIR = ROOT.parent                              # experimental/myc (AMF*.py)
CPLANTBOX_DIR = ROOT.parents[2]                    # for `import plantbox`
TEMPLATE_XML = MYC_DIR / "tomatoparameters" / "TwoHyphaePlusBAS.xml"   # nominal model
RUNS_DIR = ROOT / "runs"          # one folder per run with XML, raw output, meta.json
RESULTS_DIR = ROOT / "results"    # designs, output tables, analyses, plots

# --------------------------------------------------------------------------
# Handling of `dev`
# --------------------------------------------------------------------------
# "relative": dev is scaled with the mean (CV = dev/value stays constant).
# "absolute": dev stays at the nominal value from the template XML.
# The factor `v_cv` sets dev of v to CV × v regardless.
DEV_MODE = "relative"

# --------------------------------------------------------------------------
# Factors
# --------------------------------------------------------------------------
RUNNER, BAS = 1, 2   # subTypes in TwoHyphaePlusBAS.xml
SUBTYPE_NAMES = {RUNNER: "Runner", BAS: "BAS"}


@dataclass
class Factor:
    name: str               # unique name in design/analysis, e.g. "R_v"
    subtypes: tuple         # subTypes of the <hyphae> the value is written to (several = shared factor)
    # name attribute of the <parameter> element, or a special one (xml_tools.write_xml): "v_cv" (dev of v = CV × v),
    # "successor" (share of the first successor), "ana_L_scale" (anaMinDist = factor × derived), "cv_scale" (all dev but v × factor)
    param: str
    low: float              # range: Morris range in the SA = search range in the calibration
    high: float
    log: bool = False       # True -> log-uniform between low and high
    source: str | None = None   # justification of the range (literature, data, assumption)
    role: str = "factor"        # "factor" or "control" (should have no influence)
    # where the value goes: "hyphae" (<parameter> of the <hyphae> subtypes), "root_scale" (factor on the <parameter> `param`
    # of all <root> elements), "sim" (keyword argument `param` of model_hook.simulate / makeSimulation)
    target: str = "hyphae"
    nominal: float | None = None    # nominal value for special params, "root_scale" and "sim" (else: from the template XML)
    also: tuple = ()                # further <parameter>s of the same subtypes that get the same value (lb = ln)
    prior: tuple | None = None      # calibration: (lo, hi) literature range as ±2σ (log parameters: in log); None = weak prior
    calibrate: bool = False         # calibration: free parameter (True) or fixed at the nominal value (False)


def _pm(rel, nominal):
    """Placeholder range nominal*(1±rel)."""
    return nominal * (1 - rel), nominal * (1 + rel)


PH = 0.25  # placeholder span ±25 %
V_SOURCE = "240 µm/h (Oyarte Galvez et al. 2025), range 180–300 µm/h around the centre of the measured tip speeds"

deg = math.radians
THETA_SRC = ("emergence angle of the model (3D); the measured angles (Oyarte Galvez 2025, Ext. Data Fig. 11b: RH 73°, "
             "BAS 80°, projected, from the initial growth direction) are calibration targets, so only a weak prior; range 40–100°")

# --------------------------------------------------------------------------
# PARAMETERS: ONE list for SA and calibration (B13). Range = Morris range = calibration search range.
# calibrate=True: free in the calibration, False: fixed at the nominal value there (still varied in the SA).
# Choose by hand after the Morris screening (analyze_morris.py writes a suggestion, morris_calibration_suggestion.csv).
# --------------------------------------------------------------------------
PARAMETERS: list[Factor] = [
    # ---- RunnerHyphae (subType 1, lateral branching, b_prob = 0) ----
    Factor("R_v", (RUNNER,), "v", 0.432, 0.72, prior=(0.432, 0.72), calibrate=True, source=V_SOURCE),
    Factor("R_hlt", (RUNNER,), "hlt", 4.0, 12.5, prior=(4.0, 12.5), calibrate=True,
           source="SA range Schnepf 2016 4–10 d, death rates Schnepf 2008 4–12.5 d; range = prior"),
    Factor("R_lmax", (RUNNER,), "lmax", 1.5, 5.0, calibrate=True,
           source="not in the literature; estimate v × lifetime (≈ 1.3–5.4 cm)"),
    Factor("R_ln", (RUNNER,), "ln", 0.05, 0.6, log=True, also=("lb",), prior=(0.1, 0.6), calibrate=True,
           source="~1 mm (Giovannetti, in Schnepf 2008) to ~6 mm new growing tips (Oyarte Galvez 2025); lb coupled to ln"),
    Factor("R_theta", (RUNNER,), "theta", deg(40), deg(100), calibrate=True, source=THETA_SRC),
    Factor("R_tropismN", (RUNNER,), "tropismN", 1.0, 8.0, calibrate=True,
           source="estimate; the persistence length of the runners (≥ 10 cm, Ext. Data Fig. 11a) is a calibration target"),
    Factor("R_succ_BAS", (RUNNER,), "successor", 0.5, 0.95, calibrate=True,
           source="estimate: share of BAS among the runner laterals, nominal 0.75"),
    Factor("R_a", (RUNNER,), "a", 0.001, 0.01, log=True, prior=(0.0025, 0.01), calibrate=True,
           source="radius 50 µm (Schnepf et al. 2017), prior factor 2; lower bound 10 µm (Bago et al. 1998: diameter > 20 µm)"),
    # ---- BAS (subType 2, tip splitting, b_prob = 1) ----
    Factor("B_v", (BAS,), "v", 0.432, 0.72, prior=(0.432, 0.72), calibrate=True, source=V_SOURCE),
    # B_b without literature prior: the measurement (Oyarte Galvez 2025) is already the target beta, it would count twice
    Factor("B_b", (BAS,), "b", 0.3, 2.0, log=True, calibrate=True,
           source="tip splitting rate; branching rate ≈ 0.044 h⁻¹ ≈ 1.05 /d (Oyarte Galvez et al. 2025) is the target beta"),
    Factor("B_hlt", (BAS,), "hlt", 3.0, 12.0, prior=(5.0, 9.0), calibrate=True,
           source="BAS degenerate after ~7 d (Bago 1998; Kokkoris 2026)"),
    Factor("B_theta", (BAS,), "theta", deg(40), deg(100), calibrate=True, source=THETA_SRC),
    Factor("B_tropismN", (BAS,), "tropismN", 0.5, 4.0, calibrate=True, source="estimate"),
    Factor("B_a", (BAS,), "a", 0.00005, 0.005, log=True, prior=(0.0001, 0.0005), calibrate=True,
           source="diameter 2–10 µm = radius 1–5 µm (Bago et al. 1998, Mycologia, doi 10.1080/00275514.1998.12026878)"),
    # a BAS segment splits only while it grows (until lmax): b lmax / v sets the size of the BAS (fewer than 1 child, intended)
    Factor("B_lmax", (BAS,), "lmax", 0.01, 0.1, log=True, calibrate=True,
           source="no literature value, nominal 0.03 cm; lower bound = dx (one segment)"),
    # ---- shared by runner and BAS ----
    Factor("ana", (RUNNER, BAS), "ana", 0.3, 1.0, calibrate=True,
           source="unknown probability of fusion at an encounter (it does not always happen); nominal 0.7, inside the range"),
    # capture distance tip -> hypha: larger than the hyphal diameter, tips appear to sense other hyphae; ~100 µm plausible,
    # possibly more (Richter et al. 2024). Prior 50–300 µm as ±2σ (log, centre ≈ 120 µm)
    Factor("distTH", (RUNNER, BAS), "distTH", 0.003, 0.05, log=True, prior=(0.005, 0.03), calibrate=True,
           source="tips sense other hyphae beyond contact; ~100 µm plausible, possibly more (Richter et al. 2024)"),
    # minimal path length for anastomosis: factor on the derived value distTH / sin(theta) + 2 dx (A11). Nominal 3
    # (ANA_L_SCALE_NOMINAL): with the derived value (×1) young laterals fuse with their own parent/child near the branching
    # point (they bend back parallel to it, effective angle ~6° instead of theta) - 38 % of the anastomoses at dt = 1 h,
    # 91 % at 0.25 h. From ×3 (≈ 0.095 cm nominal) the anastomosis rate per tip no longer depends on dt (test 2026-11, 8 seeds)
    Factor("ana_L_scale", (RUNNER, BAS), "ana_L_scale", 1.5, 6.0, log=True, nominal=3.0, prior=(1.5, 6.0),
           calibrate=True, source="dt test: ×3 the derived value suppresses fusions with the own parent/child; ×/÷ 2"),
    # CV of the tip speed (dev of v = CV × v, both sub types): the front is carried by the fast tips, so v_cv together with
    # v sets the speed at the front (≈ 300 µm/h, target front_tip_speed) and the wave speed (≈ 280 µm/h) above the mean 240
    Factor("v_cv", (RUNNER, BAS), "v_cv", 0.15, 0.40, prior=(0.2, 0.4), calibrate=True,
           source="spread of the measured tip speeds (Oyarte Galvez et al. 2025), nominal CV ≈ 0.29; prior 0.2–0.4"),
]

# Only in the SA (no target acts on them, or a check)
SA_ONLY: list[Factor] = [
    # all dev except v together (A9): the CVs are estimates without ranges
    Factor("cv_scale", (RUNNER, BAS), "cv_scale", 0.5, 2.0, log=True, nominal=1.0,
           source="the dev of the parameters are estimates; all of them (except v) ×/÷ 2"),
    # control factor: b is only read during tip splitting (Hyphae.cpp), no effect with b_prob = 0 -> mu* = 0
    Factor("R_b_ctrl", (RUNNER,), "b", 0.7125, 1.1875, role="control"),
]

FACTORS: list[Factor] = PARAMETERS + SA_ONLY
FACTOR_SET = "v2"   # stored in every new design (column factor_set); designs without it use FACTORS_V1

# --------------------------------------------------------------------------
# FACTORS_V1: the factors of the first Morris round (13), kept so that its runs can still be evaluated and resubmitted
# --------------------------------------------------------------------------
FACTORS_V1: list[Factor] = [
    # ---- RunnerHyphae (subType 1, lateral branching, b_prob = 0) ----
    Factor("R_v",      (RUNNER,), "v",     0.432, 0.72, source=V_SOURCE),
    Factor("R_hlt",    (RUNNER,), "hlt",   4.0, 10.0,
           source="sensitivity range Schnepf 2016; death rates Schnepf 2008 give 4–12.5 d"),
    Factor("R_lmax",   (RUNNER,), "lmax",  1.5, 5.0,
           source="not in the literature; estimate v × lifetime (≈ 1.3–5.4 cm)"),
    Factor("R_ln",     (RUNNER,), "ln",    0.1, 0.6, log=True,
           source="~1 mm all branches (Giovannetti, in Schnepf 2008) to ~6 mm new growing tips "
                  "(Oyarte Galvez 2025); log-uniform, since factor 6 and lateral rate ∝ 1/ln"),
    Factor("R_theta",  (RUNNER,), "theta", math.radians(50), math.radians(80),
           source="60° (Schnepf 2016) to ≈ 70° (Oyarte Galvez 2025), angle definitions differ"),
    # control factor: b is only read during tip splitting (Hyphae.cpp), no effect with b_prob = 0 -> mu* = 0
    Factor("R_b_ctrl", (RUNNER,), "b",     *_pm(PH, 0.95), role="control"),

    # ---- BAS (subType 2, tip splitting, b_prob = 1) ----
    Factor("B_v",      (BAS,), "v",     0.432, 0.72, source=V_SOURCE),
    Factor("B_b",      (BAS,), "b",     0.55, 1.35,
           source="branching rate / tip density (Oyarte Galvez et al. 2025, Fig. b): ≈ 0.044 h⁻¹ ≈ 1.05 /d at the maximum, "
                  "≈ 0.65 /d at the flank; band 2 SEM ≈ ±40 %. Counts all branchings there (including runner laterals)"),
    Factor("B_hlt",    (BAS,), "hlt",   5.0, 9.0,
           source="BAS degenerate after ~7 d (Bago 1998; Kokkoris 2026)"),
    Factor("B_theta",  (BAS,), "theta", math.radians(50), math.radians(80),
           source="as R_theta; measured in Oyarte Galvez 2025 with 881 BAS branching points"),

    # ---- shared by runner and BAS ----
    Factor("distTH", (RUNNER, BAS), "distTH", *_pm(PH, 0.01)),   # nominal 0.01 cm = dx
    Factor("ana",    (RUNNER, BAS), "ana",    0.5, 1.0),        # probability, nominal 1 = bound
    Factor("v_cv",   (RUNNER, BAS), "v_cv",   0.15, 0.40,
           source="spread of the measured tip speeds (Oyarte Galvez et al. 2025), nominal CV ≈ 0.29"),
]

# Inoculum / root factors (python design.py inoculum): root phase, root architecture and placement have no source or data
# (tomato as a stand-in for the root organ culture of D. carota). Test whether the normalisation to t* removes their effect:
# they may change wait_days and the supply of hyphae, but should hardly change the outputs after t*. Hyphal parameters nominal.
# The emergence density (Schnepf et al.) is included for comparison with a documented parameter.
INOCULUM_FACTORS: list[Factor] = [
    Factor("root_days", (), "root_days", 3.0, 10.0, target="sim", nominal=5.0,
           source="assumption: the roots should grow a bit before the hyphae emerge"),
    Factor("root_r_scale", (), "r", 0.5, 2.0, log=True, target="root_scale", nominal=1.0,
           source="no data on the root organ culture: growth rate of all roots ×/÷ 2"),
    Factor("root_ln_scale", (), "ln", 0.5, 2.0, log=True, target="root_scale", nominal=1.0,
           source="no data on the root organ culture: branching distance of all roots ×/÷ 2 (ROC strongly branched)"),
    Factor("seed_pos_x", (), "seed_pos_x", -4.2, -1.2, target="sim", nominal=-4.2,
           source="placement of the root piece: 5 mm from the edge (nominal) to ~1 cm before the barrier"),
    Factor("emergence_density", (), "hyphal_emergence_density", 2.0, 6.0, target="sim", nominal=4.0,
           source="Schnepf et al. 4 /cm, ±50 % for comparison"),
]
INOCULUM_TRAJECTORIES = 10
INOCULUM_REPLICATES = 3
# in addition to the SA outputs; frozen_share: tips frozen by the death of an ancestor (A6), grows with the waiting time
INOCULUM_OUTPUTS = ["wait_days", "influx_tstar", "influx_window", "influx_decay", "frozen_share"]
# state of the colony at t* (D23), SA outputs in every Morris: the hyphal parameters act already before t*; with them one
# sees which parameters act mainly via the state at t* (waiting time, entry points, supply, frozen tips)
STATE_OUTPUTS = ["wait_days", "n_crossings_tstar", "influx_tstar", "frozen_share"]


def factor_set(name):
    """Factor list of a design: "v2" (FACTORS, default of new designs) or "v1" (first Morris round)."""
    return {"v1": FACTORS_V1, "v2": FACTORS}[name]


# Parameters that are NEVER varied (documentation / checks only)
FIXED = ["la", "lb of the BAS", "tropismT", "tropismS", "dx", "dxMin", "ldelay", "b_prob"]

# --------------------------------------------------------------------------
# Simulation (as petridishsimauto.py)
# --------------------------------------------------------------------------
N_RINGS = 15
# rings as in Oyarte Galvez et al. 2025: densities against the distance r from the barrier crossing point (per replicate),
# 15 rings of equal area; their positions in Fig. 2c (r ≈ 11 … 31 mm) fit 15 rings up to ≈ 32 mm (derived, not stated)
RING_RADIUS = 3.2           # [cm] outer radius of the last ring (former: the dish radius 4.7)
RING_CENTER = "crossing"    # "crossing": median of the barrier crossings until t* of each run; "dish": centre of the dish
# warnings (analyze_phase0.py -> results/phase0_warnings.txt, calibration/report.txt)
WARN_MIN_CROSSINGS = 5          # fewer barrier crossings until t* (median over runs): the crossing point is uncertain
WARN_MAX_OUTSIDE_RINGS = 0.10   # larger share of the growing tips beyond RING_RADIUS at the end: the rings miss the front
HYPHAE_LAYER = 0.1          # hyphae grow on the cellophane, thin layer on top [cm]
ROOT_DAYS = 5               # root phase [d]: assumption without source (the roots should grow a bit before the hyphae), see INOCULUM_FACTORS
# dt = 1 h (image interval of the experiment) is part of the model definition: the parameters hold for dt = 1 h.
# Isolated dt tests (2026-11, 8 seeds, ana_L_scale 3): anastomoses per tip independent of dt, but the network after 2 d is
# ~15 % smaller at dt = 0.25 h without anastomosis (tip growth per time −6 %, branching follows) and ~25 % with it
STEPS_PER_DAY = 24
# simulation after t* [h]: the wave fits need the longer window. Long test (Phase 0): at 60 h K1 763 / K2 0.071, at 120 h
# K1 1422 / K2 0.151 (targets 1400 / 0.15); Fig. 2 of the paper spans about 3 days. Wave quantities over 0..120 h
HOURS_HYPHAE = 120
MAX_WAIT_DAYS = 30          # without t* after this many days: the run counts as failed
# as the experiment: plates with no more than 1 cm hyphal length in the hyphal compartment were excluded (unhealthy AMF).
# Runs without t* with at most this length on the fungal side get the status "excluded", the others "failed" (slow)
EXCLUDE_CROSSED_CM = 1.0
# No abort on a stalling front: in the SA that is a result (small v), not an error
MIN_GROWING_TIPS = None
STALL_FRACTION = None
SAVE_VTP = False            # save the final geometry as .vtp (for network evaluations, large)
ROOT_EFFICIENCY = True      # relative root efficiency at the end (AMFNetwork.relativeRootEfficiency, 2D) in raw.npz

# --------------------------------------------------------------------------
# Outputs (evaluate.py)
# --------------------------------------------------------------------------
# evaluation of quantities after t* [h] (paper: 0, 10, 30, 60 h; 20, 40 for the time courses; up to HOURS_HYPHAE).
# The LAST one is the window of the wave quantities and of the pooled targets (= HOURS_HYPHAE)
SNAPSHOT_HOURS = [10, 20, 30, 40, 60, 90, 120]
# times after t* [h] at which the quantities of the image graph (lengths, BAS share, root efficiency) are computed (model_hook)
# and available as snapshot outputs (IMAGE_OUTPUTS); the calibration uses their time courses (figure with L_RH, ρ_BAS, Ê_r)
IMAGE_HOURS = [20, 40, 60]
# tip speed at the front (model_hook.tip_speed_diagnostics, Fig. 2d): growing runner tips whose distance from the centre of
# the rings is at least this quantile of all growing tips at that time count as "at the front"; averaged from SPEED_FROM h on
FRONT_QUANTILE = 0.9
SPEED_FROM = 10
# branching angle (model_hook.branching_angle_diagnostics, Extended Data Fig. 11b): directions over these arc lengths [cm]
# (the paper does not say how the angle was measured: "from the initial growth direction")
ANGLE_L_PARENT = 0.01       # parent: tangent over ±100 µm around the branching point
ANGLE_L_CHILD = 0.02        # child: first 200 µm
# persistence length (model_hook.persistence_diagnostics, Extended Data Fig. 11a): lags up to 1 cm as in the figure
PERSISTENCE_MAX_LAG = 1.0   # [cm]
PERSISTENCE_MIN_LENGTH = 0.3   # [cm] runners shorter than this are left out
# Test with a longer simulation (python design.py long): LONG_SEEDS runs at the nominal point over LONG_HOURS h,
# seeds as p0_noise_000.. Additionally evaluated at EXTRA_HOURS, but only in runs simulated that long
LONG_HOURS = 180           # beyond HOURS_HYPHAE: are the wave quantities converged at 120 h?
LONG_SEEDS = 10
EXTRA_HOURS = [150, 180]
# Quantities at every time in SNAPSHOT_HOURS
SNAPSHOT_OUTPUTS = ["total_length", "active_tips", "bas_share", "anastomoses", "front_position",
                    "spread_r90", "max_tip_density"]
# Quantities over the whole window up to max(SNAPSHOT_HOURS), stored at the last time (definitions: evaluate.WAVE_QUANTITIES).
# Per run in the SA: only quantities without a fit. The wave fits are defined per run in Phase 0 at 60 h in only 12–33 of 47
# runs, too rarely for Morris (all points of a trajectory); per run they are still computed (diagnostics).
WAVE_OUTPUTS = ["front_speed", "encounter_C", "front_tip_speed", "tip_speed_mean", "branch_angle_rh", "branch_angle_bas",
                "persistence_length"]
# Pooled per design point over its replicates, as in the calibration (evaluate_pooled.py): all calibration targets that are
# wave quantities (bas_share is a snapshot output). alpha = anastomosis, beta = branching (Oyarte Galvez et al. 2025: the
# other way round). The noise level comes from a bootstrap over the replicates
POOLED_OUTPUTS = ["alpha", "beta", "K1", "K2", "front_velocity", "wave_width"]
# Image graph at IMAGE_HOURS (evaluate.IMAGE_QUANTITIES), per run like the snapshot outputs: calibration targets with their
# time course (root_efficiency_rel, bas_share_image, L_RH_image) plus comparison quantities
IMAGE_OUTPUTS = ["root_efficiency_rel", "bas_share_image", "L_RH_image", "L_BAS_image", "bas_share_subtype_image",
                 "ana_per_length_rh", "ana_per_length_total", "global_efficiency_rel", "cost_rel"]
# held back from the calibration for the validation of the calibrated model (Oyarte Galvez et al. 2025, figure e:
# Ê_g ≈ 0.33, 0.44, 0.53 and Ĉ ≈ 0.1 at 20, 40, 60 h after t*); normalisations assumed: (X - X_MST) / (X_DT - X_MST)
VALIDATION_OUTPUTS = {"global_efficiency_rel": {20: 0.33, 40: 0.44, 60: 0.53}, "cost_rel": {20: 0.10, 40: 0.10, 60: 0.10}}
GLOBAL_EFF_MAX_SOURCES = 300   # random subset of the terminal and root nodes for Ê_g (all pairs would be slow)
POOL_MIN_RUNS = 3       # fewer finished replicates at a design point: pooled outputs undefined there
POOL_BOOT = 100
# Quality criteria of the wave fits as petridishstats.py
R2_MIN_LENGTH = 0.9
R2_MIN_TIPS = 0.5
SATURATION_MIN = 0.8
MIN_RINGS = 3
ENCOUNTER_MIN_TIPS = 5

# --------------------------------------------------------------------------
# Phase 0
# --------------------------------------------------------------------------
PHASE0_NOISE_SEEDS = 50
PHASE0_DX_VALUES = [0.01, 0.005, 0.0025]   # nominal, half, quarter (runtime grows roughly with 1/dx!)
PHASE0_DX_SEEDS = 10
# dt convergence: time steps per day (24 = 1 h, fixed as part of the model, see STEPS_PER_DAY). Documents the remaining
# dt dependence at the nominal point. Runtime grows roughly with the number of steps
PHASE0_DT_STEPS_PER_DAY = [24, 48, 96]   # dt = 1 h, 0.5 h, 0.25 h
PHASE0_DT_SEEDS = 10
# Thickness of the hyphal layer on the fungal side: the cellophane of the experiment is quasi 2D, the model uses a 1 mm
# layer in 3D. With distTH = 100 µm the layer thickness changes the encounter probability (anastomoses, beta, encounter_C).
# Compared with the thinnest layer (closest to 2D), same seeds. Thin layers: tips hit the boundary more often
PHASE0_LAYER_VALUES = [0.1, 0.05, 0.02]   # [cm], nominal HYPHAE_LAYER first
PHASE0_LAYER_SEEDS = 10
# Minimal path length along the mycelium for anastomosis (anaMinDist) as a factor on the derived value
# distTH / sin(theta) + 2 dx (≈ 0.032 cm nominal). Nominal ANA_L_SCALE_NOMINAL = 3 (≈ 0.095 cm), compared with it, same
# seeds: 1 = the derived value (fusions with the own parent/child), 1.5 and 6 = the ends of the factor range
ANA_L_SCALE_NOMINAL = 3.0
PHASE0_ANAL_SCALES = [ANA_L_SCALE_NOMINAL, 1.0, 1.5, 6.0]
PHASE0_ANAL_SEEDS = 10
# Common random numbers (C19): how much of the seed noise cancels in the difference between two nearby points with the
# SAME seeds (as the calibration differences and the Morris steps) compared with INDEPENDENT seeds? One factor is moved
# by PHASE0_CRN_STEP in the normalised space (the calibration step DIFF_STEP); the reference are the noise runs
PHASE0_CRN_PARAM = "R_v"
PHASE0_CRN_STEP = 0.15
PHASE0_CRN_SEEDS = 10
PHASE0_RANGE_POINTS = 20    # random points in the factor space, one run each: abort rate and runtime at the edges

# --------------------------------------------------------------------------
# Phase 1 (Morris)
# --------------------------------------------------------------------------
MORRIS_LEVELS = 4            # p
MORRIS_TRAJECTORIES = 20     # r (selected from MORRIS_CANDIDATES by spread)
MORRIS_CANDIDATES = 300
MORRIS_REPLICATES = 5        # Phase 0: seed CV 0.3–0.97 per output -> 5 instead of 3
BASE_SEED = 20261021
