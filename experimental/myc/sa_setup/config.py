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
    param: str              # name attribute of the <parameter> element, or "v_cv" (dev of v = CV × v)
    low: float
    high: float
    log: bool = False       # True -> log-uniform between low and high
    source: str | None = None   # justification of the range (literature, data, assumption)
    role: str = "factor"        # "factor" or "control" (should have no influence)


def _pm(rel, nominal):
    """Placeholder range nominal*(1±rel)."""
    return nominal * (1 - rel), nominal * (1 + rel)


PH = 0.25  # placeholder span ±25 %
V_SOURCE = "240 µm/h (Oyarte Galvez et al. 2025), range 180–300 µm/h around the centre of the measured tip speeds"

FACTORS: list[Factor] = [
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

# Parameters that are NEVER varied (documentation / checks only)
FIXED = ["successor percentage (0.75 BAS, 0.25 Runner)", "la", "lb", "BAS lmax (≈ 3 × dx)", "a",
         "tropismT", "tropismN", "tropismS", "dx", "dxMin", "ldelay", "b_prob"]

# --------------------------------------------------------------------------
# Simulation (as petridishsimauto.py)
# --------------------------------------------------------------------------
N_RINGS = 15
HYPHAE_LAYER = 0.1          # hyphae grow on the cellophane, thin layer on top [cm]
ROOT_DAYS = 5               # root phase [d]
STEPS_PER_DAY = 24          # dt = 1 h
HOURS_HYPHAE = 60           # simulation after t* [h], as the figures in the paper
MAX_WAIT_DAYS = 30          # without t* after this many days: the run counts as failed
# No abort on a stalling front: in the SA that is a result (small v), not an error
MIN_GROWING_TIPS = None
STALL_FRACTION = None
SAVE_VTP = False            # save the final geometry as .vtp (for network evaluations, large)
ROOT_EFFICIENCY = True      # relative root efficiency at the end (AMFNetwork.relativeRootEfficiency, 2D) in raw.npz

# --------------------------------------------------------------------------
# Outputs (evaluate.py)
# --------------------------------------------------------------------------
SNAPSHOT_HOURS = [10, 30, 60]   # evaluation of quantities after t* [h], as in the paper (0, 10, 30, 60 h)
# Test with a longer simulation (python design.py long): LONG_SEEDS runs at the nominal point over LONG_HOURS h,
# seeds as p0_noise_000.. Additionally evaluated at EXTRA_HOURS, but only in runs simulated that long
LONG_HOURS = 120
LONG_SEEDS = 10
EXTRA_HOURS = [90, 120]
# Quantities at every time in SNAPSHOT_HOURS
SNAPSHOT_OUTPUTS = ["total_length", "active_tips", "bas_share", "anastomoses", "front_position",
                    "spread_r90", "max_tip_density"]
# Quantities over the whole window up to max(SNAPSHOT_HOURS), stored at the last time.
# Only quantities without a fit: the wave fits (front_velocity, K1, K2, wave_width) are defined in Phase 0 at 60 h in only
# 12–33 of 47 runs, too rarely for Morris (all 14 points of a trajectory). They are still computed (diagnostics).
# beta: anastomosis constant from a = beta rho n [µm/h], calibration target ≈ 23 µm/h (Oyarte Galvez et al. 2025)
WAVE_OUTPUTS = ["front_speed", "encounter_C", "beta"]
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
# dt convergence: time steps per day (24 = 1 h as the image interval of the experiment, which only fixes the evaluation
# interval). Anastomosis is checked once per step at the tip node, ana and the BAS tip splitting (b dt) are per step,
# so these outputs may depend on dt. Runtime grows roughly with the number of steps
PHASE0_DT_STEPS_PER_DAY = [24, 48, 96]   # dt = 1 h, 0.5 h, 0.25 h
PHASE0_DT_SEEDS = 10
PHASE0_RANGE_POINTS = 20    # random points in the factor space, one run each: abort rate and runtime at the edges

# --------------------------------------------------------------------------
# Phase 1 (Morris)
# --------------------------------------------------------------------------
MORRIS_LEVELS = 4            # p
MORRIS_TRAJECTORIES = 20     # r (selected from MORRIS_CANDIDATES by spread)
MORRIS_CANDIDATES = 300
MORRIS_REPLICATES = 5        # Phase 0: seed CV 0.3–0.97 per output -> 5 instead of 3
BASE_SEED = 20261021
