"""Optimize root length and diameter parameters."""

import re
import sys
sys.path.append("../../../")
sys.path.append("../../../src")
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plantbox as pb
from scipy.optimize import differential_evolution


# =============================================================================
# Settings
# =============================================================================

PATH = Path("rootsystem")
BASE_NAME = "RS_base"
TARGET_NAME = "RS_optimized_S_WT"

SEED = 42
TARGET_ERROR = 0.001


# =============================================================================
# Callback
# =============================================================================

def make_callback(name, target=TARGET_ERROR):

    def callback(intermediate_result):
        error = intermediate_result.fun

        print(f"[{name}] Best error: {error:.6f}", flush=True)

        if error <= target:
            print(f"[{name}] Target error reached!", flush=True)
            return True

        return False

    return callback


# =============================================================================
# XML
# =============================================================================

def modify_xml_parameters(input_file, output_file, parameters):

    text = Path(input_file).read_text(encoding="utf-8")

    for section, values in parameters.items():

        pattern = (
            rf'(<([A-Za-z_][\w.-]*)\b'
            rf'(?=[^>]*\bname="{re.escape(section)}")'
            rf'[^>]*>)'
            rf'(.*?)'
            rf'(</\2>)'
        )

        match = re.search(pattern, text, re.DOTALL)

        if match is None:
            raise ValueError(f"Section '{section}' not found.")

        section_text = match.group(3)

        for parameter, value in values.items():

            pattern = (
                rf'(<parameter\b'
                rf'(?=[^>]*\bname="{re.escape(parameter)}")'
                rf'[^>]*\bvalue=")'
                rf'[^"]*(")'
            )

            section_text, n = re.subn(
                pattern,
                rf'\g<1>{value}\g<2>',
                section_text,
            )

            if n == 0:
                raise ValueError(
                    f"Parameter '{parameter}' not found "
                    f"in section '{section}'."
                )

        text = (
            text[:match.start(3)]
            + section_text
            + text[match.end(3):]
        )

    Path(output_file).write_text(text, encoding="utf-8")
    print(f"Saved: {output_file}", flush=True)


# =============================================================================
# Simulation
# =============================================================================

def simulate(xml_file, parameter_function):

    rs = pb.RootSystem()
    rs.readParameters(str(xml_file))

    # Apply optimized parameters.
    parameter_function(rs)

    rs.setSeed(0)
    rs.initializeLB(5, 4)
    
    r, depth = 2.5, 20  # Soil core analysis
    sc_ = pb.SDF_PlantContainer(r, r, depth, False)  # in the center of the root
    sc1 = (pb.SDF_RotateTranslate(sc_,pb.Vector3d(10, 0, 0)))  # shift 10 cm perpendicular to the row
    sc2 = (pb.SDF_RotateTranslate(sc_,pb.Vector3d(10, 0, -20)))  # shift 10 cm perpendicular to the row, 20cm downwards 
    sc3 = (pb.SDF_RotateTranslate(sc_,pb.Vector3d(10, 0, -40)))  # shift 10 cm perpendicular to the row, 40cm downwards 
    sc = [sc1, sc2, sc3]
    scVol = depth * r * r * np.pi
    interrow = 45 #cm
    interplant = 20 #cm 

    simulated_RLD = np.zeros((3,len(times)))
    simulated_diameter = np.zeros((3,len(times)))

    for day in range(1, times[-1] + 1):

        rs.simulate(1, True)

        if day not in times:
            continue

        i = np.where(times == day)[0][0]
        
        for j in range(0, len(sc)): 
            ana = pb.SegmentAnalyser(rs)
            ana.mapPeriodic(interrow, interplant)
            ana.crop(sc[j])
            
            length = np.asarray(ana.getParameter("length"))
            simulated_RLD[j,i] = np.sum(length)/scVol

            radius = np.asarray(ana.getParameter("radius"))
            simulated_diameter[j,i] = (np.sum(length * radius) / (np.sum(length)+1e-14) * 20)

    return simulated_RLD, simulated_diameter


def rmse(simulated, measured):
    return np.sqrt(np.nanmean(((simulated - measured) / measured) ** 2))


# =============================================================================
# Objective functions
# =============================================================================

def objective_length(params):

    def set_parameters(rs):

        p = rs.getRootRandomParameter()

        p[2].ln   = params[0]
        p[2].r    = params[1]

        p[4].ln = params[2]
        p[4].r  = params[3]
        p[4].tropismN = params[4]
        p[4].tropismS  = params[5]

        p[5].r = params[6]
        p[5].tropismN = params[7]
        p[5].tropismS  = params[8]

    simulated, _ = simulate(
        PATH / f"{BASE_NAME}.xml",
        set_parameters,
    )

    return rmse(simulated, real_RLD)


def objective_diameter(params):

    def set_parameters(rs):

        p = rs.getRootRandomParameter()

        p[1].a = params[0]
        p[2].a = params[1]
        p[3].a = params[2]
        p[4].a = params[3]
        p[5].a = params[4]

    _, simulated = simulate(
        PATH / f"{TARGET_NAME}.xml",
        set_parameters,
    )

    return rmse(simulated, real_diam)


# =============================================================================
# Experimental data
# =============================================================================

df = pd.read_csv("data/RLD_RD_field_experiment.csv")

data = df[
    (df["substrate"] == "S")
    & (df["genotype"] == "WT")
]

DAS = data["DAS"].values
times = np.unique(DAS[:3])

depth = np.unique(data["depth"].values)

real_RLD = (
    data
    .pivot(index="DAS", columns="depth", values="RLD_mean")
    .reindex(index=times, columns=depth)
    .fillna(np.nan)
    .to_numpy()
).T


real_RLD_SE = (
    data
    .pivot(index="DAS", columns="depth", values="RLD_SE")
    .reindex(index=times, columns=depth)
    .fillna(np.nan)
    .to_numpy()
).T

real_diam = (
    data
    .pivot(index="DAS", columns="depth", values="RD_mean")
    .reindex(index=times, columns=depth)
    .fillna(np.nan)
    .to_numpy()
).T*1e-3    #microm --> mm


real_diam_SE = (
    data
    .pivot(index="DAS", columns="depth", values="RD_SE")
    .reindex(index=times, columns=depth)
    .fillna(np.nan)
    .to_numpy()
).T*1e-3    #microm --> mm


# =============================================================================
# 1. Optimize root length of soil cores in different depths 
# =============================================================================

length_bounds = [
    (0.5, 2),    # lateral ln
    (0.1, 4),    # lateral r
    (0.05, 0.5), # basal ln
    (1.5, 4),    # basal r
    (0.05, 1),   # basal tropismN
    (0.05, 0.2),   # basal tropismS
    (1.5, 4),    # shoot-borne r
    (0.05, 1),   # shoot-borne tropismN
    (0.05, 0.2),   # shoot-borne tropismS
]

print("\nOptimizing root length density...", flush=True)

length_result = differential_evolution(
    objective_length,
    bounds=length_bounds,
    seed=SEED,
    maxiter=100,
    callback=make_callback("Length"),
    polish=False,
)

print("Root length density optimization finished.", flush=True)


# Save optimized length parameters.
length_parameters = {
    "lateral": {
        "ln": length_result.x[0],
        "r": length_result.x[1],                  
    },
    "basal": {
        "ln": length_result.x[2],
        "r": length_result.x[3],
        "tropismN": length_result.x[4],
        "tropismS": length_result.x[5],
    },
    "shootborne": {
        "r": length_result.x[6],
        "tropismN": length_result.x[7],
        "tropismS": length_result.x[8],
    },
}

modify_xml_parameters(
    PATH / f"{BASE_NAME}.xml",
    PATH / f"{TARGET_NAME}.xml",
    length_parameters,
)


# =============================================================================
# 2. Optimize root diameter of soil cores in different depths 
# =============================================================================

diameter_bounds = [
    (0.03, 0.05),    # taproot a
    (0.001, 0.03),  # lateral a
    (0.001, 0.01),  # tertiary a
    (0.03, 0.05),    # basal a
    (0.03, 0.15),    # shoot-borne a
]

print("\nOptimizing root diameter...", flush=True)

diameter_result = differential_evolution(
    objective_diameter,
    bounds=diameter_bounds,
    seed=SEED,
    maxiter=100,
    callback=make_callback("Diameter"),
    polish=False,
)

print("Root diameter optimization finished.", flush=True)


# Save optimized diameter parameters.
diameter_parameters = {
    "taproot": {"a": diameter_result.x[0]},
    "lateral": {"a": diameter_result.x[1]},
    "tertiaryroots": {"a": diameter_result.x[2]},
    "basal": {"a": diameter_result.x[3]},
    "shootborne": {"a": diameter_result.x[4]},
}


modify_xml_parameters(
    PATH / f"{TARGET_NAME}.xml",
    PATH / f"{TARGET_NAME}.xml",
    diameter_parameters,
)


# =============================================================================
# Results
# =============================================================================

print("\n" + "=" * 50)
print("OPTIMIZATION RESULTS")
print("=" * 50)
print(TARGET_NAME)
for name, result in [
    ("Length", length_result),
    ("Diameter", diameter_result),
]:
    print(f"\n{name}:")
    print(f"  Success:        {result.success}")
    print(f"  Iterations:     {result.nit}")
    print(f"  Evaluations:    {result.nfev}")
    print(f"  Best error:     {result.fun:.6f}")
    print(f"  Message:        {result.message}")

