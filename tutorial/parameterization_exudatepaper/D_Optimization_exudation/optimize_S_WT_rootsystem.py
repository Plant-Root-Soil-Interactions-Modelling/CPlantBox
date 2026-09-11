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
BASE_NAME = "RS_base_sand"
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
    tube_ = pb.SDF_PlantContainer(10, 10, 60, False) #tube with diameter 20, length 60 cm 
    tube = pb.SDF_RotateTranslate(tube_, 0, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 10))
    rs.setGeometry(tube)
    rs.initializeLB(5, 4)

    simulated_length = np.zeros(len(times))
    simulated_diameter = np.zeros(len(times))

    for day in range(1, times[-1] + 1):

        rs.simulate(1, True)

        if day not in times:
            continue

        i = np.where(times == day)[0][0]

        length = np.asarray(rs.getParameter("length"))
        radius = np.asarray(rs.getParameter("radius"))

        simulated_length[i] = np.sum(length)
        simulated_diameter[i] = (
            np.sum(length * radius) / np.sum(length) * 20
        )

    return simulated_length, simulated_diameter


def rmse(simulated, measured):
    return np.sqrt(np.mean(((simulated - measured) / measured) ** 2))


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

        p[5].r = params[4]

    simulated, _ = simulate(
        PATH / f"{BASE_NAME}.xml",
        set_parameters,
    )

    return rmse(simulated, real_length)


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

df = pd.read_csv("data/column_experiment_mean.csv")

data = df[
    (df["substrate"] == "S")
    & (df["genotype"] == "WT")
]

DAS = data["DAS"].values

times = DAS[:2]

real_length = data["RL_mean"].values[:len(times)]
real_length_SE = data["RL_SE"].values[:len(times)]

real_diam = data["RD_mean"].values[:len(times)]
real_diam_SE = data["RD_SE"].values[:len(times)]


# =============================================================================
# 1. Optimize root length
# =============================================================================

length_bounds = [                              
    (0.5, 2),    # lateral ln
    (0.1, 4),    # lateral r 
    (0.05, 0.5), # basal ln
    (1.5, 4),    # basal r
    (1.5, 4),    # shoot-borne r                
]

print("\nOptimizing root length...", flush=True)

length_result = differential_evolution(
    objective_length,
    bounds=length_bounds,
    seed=SEED,
    maxiter=100,
    callback=make_callback("Length"),
    polish=False,
)

print("Root length optimization finished.", flush=True)


# Save optimized length parameters.
length_parameters = {
    "lateral": {
        "ln": length_result.x[0],
        "r": length_result.x[1],                                                                                       
    },
    "basal": {
        "ln": length_result.x[2],
        "r": length_result.x[3],
    },
    "shootborne": {
        "r": length_result.x[4],
    },
}
modify_xml_parameters(
    PATH / f"{BASE_NAME}.xml",
    PATH / f"{TARGET_NAME}.xml",
    length_parameters,
)


# =============================================================================
# 2. Optimize root diameter
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

sys.exit()
# =============================================================================
# Final simulation
# =============================================================================

rs = pb.RootSystem()
rs.readParameters(str(xml_file))
rs.setSeed(0)
rs.setGeometry(pb.SDF_PlantContainer(10, 10, 60, False))
rs.initializeLB(5, 4)

plot_times = np.arange(times[-1] + 1)

sim_length = np.zeros(len(plot_times))
sim_diam = np.zeros(len(plot_times))

for day in range(1, plot_times[-1] + 1):

    rs.simulate(1, True)

    length = np.asarray(rs.getParameter("length"))
    radius = np.asarray(rs.getParameter("radius"))

    sim_length[day] = np.sum(length)
    sim_diam[day] = (
        np.sum(length * radius) / np.sum(length) * 20
    )


# =============================================================================
# Plot
# =============================================================================

fig, axes = plt.subplots(
    2, 1,
    figsize=(10, 10),
    sharex=True,
)

axes[0].plot(
    plot_times,
    sim_length,
    "--",
    label="Simulation",
)

axes[0].errorbar(
    times,
    real_length,
    yerr=real_length_SE,
    fmt="o",
    color="r",
    label="Measurement",
)

axes[0].set_ylabel("Root length")
axes[0].legend(loc="upper left")


axes[1].plot(
    plot_times,
    sim_diam,
    "--",
    label="Simulation",
)

axes[1].errorbar(
    times,
    real_diam,
    yerr=real_diam_SE,
    fmt="o",
    color="r",
    label="Measurement",
)

axes[1].set_xlabel("Days after sowing")
axes[1].set_ylabel("Root diameter [mm]")
axes[1].set_ylim(0, 0.5)
axes[1].legend(loc="upper left")

plt.tight_layout()
plt.show()
