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
import csv
import os
from scipy import interpolate

# =============================================================================
# Settings
# =============================================================================

PATH = Path("rootsystem")
PATH_EXU = Path("data")
NAME = "RS_optimized_S_WT"
EXU_NAME = "exudation_rates"

SEED = 42
TARGET_ERROR = 0.001
tip = 3.5 #cm 

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


def write_csv(filename, headers, first_column, second_column=None, third_column=None):
    """
    Create or update a CSV file with 3 columns.

    Parameters
    ----------
    filename : str
        Path to the CSV file.

    headers : list[str]
        Three column headers, e.g. ["Name", "Value 1", "Value 2"].

    first_column : list
        Values for the first column when creating the file.

    second_column : list, optional
        Values to write/update in the second column.

    third_column : list, optional
        Values to write/update in the third column.
    """

    if len(headers) != 3:
        raise ValueError("headers must contain exactly 3 strings.")

    # ---------------------------------------------------------
    # CASE 1: CSV does not exist -> create it
    # ---------------------------------------------------------
    if not os.path.exists(filename):

        if second_column is None:
            second_column = [""] * len(first_column)

        if third_column is None:
            third_column = [""] * len(first_column)

        if len(second_column) != len(first_column):
            raise ValueError("second_column must have the same length as first_column.")

        if len(third_column) != len(first_column):
            raise ValueError("third_column must have the same length as first_column.")

        with open(filename, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)

            # Write headers
            writer.writerow(headers)

            # Write data
            for a, b, c in zip(first_column, second_column, third_column):
                writer.writerow([a, b, c])

        print(f"Created: {filename}")

    # ---------------------------------------------------------
    # CASE 2: CSV already exists -> update it
    # ---------------------------------------------------------
    else:

        with open(filename, "r", newline="", encoding="utf-8") as file:
            reader = csv.reader(file)
            rows = list(reader)

        if not rows:
            raise ValueError("CSV file is empty.")

        # Check that the CSV has 3 columns
        if len(rows[0]) != 3:
            raise ValueError("Existing CSV does not have exactly 3 columns.")

        # Number of data rows
        number_of_rows = len(rows) - 1

        # Update second column
        if second_column is not None:
            if len(second_column) != number_of_rows:
                raise ValueError(
                    f"second_column must contain {number_of_rows} values."
                )

            for i, value in enumerate(second_column, start=1):
                rows[i][1] = value

        # Update third column
        if third_column is not None:
            if len(third_column) != number_of_rows:
                raise ValueError(
                    f"third_column must contain {number_of_rows} values."
                )

            for i, value in enumerate(third_column, start=1):
                rows[i][2] = value

        # Write the updated CSV
        with open(filename, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerows(rows)

        print(f"Updated: {filename}")

# =============================================================================
# Simulation
# =============================================================================

def simulate(xml_file, params):

    rs = pb.RootSystem()
    rs.readParameters(str(xml_file))
    
    #get lmax of the different root types
    lmax = []
    for pp in rs.getRootRandomParameter():
        lmax.append(pp.lmax)


    rs.setSeed(0)
    tube_ = pb.SDF_PlantContainer(10, 10, 70, False) #tube with diameter 20, length 60 cm 
    tube = pb.SDF_RotateTranslate(tube_, 0, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 10))
    rs.setGeometry(tube)
    rs.setGeometry(pb.SDF_PlantContainer(10, 10, 60, False))
    rs.initializeLB(5, 4)

    simulated_exu = np.zeros(len(times))

    for day in range(1, times[-1] + 1):

        rs.simulate(1, True)

        if day not in times:
            continue

        i = np.where(times == day)[0][0]
        
        polylengths = np.asarray(rs.getParameter("length"))
        radii = np.asarray(rs.getParameter("radius"))
        types = np.asarray(rs.getParameter("type"))
        polylines = rs.getPolylines()
        
        kex_tip = params[i]
        kex_base = kex_tip/2

        sf = []
        for j in range(0, len(polylines)):
            a = radii[j]
            roottype = int(types[j])
            l_ = 0
            for k in range(0, len(polylines[j])-1):

                m = polylines[j][-1-k]
                n = polylines[j][-2-k]
                p0 = np.array([m.x, m.y, m.z])
                p1 = np.array([n.x, n.y, n.z])
                l  = np.linalg.norm(p0 - p1)
                l_ = l_+l
                #tip exudation rate 
                if (polylengths[j]-l_)<tip: #3.5cm - tip exudation 
                    kexu = kex_tip
                    c = 2
                #base exudation rate 
                else:
                    kexu = kex_base
                    c = 1
                #if growth has already stopped (95% of total length reached) 
                if polylengths[j]>=lmax[roottype]*0.99:
                    #print('REACHED')
                    kexu = kex_base
                    c = 1
                #if artificial shoot 
                if roottype == 0:
                    kexu = 0
                    c = 0

                sf.append(2 * np.pi * a * l * kexu) # mol/root segment / day


        simulated_exu[i] = np.sum(sf) # mol/plant/day

    return simulated_exu


def rmse(simulated, measured):
    return np.sqrt(np.mean(((simulated - measured) / measured) ** 2))


# =============================================================================
# Objective functions
# =============================================================================

def objective_exurates(params):

    simulated = simulate(
        PATH / f"{NAME}.xml", 
        params
    )

    return rmse(simulated, real_exu)


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

real_exu = data["Total exudation"].values[:len(times)]*24/10**3 #mmol/plant/h --> mol/plant/day
real_exu_SE = data["Total exudation SE"].values[:len(times)]*24/10**3 #mmol/plant/h --> mol/plant/day


# =============================================================================
# 1. Optimize root exudation rate
# =============================================================================

exu_bounds = [
    (1E-8, 1E-4),   # mol/cm^2/day
    (1E-8, 1E-4),   # mol/cm^2/day
]

print("\nOptimizing exudation rate...", flush=True)

exu_result = differential_evolution(
    objective_exurates,
    bounds=exu_bounds,
    seed=SEED,
    maxiter=500,
    callback=make_callback("Exudation"),
    polish=False,
)

print("Root exudation optimization finished.", flush=True)

# Save optimized exudation parameters.
exu_rates = np.asarray([exu_result.x[0],exu_result.x[1]])
write_csv(str(PATH_EXU / f"{EXU_NAME}.csv"), ['DAS','Loam', 'Sand'], times, None, exu_rates)


# =============================================================================
# Results
# =============================================================================

print("\n" + "=" * 50)
print("OPTIMIZATION RESULTS")
print("=" * 50)

for name, result in [
    ("Exudation", exu_result),
]:
    print(f"\n{name}:")
    print(f"  Success:        {result.success}")
    print(f"  Iterations:     {result.nit}")
    print(f"  Evaluations:    {result.nfev}")
    print(f"  Best error:     {result.fun:.6f}")
    print(f"  Message:        {result.message}")


# =============================================================================
# Final simulation
# =============================================================================
times = np.insert(times, 0, 0)
exu_rates= np.insert(exu_rates, 0,exu_rates[0])
f = interpolate.interp1d(times, exu_rates)

rs = pb.RootSystem()
rs.readParameters(str(PATH / f"{NAME}.xml"))
rs.setSeed(0)
tube_ = pb.SDF_PlantContainer(10, 10, 70, False) #tube with diameter 20, length 60 cm 
tube = pb.SDF_RotateTranslate(tube_, 0, pb.SDF_Axis.zaxis, pb.Vector3d(0, 0, 10))
rs.setGeometry(tube)
rs.initializeLB(5, 4)

#get lmax of the different root types
lmax = []
for pp in rs.getRootRandomParameter():
    lmax.append(pp.lmax)

plot_times = np.arange(times[-1] + 1)
sim_exu = np.zeros(len(plot_times))

for day in range(1, plot_times[-1] + 1):

    kex_tip = f(day)
    kex_base = kex_tip/2
    rs.simulate(1, True)
    
    polylengths = np.asarray(rs.getParameter("length"))
    radii = np.asarray(rs.getParameter("radius"))
    types = np.asarray(rs.getParameter("type"))
    polylines = rs.getPolylines()
    

    sf = []
    for j in range(0, len(polylines)):
        a = radii[j]
        roottype = int(types[j])
        l_ = 0
        for k in range(0, len(polylines[j])-1):

            m = polylines[j][-1-k]
            n = polylines[j][-2-k]
            p0 = np.array([m.x, m.y, m.z])
            p1 = np.array([n.x, n.y, n.z])
            l  = np.linalg.norm(p0 - p1)
            l_ = l_+l
            #tip exudation rate 
            if (polylengths[j]-l_)<tip:
                kexu = kex_tip
                c = 2
            #base exudation rate 
            else:
                kexu = kex_base
                c = 1
            #if growth has already stopped (95% of total length reached) 
            if polylengths[j]>=lmax[roottype]*0.99:
                #print('REACHED')
                kexu = kex_base
                c = 1
            #if artificial shoot 
            if roottype == 0:
                kexu = 0
                c = 0

            sf.append(2 * np.pi * a * l * kexu) # mol/root segment / day


    sim_exu[day] = np.sum(sf) # mol/plant/day


# =============================================================================
# Plot
# =============================================================================

fig, axes = plt.subplots(
    1, 1,
    figsize=(10, 10),
    sharex=True,
)

axes.plot(
    plot_times,
    sim_exu,
    "--",
    label="Simulation",
)

axes.errorbar(
    times[1:],
    real_exu,
    yerr=real_exu_SE,
    fmt="o",
    color="r",
    label="Measurement",
)

axes.set_ylabel("Root exudation")
axes.legend(loc="upper left")


plt.tight_layout()
plt.show()
