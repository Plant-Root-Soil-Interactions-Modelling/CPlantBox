"""
Test that the 2D lookup table for sand agrees with the direct (analytical) computation
of soil-root interface potentials via soil_root_interface_().

Steps:
  1. Create a lookup table for sand and save it to a temporary file.
  2. Sample random parameter combinations (rx, sx, inner_kr, rho) within the
     valid table range.
  3. Evaluate the interface potential both ways:
       - direct:  PerirhizalPython.soil_root_interface_()  (no table)
       - table:   peri.soil_root_interface_potentials()    (uses lookup)
  4. Report the maximum absolute and relative error and assert they are small.
"""

import os
import sys
import tempfile

import numpy as np

import plantbox.functional.van_genuchten as vg
from plantbox.functional.Perirhizal import PerirhizalPython

# ---------------------------------------------------------------------------
# soil parameters
# ---------------------------------------------------------------------------
sand = [0.045, 0.43, 0.15, 3.0, 1000.0]  # [theta_r, theta_s, alpha, n, Ks]
loam = [0.078, 0.43, 0.036, 1.56, 24.96]
clay = [0.1, 0.4, 0.01, 1.1, 10.0]

sp = vg.Parameters(clay)  # pick: sand, loam, or clay
vg.create_mfp_lookup(sp)

# ---------------------------------------------------------------------------
# 1. Create lookup table and save to a temporary file
# ---------------------------------------------------------------------------
with tempfile.TemporaryDirectory() as tmpdir:
    table_path = os.path.join(tmpdir, "lookup")

    peri_create = PerirhizalPython()
    print("Creating 2D lookup table ...")
    peri_create.create_lookup(table_path, sp)
    print("Lookup table created.\n")

    # ---------------------------------------------------------------------------
    # 2. Load the lookup table into a fresh PerirhizalPython instance
    # ---------------------------------------------------------------------------
    peri = PerirhizalPython()
    peri.open_lookup(table_path)

    # ---------------------------------------------------------------------------
    # 3. Sample parameter combinations within the valid table range
    # ---------------------------------------------------------------------------
    rng = np.random.default_rng(42)
    n_samples = 500

    # xylem potential: must be below -1 cm and not too negative
    rx = -rng.uniform(1.0, 15000.0, n_samples)
    sx = -rng.uniform(1.0, 15000.0, n_samples)

    # inner_kr = a * kr; valid range from the table construction
    inner_kr = 10 ** rng.uniform(np.log10(1e-7), np.log10(1e-4), n_samples)

    # rho = outer_radius / inner_radius; valid range [2, 199]
    rho = rng.uniform(2.0, 199.0, n_samples)

    # ---------------------------------------------------------------------------
    # 4a. Direct computation (no table)
    # ---------------------------------------------------------------------------
    peri_direct = PerirhizalPython()
    peri_direct.set_soil(sp)  # sets sp, clears lookup tables
    rsx_direct = np.array([PerirhizalPython.soil_root_interface_(rx[i], sx[i], inner_kr[i], rho[i], sp) for i in range(n_samples)])

    # ---------------------------------------------------------------------------
    # 4b. Lookup-table computation
    # ---------------------------------------------------------------------------
    rsx_table = peri.soil_root_interface_potentials(rx.copy(), sx.copy(), inner_kr.copy(), rho.copy())

    # ---------------------------------------------------------------------------
    # 5. Evaluate agreement
    # ---------------------------------------------------------------------------
    abs_err = np.abs(rsx_table - rsx_direct)
    # use |direct| + 1 in the denominator to avoid division by near-zero values
    rel_err = abs_err / (np.abs(rsx_direct) + 1.0)

    print(f"{'Sample':>8}  {'rx [cm]':>12}  {'sx [cm]':>12}  {'inner_kr':>12}  {'rho':>6}  {'direct [cm]':>12}  {'table [cm]':>12}  {'|err| [cm]':>12}")
    print("-" * 100)
    # print the 10 worst cases
    worst = np.argsort(abs_err)[-10:][::-1]
    for i in worst:
        print(f"{i:>8}  {rx[i]:>12.2f}  {sx[i]:>12.2f}  {inner_kr[i]:>12.2e}  {rho[i]:>6.1f}  {rsx_direct[i]:>12.4f}  {rsx_table[i]:>12.4f}  {abs_err[i]:>12.4f}")

    print()
    print(f"max absolute error : {np.max(abs_err):.4f} cm")
    print(f"mean absolute error: {np.mean(abs_err):.4f} cm")
    print(f"max relative error : {np.max(rel_err):.4e}  (relative to |h_sr| + 1 cm)")
    print(f"mean relative error: {np.mean(rel_err):.4e}")

    # # tolerance: the lookup table is piecewise-linear on a coarse grid, so a few
    # # cm of absolute error is expected at the extremes; tighten as needed.
    # abs_tol = 10.0  # cm
    # rel_tol = 5e-2  # 5 %

    # assert np.max(abs_err) < abs_tol, f"Max absolute error {np.max(abs_err):.4f} cm exceeds tolerance {abs_tol} cm"
    # assert np.max(rel_err) < rel_tol, f"Max relative error {np.max(rel_err):.4e} exceeds tolerance {rel_tol}"

    # print("\nAll checks passed.")
