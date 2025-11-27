import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. Load data
# ============================================================

excel_path = "Input_Data_VR_eff_calibration.xlsx"

kom = pd.read_excel(excel_path, sheet_name="Kompost",   engine="openpyxl")
kon = pd.read_excel(excel_path, sheet_name="Kontrolle", engine="openpyxl")


# ============================================================
# 1a. van Genuchten parameterization for depth-dependent θ_p(z)
# ============================================================

# Suction at plastic limit (cm)
H_P = -15000.0

# Layer boundaries (top to bottom, absolute depth in cm)
z_breaks = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]

# van Genuchten parameters per layer (top to bottom) for Kompost ("AF")
af_params_top_to_bot = [
    [0.056, 0.372, 0.0162, 4.48, 218.0],
    [0.056, 0.372, 0.0162, 4.48, 218.0],
    [0.047, 0.386, 0.0191, 3.82, 529.0],
    [0.001, 0.707, 0.0452, 1.41, 538.0],
    [0.033, 0.340, 0.0127, 2.046, 96.0],
    [0.019, 0.253, 0.0120, 1.922, 36.67],
]

# van Genuchten parameters per layer (top to bottom) for Kontrolle ("NF")
nf_params_top_to_bot = [
    [0.048, 0.351, 0.0190, 4.887, 421.67],
    [0.048, 0.351, 0.0190, 4.887, 421.67],
    [0.062, 0.371, 0.0160, 3.573, 373.3],
    [0.046, 0.368, 0.0152, 2.736, 105.67],
    [0.021, 0.304, 0.01590, 1.89, 127.67],
    [0.028, 0.243, 0.0124, 1.999, 7.89],
]


def vg_theta(h, theta_r, theta_s, alpha, n):
    """
    Compute volumetric water content θ(h) using van Genuchten equation.
    h must be negative (suction, cm).
    """
    h = np.abs(h)
    m = 1.0 - 1.0 / n
    return theta_r + (theta_s - theta_r) / np.power((1.0 + (alpha * h)**n), m)


def compute_theta_p_profile(depth_grid, z_breaks, vg_params):
    """
    Compute depth-dependent plastic-limit θ_p(z) using van Genuchten parameters.

    We map each depth to a layer defined by ``z_breaks`` (top to bottom) and
    assign a constant θ_p per layer:
        θ_p = θ(h = H_P) with H_P given in cm (negative suction).

    Parameters
    ----------
    depth_grid : array-like
        Depths where PR is defined (e.g. 0, 2, 4, ... or 0, -2, -4, ...).
    z_breaks : list of float
        Layer boundaries (top to bottom) in cm. Can be positive or negative;
        the mapping is done on absolute depth.
    vg_params : list of [θr, θs, α, n, Ks]
        van Genuchten parameters for each layer (top to bottom).
    """
    depth = np.array(depth_grid, dtype=float)
    depth_abs = np.abs(depth)
    theta_p = np.zeros_like(depth, dtype=float)

    for i in range(len(z_breaks) - 1):
        z_top = z_breaks[i]
        z_bot = z_breaks[i+1]

        # Work with absolute depths (cm), independent of sign convention
        z_top_abs = min(abs(z_top), abs(z_bot))
        z_bot_abs = max(abs(z_top), abs(z_bot))

        mask = (depth_abs >= z_top_abs) & (depth_abs < z_bot_abs)
        if not np.any(mask):
            continue

        theta_r, theta_s, alpha, n, Ks = vg_params[i]
        theta_p_layer = vg_theta(H_P, theta_r, theta_s, alpha, n)
        theta_p[mask] = theta_p_layer

    return theta_p


def build_depth_grid_profiles(df):
    """
    Build depth-resolved profiles directly from the Excel sheet,
    mapping SWC and BD as piecewise-constant layers and BD intervals.

    Returns depth, PR0, theta0 (SWC), rho0 (BD) on the PR depth grid.
    """
    # high-resolution PR vs depth arrays
    depth = df["Tiefe Penetration Resistance (cm)"].to_numpy()
    pr0   = df["Penetration Resistance (MPa)"].to_numpy()

    # extract SWC and BD dicts from first row
    swc_cols = [c for c in df.columns if c.startswith("SWC_")]
    bd_cols  = [c for c in df.columns if c.startswith("BD_")]

    swc_dict = {c: float(df[c].iloc[0]) for c in swc_cols}
    bd_dict  = {c: float(df[c].iloc[0]) for c in bd_cols}

    # ---------- SWC(z) as stepwise profile (Option A) ----------
    theta0 = np.zeros_like(depth, dtype=float)

    # keys like SWC_10cm, SWC_20cm, ...
    swc_depths = sorted(
        (int(k.split("_")[1].replace("cm", "")), v)
        for k, v in swc_dict.items()
    )
    dvals = [d for d, _ in swc_depths]    # [10, 20, 30, 40, 60, 100]
    sval  = {d: v for d, v in swc_depths} # mapping depth -> SWC value

    for i, z in enumerate(depth):
        # Assume z >= 0 => depth
        if z < dvals[0]:
            theta0[i] = sval[dvals[0]]       # above 10 cm -> SWC_10
        elif z < dvals[1]:
            theta0[i] = sval[dvals[0]]       # 10–20 cm -> SWC_10
        elif z < dvals[2]:
            theta0[i] = sval[dvals[1]]       # 20–30 -> SWC_20
        elif z < dvals[3]:
            theta0[i] = sval[dvals[2]]       # 30–40 -> SWC_30
        elif z < dvals[4]:
            theta0[i] = sval[dvals[3]]       # 40–60 -> SWC_40
        elif z < dvals[5]:
            theta0[i] = sval[dvals[4]]       # 60–100-> SWC_60
        else:
            theta0[i] = sval[dvals[5]]       # >100  -> SWC_100

    # ---------- BD(z) as stepwise profile from BD intervals ----------
    rho0 = np.zeros_like(depth, dtype=float)

    # keys like BD_0_30_gcm3
    intervals = []
    for k, v in bd_dict.items():
        parts = k.split("_")  # ["BD", "0", "30", "gcm3"]
        z1 = int(parts[1])
        z2 = int(parts[2])
        intervals.append((z1, z2, v))

    # Sort by top depth
    intervals = sorted(intervals, key=lambda x: x[0])

    for i, z in enumerate(depth):
        for (z1, z2, v) in intervals:
            if (z >= z1) and (z < z2):
                rho0[i] = v
                break
        else:
            # if no interval matches, use last or first as fallback
            if z < intervals[0][0]:
                rho0[i] = intervals[0][2]
            else:
                rho0[i] = intervals[-1][2]

    return depth, pr0, theta0, rho0


def build_depth_grid_profiles_interp(df):
    """
    Same as build_depth_grid_profiles, but:
      - SWC is linearly interpolated between SWC_xxcm depths
      - BD  is linearly interpolated between layer midpoints

    This yields smoother profiles while respecting the original depth structure:
      depth(z), PR0(z), theta0(z), rho0(z).

    Assumptions:
      - SWC_xxcm columns contain ONE real SWC value per depth,
        repeated in all rows. We use the first row.
      - BD_a_b_gcm3 columns contain ONE real BD per interval [a,b),
        repeated in all rows. We use the first row.
      - SWC is mapped to depth as piecewise-constant layers (Option A).
      - BD is mapped stepwise to the PR depth grid according to its intervals.
    """
    # high-resolution PR vs depth arrays
    depth = df["Tiefe Penetration Resistance (cm)"].to_numpy()
    pr0   = df["Penetration Resistance (MPa)"].to_numpy()

    # extract SWC and BD dicts from first row
    swc_cols = [c for c in df.columns if c.startswith("SWC_")]
    bd_cols  = [c for c in df.columns if c.startswith("BD_")]

    swc_dict = {c: float(df[c].iloc[0]) for c in swc_cols}
    bd_dict  = {c: float(df[c].iloc[0]) for c in bd_cols}

    # ---------- SWC: interpolate between SWC_xxcm depths ----------
    # keys like SWC_10cm, SWC_20cm, ...
    swc_depths = sorted(
        (int(c.split("_")[1].replace("cm", "")), float(df[c].iloc[0]))
        for c in swc_cols
    )
    swc_z = np.array([z for z, _ in swc_depths], dtype=float)
    swc_v = np.array([v for _, v in swc_depths], dtype=float)
    theta0 = np.interp(depth, swc_z, swc_v)   # extrapolation = end values

    # --- BD: interpolation over layer midpoints (0–30,30–45,...) ---
    bd_layers = []
    for c in bd_cols:
        parts = c.split("_")   # ["BD", "0", "30", "gcm3"]
        z1 = int(parts[1])
        z2 = int(parts[2])
        v  = float(df[c].iloc[0])
        zmid = 0.5 * (z1 + z2)
        bd_layers.append((zmid, v))

    bd_layers = sorted(bd_layers, key=lambda x: x[0])
    bd_z = np.array([z for z, _ in bd_layers], dtype=float)
    bd_v = np.array([v for _, v in bd_layers], dtype=float)
    rho0 = np.interp(depth, bd_z, bd_v)

    return depth, pr0, theta0, rho0


# --- Option A: stepwise layers (old method) ---
# depth_kom, pr0_kom, theta0_kom, rho0_kom = build_depth_grid_profiles(kom)
# depth_kon, pr0_kon, theta0_kon, rho0_kon = build_depth_grid_profiles(kon)

# --- Option B: interpolated SWC & BD (new method) ---
depth_kom, pr0_kom, theta0_kom, rho0_kom = build_depth_grid_profiles_interp(kom)
depth_kon, pr0_kon, theta0_kon, rho0_kon = build_depth_grid_profiles_interp(kon)

# ============================================================
# 2a. Depth-dependent plastic-limit θ_p(z) from VG parameters
# ============================================================

# Kompost (amended, "AF") and Kontrolle ("NF") θ_p(z) profiles
theta_p_kom = compute_theta_p_profile(depth_kom, z_breaks, af_params_top_to_bot)
theta_p_kon = compute_theta_p_profile(depth_kon, z_breaks, nf_params_top_to_bot)


# ============================================================
# 3. Moisture-only scaling (your original simple model)
# ============================================================

def scale_pr_profile(pr0, theta0, theta_new, c_theta=3.0):
    """
    Scale baseline PR profile pr0 from baseline theta0
    to a new water content theta_new using:

      PR(z,t) = PR0(z) * (theta0(z) / theta_new(z))^c_theta
    """
    pr0       = np.asarray(pr0, dtype=float)
    theta0    = np.asarray(theta0, dtype=float)
    theta_new = np.asarray(theta_new, dtype=float)

    theta_new_safe = np.maximum(theta_new, 1e-6)
    moisture_factor = (theta0 / theta_new_safe) ** c_theta
    pr_scaled = pr0 * moisture_factor
    return pr_scaled


# Define simple scenarios (still useful to see)
theta_dry_kom = theta0_kom * 0.7   # 30% drier
theta_wet_kom = theta0_kom * 1.3   # 30% wetter

theta_dry_kon = theta0_kon * 0.7
theta_wet_kon = theta0_kon * 1.3

c_theta = 3.0  # moisture sensitivity exponent

pr_dry_kom = scale_pr_profile(pr0_kom, theta0_kom, theta_dry_kom, c_theta=c_theta)
pr_wet_kom = scale_pr_profile(pr0_kom, theta0_kom, theta_wet_kom, c_theta=c_theta)
pr_dry_kon = scale_pr_profile(pr0_kon, theta0_kon, theta_dry_kon, c_theta=c_theta)
pr_wet_kon = scale_pr_profile(pr0_kon, theta0_kon, theta_wet_kon, c_theta=c_theta)


# ============================================================
# 4. Hybrid Vaz-based model (BD + moisture)
# ============================================================

def vaz_pr(theta_v, rho_b,
           theta_p=0.05,  # PWP (approx), adjust if you have better data
           rho_s=2.65,    # particle density [g/cm3]
           rho_bmin=None,
           rho_bmax=None):
    """
    Vaz et al. / Morandage et al. PR model:

      PR = exp(1.5 + 2.18 * rho_b_star - 4 * s_p)

      s_p  = (theta_v - theta_p) / (theta_s - theta_p)
      theta_s = 1 - rho_b / rho_s
      rho_b_star = (rho_b - rho_bmin) / (rho_bmax - rho_bmin)
    """
    theta_v = np.asarray(theta_v, dtype=float)
    rho_b   = np.asarray(rho_b,   dtype=float)

    # Saturated water content based on BD and particle density
    theta_s = 1.0 - rho_b / rho_s

    # Plasticity saturation
    theta_p_arr = np.asarray(theta_p, dtype=float)
    denom = np.maximum(theta_s - theta_p_arr, 1e-6)
    s_p = (theta_v - theta_p_arr) / denom
    s_p = np.clip(s_p, 0.0, 1.0)

    # BD normalization
    if rho_bmin is None:
        rho_bmin = float(np.nanmin(rho_b))
    if rho_bmax is None:
        rho_bmax = float(np.nanmax(rho_b))
    drho = max(rho_bmax - rho_bmin, 1e-6)
    rho_b_star = (rho_b - rho_bmin) / drho

    PR = np.exp(1.5 + 2.18 * rho_b_star - 4.0 * s_p)
    return PR


def smooth_moving_average(x, window=5):
    """
    Simple centered moving average for smoothing k(z).
    """
    x = np.asarray(x, dtype=float)
    if window <= 1:
        return x
    pad = window // 2
    x_pad = np.pad(x, pad, mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(x_pad, kernel, mode="valid")


def hybrid_pr(theta_v, rho_b, k_depth,
              theta_p=0.05, rho_s=2.65,
              rho_bmin=None, rho_bmax=None):
    """
    Hybrid model:
        PR_hybrid(z) = k(z) * PR_Vaz(theta_v(z), rho_b(z))
    """
    pr_theory = vaz_pr(theta_v, rho_b,
                       theta_p=theta_p, rho_s=rho_s,
                       rho_bmin=rho_bmin, rho_bmax=rho_bmax)
    return k_depth * pr_theory


# ----- 4.1 Compute Vaz-based PR and depth-correction k(z) -----

# Use global BD min/max for normalization across both profiles
all_bd_values = np.concatenate([rho0_kom, rho0_kon])
rho_bmin_global = float(np.nanmin(all_bd_values))
rho_bmax_global = float(np.nanmax(all_bd_values))

eps = 1e-6

# Kompost
pr_theory_kom = vaz_pr(theta0_kom, rho0_kom,
                       theta_p=theta_p_kom, rho_s=2.65,
                       rho_bmin=rho_bmin_global,
                       rho_bmax=rho_bmax_global)
k_kom = pr0_kom / (pr_theory_kom + eps)
k_kom_smooth = smooth_moving_average(k_kom, window=3)

# Kontrolle
pr_theory_kon = vaz_pr(theta0_kon, rho0_kon,
                       theta_p=theta_p_kon, rho_s=2.65,
                       rho_bmin=rho_bmin_global,
                       rho_bmax=rho_bmax_global)
k_kon = pr0_kon / (pr_theory_kon + eps)
k_kon_smooth = smooth_moving_average(k_kon, window=3)

# Hybrid baseline reconstruction (should be close to measured PR)
pr_hybrid_kom = hybrid_pr(theta0_kom, rho0_kom, k_kom_smooth,
                          theta_p=theta_p_kom, rho_s=2.65,
                          rho_bmin=rho_bmin_global,
                          rho_bmax=rho_bmax_global)

pr_hybrid_kon = hybrid_pr(theta0_kon, rho0_kon, k_kon_smooth,
                          theta_p=theta_p_kon, rho_s=2.65,
                          rho_bmin=rho_bmin_global,
                          rho_bmax=rho_bmax_global)


import json  # at top of file

# ============================================================
# 4.1c Export high-resolution k(z) profiles for coupling
#      → Kompost and Kontrolle separately
# ============================================================

max_sim_depth = 90.0  # cm

def extend_profile_to_90cm(depth, k_profile, max_depth=90.0):
    """
    Extend a depth–k(z) profile to max_depth (cm) by holding the
    last k-value constant below the deepest calibrated depth.
    """
    depth = np.asarray(depth, dtype=float)
    k_profile = np.asarray(k_profile, dtype=float)

    last_depth = float(depth[-1])
    last_k     = float(k_profile[-1])

    if last_depth >= max_depth:
        # already deep enough, just truncate if necessary
        mask = depth <= max_depth
        return depth[mask], k_profile[mask]

    extra_depths = np.arange(last_depth + 1.0, max_depth + 1.0, 1.0)
    extra_k      = np.full_like(extra_depths, last_k, dtype=float)

    depth_ext = np.concatenate([depth, extra_depths])
    k_ext     = np.concatenate([k_profile, extra_k])
    return depth_ext, k_ext


# --- Kompost profile ---
depth_kom_ext, k_kom_ext = extend_profile_to_90cm(depth_kom, k_kom_smooth, max_depth=max_sim_depth)

export_kom = {
    "depth_cm": depth_kom_ext.tolist(),
    "k_profile": k_kom_ext.tolist(),
    "note": "High-resolution k(z) derived from Kompost PR calibration and extended to 90 cm using last k-value"
}

with open("PR_calibration_kompost.json", "w") as f:
    json.dump(export_kom, f, indent=2)


# --- Kontrolle profile ---
depth_kon_ext, k_kon_ext = extend_profile_to_90cm(depth_kon, k_kon_smooth, max_depth=max_sim_depth)

export_kon = {
    "depth_cm": depth_kon_ext.tolist(),
    "k_profile": k_kon_ext.tolist(),
    "note": "High-resolution k(z) derived from Kontrolle PR calibration and extended to 90 cm using last k-value"
}

with open("PR_calibration_kontrolle.json", "w") as f:
    json.dump(export_kon, f, indent=2)






# Example BD modification: change Kompost BD in 45–60 cm to 1.4 g/cm3
# ============================================================
# 4.1b Generate BD scenarios correctly (modify BD layer, then
#      interpolate EXACTLY as for baseline BD)
# ============================================================

# First, reconstruct the BD layer structure from the Kompost sheet
bd_cols = [c for c in kom.columns if c.startswith("BD_")]

bd_layers = []
for c in bd_cols:
    parts = c.split("_")   # ["BD","0","30","gcm3"]
    z1 = int(parts[1])
    z2 = int(parts[2])
    v  = float(kom[c].iloc[0])
    bd_layers.append((z1, z2, v))

bd_layers = sorted(bd_layers, key=lambda x: x[0])

# Extract midpoints and layer values (this matches build_depth_grid_profiles_interp)
bd_z = np.array([0.5*(z1+z2) for (z1, z2, v) in bd_layers], dtype=float)
bd_v = np.array([v for (_, _, v) in bd_layers], dtype=float)

# ----------------------------------------------
# Generate 4 BD scenarios between 0.84–1.64 g/cm³
# modifying ONLY the 45–60 cm layer
# ----------------------------------------------
bd_start  = 0.84
bd_end    = 1.64
bd_values = np.linspace(bd_start, bd_end, 4)

pr_mod_scenarios = []  # store (bd_value, PR_profile)

for bd_val in bd_values:
    # copy original layer values
    bd_v_mod = bd_v.copy()

    # find the 45–60 layer index
    # (in your Kompost BD: layers are 0–30, 30–45, 45–60, 60–75, 75–90)
    idx_45_60 = None
    for i, (z1, z2, v) in enumerate(bd_layers):
        if (z1 == 45) and (z2 == 60):
            idx_45_60 = i
            break
    if idx_45_60 is None:
        raise ValueError("Could not find BD layer 45–60 cm in Kompost sheet.")

    # modify only this layer
    bd_v_mod[idx_45_60] = bd_val

    # interpolate modified BD to depth grid (same method as baseline)
    rho_mod_kom = np.interp(depth_kom, bd_z, bd_v_mod)

    # compute hybrid PR for this scenario
    pr_mod = hybrid_pr(theta0_kom, rho_mod_kom, k_kom_smooth,
                       theta_p=theta_p_kom, rho_s=2.65,
                       rho_bmin=rho_bmin_global,
                       rho_bmax=rho_bmax_global)

    pr_mod_scenarios.append((bd_val, pr_mod))


# ============================================================
# 4.2 Hybrid moisture scenarios (±30%) – NEW
# ============================================================

pr_dry_kom_hyb = hybrid_pr(theta_dry_kom, rho0_kom, k_kom_smooth,
                           theta_p=theta_p_kom, rho_s=2.65,
                           rho_bmin=rho_bmin_global,
                           rho_bmax=rho_bmax_global)

pr_wet_kom_hyb = hybrid_pr(theta_wet_kom, rho0_kom, k_kom_smooth,
                           theta_p=theta_p_kom, rho_s=2.65,
                           rho_bmin=rho_bmin_global,
                           rho_bmax=rho_bmax_global)

pr_dry_kon_hyb = hybrid_pr(theta_dry_kon, rho0_kon, k_kon_smooth,
                           theta_p=theta_p_kon, rho_s=2.65,
                           rho_bmin=rho_bmin_global,
                           rho_bmax=rho_bmax_global)

pr_wet_kon_hyb = hybrid_pr(theta_wet_kon, rho0_kon, k_kon_smooth,
                           theta_p=theta_p_kon, rho_s=2.65,
                           rho_bmin=rho_bmin_global,
                           rho_bmax=rho_bmax_global)


# ============================================================
# 5. Plot: Profiles & hybrid scaling (final version)
# ============================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# ---------------- Row 1 — Profiles (PR, SWC, BD) ----------------

# (1) PR profiles (measured + hybrid)
ax = axes[0, 0]
ax.plot(pr0_kom,        depth_kom, color="black", label="Kompost (measured)")
ax.plot(pr0_kon,        depth_kon, color="gray",  label="Kontrolle (measured)")
ax.plot(pr_hybrid_kom,  depth_kom, color="green", linestyle="--", label="Kompost hybrid")
ax.plot(pr_hybrid_kon,  depth_kon, color="lime",  linestyle="--", label="Kontrolle hybrid")
ax.axvline(2.0, linestyle="--", color="black", label="2 MPa")
ax.set_xlabel("PR (MPa)")
ax.set_ylabel("Depth (cm)")
ax.set_title("PR Profiles (measured vs hybrid)")
ax.set_ylim(0, max(depth_kom.max(), depth_kon.max()))
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.set_xlim(0, 20)
ax.legend()

# (2) SWC profiles
ax = axes[0, 1]
ax.plot(theta0_kom, depth_kom, label="Kompost SWC", color="blue")
ax.plot(theta0_kon, depth_kon, label="Kontrolle SWC", color="red")
ax.set_xlabel("SWC (cm³/cm³)")
ax.set_ylabel("Depth (cm)")
ax.set_title("Soil water content profiles")
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.legend()

# (3) BD profiles
ax = axes[0, 2]
ax.plot(rho0_kom, depth_kom, label="Kompost BD", color="blue")
ax.plot(rho0_kon, depth_kon, label="Kontrolle BD", color="red")
ax.set_xlabel("BD (g/cm³)")
ax.set_ylabel("Depth (cm)")
ax.set_title("Bulk density profiles")
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.legend()

# ---------------- Row 2 — Scenarios & k(z) ----------------

# (4) Kompost – hybrid moisture + hybrid BD scenarios
ax = axes[1, 0]
ax.plot(pr0_kom,        depth_kom, color="black", label="Kompost (measured)")
ax.plot(pr_hybrid_kom,  depth_kom, color="green", linestyle="--", label="Kompost hybrid")

# moisture scenarios (hybrid)
ax.plot(pr_dry_kom_hyb, depth_kom, color="blue",  linestyle=":", label="Kompost hybrid dry (-30%)")
ax.plot(pr_wet_kom_hyb, depth_kom, color="cyan",  linestyle=":", label="Kompost hybrid wet (+30%)")

# BD scenarios
colors_bd = ["purple", "magenta", "orange", "brown"]
for (bd_val, pr_mod), col in zip(pr_mod_scenarios, colors_bd):
    ax.plot(pr_mod, depth_kom, color=col, linestyle="-.",
            label=f"BD 45–60 = {bd_val:.2f} g/cm³")

ax.axvline(2.0, linestyle="--", color="black", label="2 MPa")
ax.set_xlabel("PR (MPa)")
ax.set_ylabel("Depth (cm)")
ax.set_title("Kompost: Hybrid PR scenarios (moisture & BD)")
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.set_xlim(0, 20)
ax.legend(fontsize=8)

# (5) Kontrolle – hybrid moisture scenarios
ax = axes[1, 1]
ax.plot(pr0_kon,        depth_kon, color="black", label="Kontrolle (measured)")
ax.plot(pr_hybrid_kon,  depth_kon, color="green", linestyle="--", label="Kontrolle hybrid")

ax.plot(pr_dry_kon_hyb, depth_kon, color="blue", linestyle=":", label="Kontrolle hybrid dry (-30%)")
ax.plot(pr_wet_kon_hyb, depth_kon, color="cyan", linestyle=":", label="Kontrolle hybrid wet (+30%)")

ax.axvline(2.0, linestyle="--", color="black", label="2 MPa")
ax.set_xlabel("PR (MPa)")
ax.set_ylabel("Depth (cm)")
ax.set_title("Kontrolle: Hybrid PR moisture scenarios")
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.set_xlim(0, 20)
ax.legend(fontsize=8)

# (6) Depth-dependent k(z)
ax = axes[1, 2]
ax.plot(k_kom,        depth_kom,  label="Kompost k (raw)", color="green", alpha=0.3)
ax.plot(k_kom_smooth, depth_kom,  label="Kompost k (smooth)",color="green")
ax.plot(k_kon,        depth_kon,  label="Kontrolle k (raw)", color="orange", alpha=0.3)
ax.plot(k_kon_smooth, depth_kon,  label="Kontrolle k (smooth)", color="orange")

ax.axvline(1.0, linestyle="--", color="black", label="k = 1")

ax.set_xlabel("Correction factor k(z)")
ax.set_ylabel("Depth (cm)")
ax.set_title("Depth-dependent correction factors k(z)")
ax.invert_yaxis()
ax.grid(True, linestyle=":", linewidth=0.5)
ax.legend()
plt.savefig("Penetration_resistance_adaption.png")
plt.tight_layout()
