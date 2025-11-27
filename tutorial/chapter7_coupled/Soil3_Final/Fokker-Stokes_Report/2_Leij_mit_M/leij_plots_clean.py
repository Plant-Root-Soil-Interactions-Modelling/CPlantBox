
"""
Lean plotting/scenario harness (analytical only).
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import leij_core_clean_drift_left as lm2
import json

def summarize_psd(r, f):
    m0 = float(np.trapz(f, r))
    if m0 <= 0:
        return dict(phi0=0.0, r_m=float("nan"), sigma=float("nan"))
    p = f / m0
    ln_r = np.log(r)
    mu = float(np.trapz(ln_r * p, r))
    var = float(np.trapz((ln_r - mu) ** 2 * p, r))
    return dict(phi0=m0, r_m=float(math.exp(mu)), sigma=float(math.sqrt(max(var, 0.0))))

def run_settling_scenarios(
    p0: lm2.SoilParams,
    *,
    r_min_um: float = 0.02,
    r_max_um: float = 300.0,
    N_plot: int = 900,
    T_list = (0.5, 1.0, 1.5),
    eps_list = (0.08, 0.15, 0.25),
    phi_drop: float = 1.00,                 # target φ at T (per_T) or at Tmax (end_anchored)
    title: str = "Forward scenarios",
    save_png: str | None = "scenarios_psd.png",
    phi_mode: str = "per_T",                # NEW: "per_T" (old behavior) or "end_anchored"
    xi_power_on_term2: bool = True,        # keep your Eq.24 variant toggle
):
    import numpy as np, math, matplotlib.pyplot as plt

    r = np.geomspace(r_min_um, r_max_um, int(N_plot))
    f0 = lm2.kosugi_psd(r, p0)
    kappa_0 = calc_permeability_proxy(r, f0) # The reference permeability
    # NEW: one global decay rate if end-anchored
    if phi_mode == "end_anchored":
        T_end_max = max(T_list)
        M_over_u0w_global = 0.0 if abs(phi_drop-1.0) < 1e-12 else math.log(phi_drop)/T_end_max
    else:
        M_over_u0w_global = None

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.semilogx(r, f0, lw=2.8, color="k", label="initial PSD")

    



    for T in T_list:
        # OLD per-T behavior vs NEW end-anchored behavior
        if phi_mode == "per_T":
            M_over_u0w = 0.0 if abs(phi_drop-1.0) < 1e-12 else math.log(phi_drop)/T
        else:
            M_over_u0w = M_over_u0w_global

        for eps0 in eps_list:
            f = lm2.analytical_solution_rdep(
                r, p0, float(T),
                eps0=float(eps0),
                M_over_u0w=float(M_over_u0w),
                r_min_um=float(r_min_um),
                xi_power_on_term2=bool(xi_power_on_term2),
            )
            phi_frac = float(np.trapezoid(f, r)) / p0.phi_0
            # build speaking name from T, eps0
            if T == T_min and eps0 == eps_min:
                name = "geringe Setzung (eng)"
            elif T == T_min and eps0 == eps_max:
                name = "geringe Setzung (breit)"
            elif T == T_max and eps0 == eps_min:
                name = "starke Setzung (eng)"
            elif T == T_max and eps0 == eps_max:
                name = "starke Setzung (breit)"
            else:
                name = "mittlere Setzung"
            ax.semilogx(r, f, lw=2.0, ls="--", label=f"{name}, φ/φ₀={phi_frac:.2f}")


    ax.set_xlabel("r (µm)")
    ax.set_ylabel("f (µm⁻¹)")
    ax.set_xscale('log')
    ax.set_xlim(10**-3, 10**2)
    #ax.set_title(title + ("" if phi_mode=="per_T" else "  (end-anchored φ path)"))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=160)


def run_settling_scenarios_end_anchored(
    p0: lm2.SoilParams,
    *,
    r_min_um: float = 0.02,
    r_max_um: float = 300.0,
    N_plot: int = 900,
    T_list = (0.5, 1.0, 1.5),
    eps_list = (0.08, 0.15, 0.25),
    phi_drop_anchor: float = 0.64,   # φ(T_max)/φ0  (≤1)
    title: str = " ",
    save_png: str | None = "scenarios_psd_end_anchored.png",
    xi_power_on_term2: bool = False,
):
    import numpy as np, math, matplotlib.pyplot as plt
    r = np.geomspace(r_min_um, r_max_um, int(N_plot))
    f0 = lm2.kosugi_psd(r, p0)

    T_end_max = max(T_list)
    M_over_u0w_anchor = math.log(phi_drop_anchor) / T_end_max  # negative if drop

    fig, ax = plt.subplots(figsize=(9.5, 5.2))
    ax.semilogx(r, f0, lw=2.8, color="k", label="initial PSD")

    for T in T_list:
        for eps0 in eps_list:
            f = lm2.analytical_solution_rdep(
                r, p0, float(T),
                eps0=float(eps0),
                M_over_u0w=float(M_over_u0w_anchor),
                r_min_um=float(r_min_um),
                xi_power_on_term2=bool(xi_power_on_term2),
            )
            ax.semilogx(r, f, lw=2.0, ls="--",
                        label=f"T={T}, eps0={eps0}")

    ax.set_xlabel("r (µm)")
    ax.set_ylabel("f (µm⁻¹)")
    ax.set_xscale('log')
    ax.set_xlim(10**-3, 10**2)
    ax.set_title(title + f"  (φ(T) = φ0·{phi_drop_anchor}^(T/Tmax))")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=160)
        
def run_settling_scenarios_pairs(
    p0,
    pairs,
    *,
    r_min_um: float,
    r_max_um: float,
    N_plot: int,
    phi_mode: str,
    phi_value: float,
    title: str,
    save_png: str,
    xi_power_on_term2: bool = False,
    Ks_measured: float,
):
    import numpy as np
    import matplotlib.pyplot as plt

    r = np.geomspace(r_min_um, r_max_um, int(N_plot))

    # initial PSD
    f0 = lm2.kosugi_psd(r, p0)
    phi0 = float(np.trapezoid(f0, r))
    
    # reference permeability proxy from initial PSD
    kappa_0 = calc_permeability_proxy(r, f0)

    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.semilogx(r, f0, "k", lw=2.5, label="initiale PSD (gemessen)")

    # mins/max for legend classification
    T_min = min(t for t, _ in pairs)
    T_max = max(t for t, _ in pairs)
    eps_min = min(e for _, e in pairs)
    eps_max = max(e for _, e in pairs)

    # end-anchored mass factor
    if phi_mode == "end_anchored":
        if abs(phi_value - 1.0) < 1e-12:
            M_over_u0w = 0.0
        else:
            M_over_u0w = np.log(phi_value) / T_max
    else:
        M_over_u0w = 0.0

    for (T, eps0) in pairs:
        f = lm2.analytical_solution_rdep(
            r, p0, float(T),
            eps0=float(eps0),
            M_over_u0w=float(M_over_u0w),
            r_min_um=float(r_min_um),
            xi_power_on_term2=bool(xi_power_on_term2),
        )
        phiT = float(np.trapezoid(f, r))
        phi_frac = phiT / phi0

        # Estimate Ks scaling using Capillary Bundle Theory (Hagen-Poiseuille): Ks ~ ∫ r² f(r) dr
        kappa_new = calc_permeability_proxy(r, f)
        # Ratio of new permeability to old
        ks_ratio = kappa_new / max(kappa_0, 1e-12) 

        #  Adapted Ks value (e.g., 712 cm/day for sand, or your measured value)
        Ks_new = Ks_measured * ks_ratio


        # --- Legende (gleiche Form wie VG) ---
        # Setzung
        if abs(T - T_min) < 1e-9:
            s_txt = "geringe Setzung"
        elif abs(T - T_max) < 1e-9:
            s_txt = "starke Setzung"
        else:
            s_txt = "mittlere Setzung"

        # Heterogenisierung
        if abs(eps0 - eps_min) < 1e-9:
            h_txt = "geringe Heterogenisierung"
        elif abs(eps0 - eps_max) < 1e-9:
            h_txt = "starke Heterogenisierung"
        else:
            h_txt = "mittlere Heterogenisierung"

        label = f"{s_txt}, {h_txt}, φ_end/φ₀={phi_frac:.2f}"
        ax.semilogx(r, f, lw=2.0, ls="--", label=label)

    ax.set_xlabel("Porenradius r (µm)")
    ax.set_ylabel("PSD f(r)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=160)

import matplotlib.pyplot as plt
import numpy as np

def plot_psd_and_vg_side_by_side(
    p0,
    pairs,
    *,
    rho_s: float,
    phi_end_frac: float,
    theta_r: float,
    theta_s_initial: float,
    alpha0: float,
    n0: float,
    r_min_um: float = 0.001,
    r_max_um: float = 3000.0,
    N_plot_psd: int = 900,
    N_plot_vg: int = 1000,
    h_fit_min: float = 1.0,
    h_fit_max: float = 1500.0,
    h_plot_max: float = 1e5,
    xi_power_on_term2: bool = False,
    save_png: str = None,
):
    import numpy as np
    import matplotlib.pyplot as plt

    # --- Grids ---
    r = np.geomspace(r_min_um, r_max_um, int(N_plot_psd))
    h_fit  = np.geomspace(h_fit_min, h_fit_max, 200)
    h_plot = np.geomspace(1.0, h_plot_max, 400)

    # --- mins/max for naming ---
    T_min = min(t for t, _ in pairs)
    T_max = max(t for t, _ in pairs)
    eps_min = min(e for _, e in pairs)
    eps_max = max(e for _, e in pairs)

    # --- end-anchored mass rate (nur max T erreicht phi_end_frac) ---
    M_over_u0w = 0.0 if abs(phi_end_frac - 1.0) < 1e-12 else np.log(phi_end_frac) / T_max

    # --- figure ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.8))

    # -------------------- LEFT: PSD --------------------
    f0 = lm2.kosugi_psd(r, p0)
    phi0 = float(np.trapezoid(f0, r))
    
    theta_r_sand = 0.046
    theta_s_sand = 0.368
    alpha_sand   = 0.0152  # 1/cm
    n_sand       = 2.736
    p_sand = lm2.SoilParams.from_van_genuchten(
        theta_s_sand, theta_r_sand, alpha_sand, n_sand
    )
    f_sand = lm2.kosugi_psd(r, p_sand)
    
    
    ax1.axvspan(0.2, 50.0, color="0.7", alpha=0.5, zorder=0)
    ax1.semilogx(r, f0, "k", lw=2.5, label="initiale PSD (gemessen)")
    # dann Sand

    # legend entries für rechts sammeln
    legend_entries = [("initiale PSD (gemessen)", None, None, 1.0)]

    for (T, eps0) in pairs:
        f = lm2.analytical_solution_rdep(
            r, p0, float(T),
            eps0=float(eps0),
            M_over_u0w=float(M_over_u0w),
            r_min_um=float(r_min_um),
            xi_power_on_term2=bool(xi_power_on_term2),
        )
        phiT = float(np.trapezoid(f, r))
        phi_frac = phiT / phi0
        ax1.semilogx(r, f, lw=2.0, ls="--")

        # deutsche Szenarien
        if abs(T - T_min) < 1e-9:
            s_txt = "geringe Setzung"
        elif abs(T - T_max) < 1e-9:
            s_txt = "starke Setzung"
        else:
            s_txt = "mittlere Setzung"

        if abs(eps0 - eps_min) < 1e-9:
            h_txt = "geringe Heterogenisierung"
        elif abs(eps0 - eps_max) < 1e-9:
            h_txt = "starke Heterogenisierung"
        else:
            h_txt = "mittlere Heterogenisierung"

        lab = fr"{s_txt}, {h_txt}, $\phi_{{\mathrm{{end}}}}/\phi_0={phi_frac:.2f}$"

        legend_entries.append((lab, T, eps0, phi_frac))

    ax1.set_xlabel("Porenradius r (µm)")
    ax1.set_ylabel("Porengrößenverteilung f(r) (µm⁻¹)")
    #ax1.set_title("Porengrößenverteilungen (BD-ankert)")
    ax1.grid(True, which="both", alpha=0.25)
    # keine Legende links

    # -------------------- RIGHT: VG --------------------
    # initiale gemessene VG
    theta_ref = vg_theta(h_plot, theta_r, theta_s_initial, alpha0, n0)
    #ax2.plot(theta_ref, h_plot, "k", lw=2.5, label="initiale VG (Messung)")

        # 2) zusätzliche Referenz: unbehandelter Sand
    theta_r_sand = 0.046
    theta_s_sand = 0.368
    alpha_sand   = 0.0152   # 1/cm
    n_sand       = 2.736
    theta_sand_curve = vg_theta(h_plot, theta_r_sand, theta_s_sand, alpha_sand, n_sand)
    ax2.plot(theta_sand_curve, h_plot, color="0.4", lw=2.0, ls="-",
             label="Sandboden neben Meliorationsstreifen (Messung)")

    # T=0 PSD → VG (optional, gleiche Logik wie Szenarien)
    # (wir nehmen hier die p0-PSD, nicht f aus Szenarien)
    # --- T=0 PSD → VG (referenz) ---
    phi0_psd = float(np.trapezoid(f0, r))
    # residuales Wasser rein aus PSD @ 15000 cm ableiten
    theta_pwp_cap0 = float(theta_from_psd(r, f0, np.array([15000.0]), theta_r=0.0)[0])
    # HIER: direkt übernehmen, nicht max(...)
    theta_r0_used = theta_pwp_cap0
    theta_s0_used = theta_r0_used + phi0_psd
    theta_target0 = theta_from_psd(r, f0, h_fit, theta_r=0.0) + theta_r0_used
    w0 = np.ones_like(h_fit); w0[(h_fit >= 30) & (h_fit <= 500)] = 1.5
    alpha0_fit, n0_fit = fit_vg_to_theta(h_fit, theta_target0, theta_r0_used, theta_s0_used, weights=w0)
    theta_ref_psd = vg_theta(h_plot, theta_r0_used, theta_s0_used, alpha0_fit, n0_fit)
    ax2.plot(theta_ref_psd, h_plot, lw=2.5, color="k", label="Meliorationsstreifen mit Kompost (Messung)")


    for (label_text, T, eps0, _) in legend_entries[1:]:
        # Szenario-PSD erneut berechnen
        f = lm2.analytical_solution_rdep(
            r, p0, float(T),
            eps0=float(eps0),
            M_over_u0w=float(M_over_u0w),
            r_min_um=float(r_min_um),
            xi_power_on_term2=bool(xi_power_on_term2),
        )
        phiT = float(np.trapezoid(f, r))

        # --- HIER Variante 2: θr nie unter gemessenes θr ---
        theta_pwp_cap = float(theta_from_psd(r, f, np.array([15000.0]), theta_r=0.0)[0])
        theta_r_T = max(theta_r, theta_pwp_cap)   # <- das war vorher nicht überall so
        theta_s_T = theta_r_T + phiT
        
        # --- NEW: scenario bulk density ---
        theta_s_initial_leij = theta_s_initial  # = 0.707
        BD_new = scale_bulk_density_physically_corrected(theta_s_initial_leij, theta_s_T, BD_initial)


        theta_target = theta_from_psd(r, f, h_fit, theta_r=0.0) + theta_r_T
        w = np.ones_like(h_fit); w[(h_fit >= 30) & (h_fit <= 500)] = 1.5
        alpha_fit, n_fit = fit_vg_to_theta(h_fit, theta_target, theta_r_T, theta_s_T, weights=w)

        theta_curve = vg_theta(h_plot, theta_r_T, theta_s_T, alpha_fit, n_fit)
        ax2.plot(theta_curve, h_plot, lw=2.0, ls="--", label=label_text)

    ax2.set_xlabel("θ (cm³ cm⁻³)")
    ax2.set_ylabel("-Ψ (cm)")
    ax2.set_yscale("log")
    ax2.set_ylim(1.0, h_plot_max)
    ax2.set_xlim(0.0, theta_s_initial * 1.02)
    #ax2.set_title("Wasserhaltekurven (abgeleitet)")
    ax2.grid(True, which="both", alpha=0.25)
    ax2.legend(frameon=False, fontsize=8, loc="upper right")
    ax1.text(-0.1, 1.05, "a)", transform=ax1.transAxes,
             fontsize=12, fontweight='bold', va='top', ha='right')
    ax2.text(-0.1, 1.05, "b)", transform=ax2.transAxes,
             fontsize=12, fontweight='bold', va='top', ha='right')

    fig.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=300)


# --- Helpers for VG conversion/fitting ---------------------------------------

def _cumtrapz_np(y, x):
    import numpy as np
    y = np.asarray(y, dtype=float); x = np.asarray(x, dtype=float)
    dx = np.diff(x); mid = 0.5*(y[:-1] + y[1:])
    return np.concatenate([[0.0], np.cumsum(dx*mid)])

def theta_from_psd(r_um, f, h_cm, theta_r, A_cm2: float = 0.149):
    """
    θ(h) = θ_r + ∫_0^{A/h} f(r) dr, with r in µm, f in µm⁻¹, h in cm.
    A is Jurin constant for water in cm²; A≈0.149.
    """
    import numpy as np
    r = np.asarray(r_um, dtype=float)
    f = np.asarray(f, dtype=float)
    order = np.argsort(r)
    r = r[order]; f = f[order]
    cum = _cumtrapz_np(f, r)
    r_lim_um = (A_cm2 / np.maximum(h_cm, 1e-12)) * 1e4  # cm²/cm = cm → ×1e4 to µm
    return theta_r + np.interp(r_lim_um, r, cum, left=0.0, right=cum[-1])

def vg_theta(h, theta_r, theta_s, alpha_cm_inv, n):
    """
    van Genuchten water retention: θ(h) = θr + (θs-θr) / [1 + (α h)^n]^m, m=1-1/n
    """
    import numpy as np, math
    m = 1.0 - 1.0/max(n, 1.0000001)
    return theta_r + (theta_s - theta_r) / np.power(1.0 + np.power(alpha_cm_inv * np.maximum(h, 0.0), n), m)

def fit_vg_to_theta(h, theta_target, theta_r, theta_s,
                    alpha_grid=None, n_grid=None, weights=None):
    """
    Simple grid-search fit for (α, n) with θr, θs fixed.
    """
    import numpy as np
    if alpha_grid is None:
        alpha_grid = np.logspace(-4.5, -0.8, 80)  # 3.16e-5 .. 0.158 cm⁻¹
    if n_grid is None:
        n_grid = np.linspace(1.1, 3.0, 80)
    theta_target = np.asarray(theta_target, dtype=float)
    if weights is None:
        weights = np.ones_like(theta_target)
    best = (None, None, float("inf"))
    for a in alpha_grid:
        preds = np.stack([vg_theta(h, theta_r, theta_s, a, nn) for nn in n_grid], axis=0)
        err = np.mean(((preds - theta_target) * weights)**2, axis=1)
        j = int(np.argmin(err)); e = float(err[j])
        if e < best[2]:
            best = (float(a), float(n_grid[j]), e)
    return best[0], best[1]

def kosugi_equiv_from_psd(r_um, f):
    """Return (phi, sigma, r_m_um) by matching log-moments of the PSD."""
    import numpy as np
    r = np.asarray(r_um, dtype=float)
    f = np.asarray(f, dtype=float)
    phi = float(np.trapezoid(f, r))
    if phi <= 0:
        return 0.0, 0.0, float(np.nan)
    p = f / phi
    ln_r = np.log(r)
    mu  = float(np.trapezoid(ln_r * p, r))
    var = float(np.trapezoid((ln_r - mu)**2 * p, r))
    sigma = float(max(0.0, var)**0.5)
    r_m   = float(np.exp(mu))
    return phi, sigma, r_m

def calc_permeability_proxy(r, f):
    """
    Calculates the integral of r^2 * f(r).
    According to Poiseuille's law, Ks is proportional to this moment.
    """
    return float(np.trapz(r**2 * f, r))

        
def plot_vg_switched_axes(
    p0: lm2.SoilParams,
    pairs,
    *,
    # initial VG (for the blue reference only)
    theta_r: float,
    theta_s_initial: float,
    alpha0: float,
    n0: float,

    # porosity end-anchoring (only max-T hits this)
    phi_end_frac: float,

    title: str,
    save_png: str,

    # grids / ranges
    r_min_um: float = 0.02,
    r_max_um: float = 300.0,
    N_plot: int = 1000,
    h_fit_min: float = 0.1,
    h_fit_max: float = 1500.0,
    h_plot_max: float = 15000.0,

    # Eq.24 toggles
    xi_power_on_term2: bool = False,
    rho_s=None,              # optional: show BD in legend

    # NEW: adaptive residual water from PSD
    adaptive_theta_r: bool = True,
    h_pwp: float = 15000.0,  # cm; θ_r(T) := θ_from_PSD(h_pwp) with θ_r=0
    # optional PAW readout (only printed)
    h_fc: float | None = 330.0,  # cm (set None to skip)
    Ks0_ref: float | None = None,  # measured Ks for initial state [cm/day]
    BD_initial=None,
):
    """
    For each scenario (T, eps0):
      1) Compute PSD f(r,T).
      2) φ(T) = ∫ f dr.
      3) If adaptive_theta_r: θ_r(T) = θ_from_PSD(h_pwp) with θ_r=0 (conservative).
         Else: θ_r(T) = provided 'theta_r' (fixed).
      4) θ_s(T) = θ_r(T) + φ(T).
      5) Build θ(h) target from PSD (fit range) and fit VG(α,n) with θ_r, θ_s fixed.
      6) Plot θ(h) with x=θ, y=-Ψ (log); legend shows φ/φ0 (+ BD).
    Also computes the T=0 PSD→VG reference with the same θ_r derivation.
    Prints VG sets + moment-matched Kosugi equivalents. If h_fc is given, prints PAW.
    """
    import numpy as np, math, matplotlib.pyplot as plt

    # Grids
    r = np.geomspace(r_min_um, r_max_um, int(N_plot))
    h_fit  = np.geomspace(h_fit_min, h_fit_max, 200)
    h_plot = np.geomspace(1.0, h_plot_max, 400)

    # End-anchored mass rate so only max T hits phi_end_frac
    T_max = max(t for t, _ in pairs)
    M_over_u0w = 0.0 if abs(phi_end_frac-1.0) < 1e-12 else math.log(phi_end_frac)/T_max

    # Scenario labels
    Ts     = [t for t,_ in pairs]; epss = [e for _,e in pairs]
    T_min, T_max_seen = min(Ts), max(Ts)
    eps_min, eps_max  = min(epss), max(epss)
    def time_tag(T):
        return "Early" if abs(T - T_min) < 1e-9 else ("Late" if abs(T - T_max_seen) < 1e-9 else "Mid")
    def spread_tag(e0):
        return "Tight" if abs(e0 - eps_min) < 1e-9 else ("Broad" if abs(e0 - eps_max) < 1e-9 else "Mid")

    # Tables
    rows_vg = []  # (kind, T, eps0, theta_s, alpha, n, theta_r_used, Ks)
    rows_k  = []  # (kind, T, eps0, phi, sigma, r_m_um)
    rows_paw = [] # optional PAW

    # Initial VG (measured reference, blue)
    theta_ref = vg_theta(h_plot, theta_r, theta_s_initial, alpha0, n0)

    fig, ax = plt.subplots(figsize=(9.6, 6.0))
    #ax.plot(theta_ref, h_plot, lw=2.5, label="initial VG (reference)")
    rows_vg.append(("initial VG", None, None, theta_s_initial, alpha0, n0, theta_r, Ks0_ref, BD_initial))


    # Initial Kosugi from mapping p0
    rows_k.append(("initial K", None, None, p0.phi_0, p0.sigma, p0.r_m))

    # T=0 PSD → VG (apple-to-apple), with adaptive θr if enabled
    f0 = lm2.kosugi_psd(r, p0)
    phi0_psd = float(np.trapezoid(f0, r))
    
    # reference permeability proxy from initial PSD
    kappa0 = calc_permeability_proxy(r, f0)

    # T=0: Residualwasser aus GENAU diesem PSD (f0) ableiten
    if adaptive_theta_r:
        theta_pwp_cap0 = float(theta_from_psd(r, f0, np.array([h_pwp]), theta_r=0.0)[0])
        # nie unter den gemessenen Wert gehen
        theta_r0_used = max(theta_r, theta_pwp_cap0)
    else:
        theta_r0_used = theta_r

    # gesättigt bei T=0
    theta_s_fit0 = theta_r0_used + phi0_psd


    theta_target0 = theta_from_psd(r, f0, h_fit, theta_r=0.0) + theta_r0_used
    w0 = np.ones_like(h_fit); w0[(h_fit>=30) & (h_fit<=500)] = 1.5
    alpha0_fit, n0_fit = fit_vg_to_theta(h_fit, theta_target0, theta_r0_used, theta_s_fit0, weights=w0)
    theta_ref_psd = vg_theta(h_plot, theta_r0_used, theta_s_fit0, alpha0_fit, n0_fit)
    ax.plot(theta_ref_psd, h_plot, lw=2.5, color="k", label="Initiale Van-Genuchten Parameter (Messung)")
    
    Ks_T0 = Ks0_ref
    rows_vg.append(("T=0 PSD→VG", 0.0, 0.0, theta_s_fit0, alpha0_fit, n0_fit, theta_r0_used, Ks_T0, BD_initial))


    phi0_k, sigma0_k, rm0_k = kosugi_equiv_from_psd(r, f0)
    rows_k.append(("T=0 PSD→K", 0.0, 0.0, phi0_k, sigma0_k, rm0_k))

    if h_fc is not None:
        # PAW for the T=0 PSD→VG reference
        th_fc0 = float(vg_theta(np.array([h_fc]), theta_r0_used, theta_s_fit0, alpha0_fit, n0_fit)[0])
        th_pwp0 = float(vg_theta(np.array([h_pwp]), theta_r0_used, theta_s_fit0, alpha0_fit, n0_fit)[0])
        rows_paw.append(("T=0 PSD→VG", 0.0, 0.0, th_fc0 - th_pwp0))

    # vor der Schleife (einmal):
    T_min = min(t for t, _ in pairs)
    T_max = max(t for t, _ in pairs)
    eps_min = min(e for _, e in pairs)
    eps_max = max(e for _, e in pairs)

    # Scenario curves
    for (T, eps0) in pairs:
        f = lm2.analytical_solution_rdep(
            r, p0, float(T),
            eps0=float(eps0), M_over_u0w=float(M_over_u0w),
            r_min_um=float(r_min_um),
            xi_power_on_term2=bool(xi_power_on_term2),
        )
        phiT = float(np.trapezoid(f, r))

        # θr(T) adaptiv aus PSD bei h_pwp
        if adaptive_theta_r:
            theta_pwp_cap = float(theta_from_psd(r, f, np.array([h_pwp]), theta_r=0.0)[0])
            theta_r_T = theta_pwp_cap
        else:
            theta_r_T = theta_r

        theta_s_T = theta_r_T + phiT
        
        # --- scenario bulk density (for root model), if BD_initial (measured) is given ---
        BD_new = None
        if BD_initial is not None:
            BD_new = scale_bulk_density_physically_corrected(theta_s_initial, theta_s_T, BD_initial)


        # Zielkurve für VG-Fit
        theta_target = theta_from_psd(r, f, h_fit, theta_r=0.0) + theta_r_T
        w = np.ones_like(h_fit); w[(h_fit >= 30) & (h_fit <= 500)] = 1.5
        alpha_fit, n_fit = fit_vg_to_theta(
            h_fit, theta_target, theta_r_T, theta_s_T, weights=w
        )

        theta_curve = vg_theta(h_plot, theta_r_T, theta_s_T, alpha_fit, n_fit)

        # Porositätsverhältnis
        phi_frac = phiT / max(phi0_psd, 1e-12)
        bd_info  = f", BD≈{rho_s*(1.0 - phiT):.2f} g/cm³" if rho_s is not None else ""

        # Setzung
        if abs(T - T_min) < 1e-9:
            s_txt = "geringe Setzung"
        elif abs(T - T_max) < 1e-9:
            s_txt = "starke Setzung"
        else:
            s_txt = "mittlere Setzung"

        # Heterogenisierung
        if abs(eps0 - eps_min) < 1e-9:
            h_txt = "geringe Heterogenisierung"
        elif abs(eps0 - eps_max) < 1e-9:
            h_txt = "starke Heterogenisierung"
        else:
            h_txt = "mittlere Heterogenisierung"

        label = f"{s_txt}, {h_txt}, φ_end/φ₀={phi_frac:.2f}"
        ax.plot(theta_curve, h_plot, lw=2.0, ls="--", label=label)

        # --- Ks scaling from PSD (Poiseuille-based proxy): Ks ~ ∫ r² f(r) dr ---
        kappa_T = calc_permeability_proxy(r, f)
        ks_ratio = kappa_T / max(kappa0, 1e-30)
        Ks_T = Ks0_ref * ks_ratio

        rows_vg.append(("scenario", T, eps0, theta_s_T, alpha_fit, n_fit, theta_r_T, Ks_T, BD_new))


        phiK, sigmaK, rmK = kosugi_equiv_from_psd(r, f)
        rows_k.append(("scenario K", T, eps0, phiK, sigmaK, rmK))

        if h_fc is not None:
            th_fc  = float(vg_theta(np.array([h_fc ]), theta_r_T, theta_s_T, alpha_fit, n_fit)[0])
            th_pwp = float(vg_theta(np.array([h_pwp]), theta_r_T, theta_s_T, alpha_fit, n_fit)[0])
            rows_paw.append(("scenario", T, eps0, th_fc - th_pwp))

    # Axes/legend
    ax.set_xlabel("θ (cm³ cm⁻³)")
    ax.set_ylabel("-Ψ (cm)")
    ax.set_yscale("log")
    ax.set_ylim(1.0, h_plot_max)
    ax.set_xlim(0.0, max(theta_s_initial, float(np.max(theta_ref))))
    #ax.set_title(title + ("  (θᵣ from PSD@%g cm)" % h_pwp if adaptive_theta_r else "  (θᵣ fixed)"))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, ncol=1, loc="best")
    fig.tight_layout()
    if save_png:
        fig.savefig(save_png, dpi=160)

    # --- Print VG table
    print(f"\nVG parameter sets for: {title}")
    print(f"(theta_r mode = {'PSD@%gcm'%h_pwp if adaptive_theta_r else 'fixed'}, "
          f"phi_end/phi0 = {phi_end_frac:.2f}; r_max_um = {r_max_um:g})")
    print(f"{'type':<16} {'T':>5} {'eps0':>6} {'theta_r':>9} {'theta_s':>9} "
          f"{'alpha[1/cm]':>13} {'n':>6} {'Ks[cm/d]':>10} {'BD[g/cm³]':>10}")
    for kind, T, eps0, ths, a, nn, thr, Ks_val, BD_val in rows_vg:
        tstr  = f"{T:>5.2f}"    if T is not None else "  -- "
        eostr = f"{eps0:>6.2f}" if eps0 is not None else "  --  "
        print(f"{kind:<16} {tstr} {eostr} "
              f"{thr:>9.3f} {ths:>9.3f} {a:>13.6g} {nn:>6.3f} {Ks_val:>10.2f} {BD_val:>10.2f}")



    # --- Print equivalent Kosugi table
    print(f"\nEquivalent Kosugi (moment-matched) for: {title}")
    print(f"{'type':<16} {'T':>5} {'eps0':>6} {'phi':>8} {'sigma':>8} {'r_m [um]':>11}")
    for kind, T, eps0, phiK, sigK, rmK in rows_k:
        tstr  = f"{T:>5.2f}"    if T is not None else "  -- "
        eostr = f"{eps0:>6.2f}" if eps0 is not None else "  --  "
        print(f"{kind:<16} {tstr} {eostr} {phiK:>8.3f} {sigK:>8.3f} {rmK:>11.3f}")

    # --- Optional PAW table
    # --- Optional PAW table (jetzt MIT Sand-Referenz) ---
    if h_fc is not None:
        # 1) Sand-Referenz zuerst ausgeben
        theta_r_sand = 0.046
        theta_s_sand = 0.368
        alpha_sand   = 0.0152
        n_sand       = 2.736
        theta_fc_sand  = float(vg_theta(np.array([h_fc]),  theta_r_sand, theta_s_sand, alpha_sand, n_sand)[0])
        theta_pwp_sand = float(vg_theta(np.array([h_pwp]), theta_r_sand, theta_s_sand, alpha_sand, n_sand)[0])
        paw_sand = max(0.0, theta_fc_sand - theta_pwp_sand)

        print(f"\nPAW (θ_FC - θ_PWP) with FC={h_fc} cm and PWP={h_pwp} cm for: {title}")
        print(f"{'type':<16} {'T':>5} {'eps0':>6} {'PAW':>8}")
        # Sand zuerst
        print(f"{'Sand (unbehandelt)':<16} {'  -- ':>5} {'  --  ':>6} {paw_sand:>8.3f}")
        # dann alle Szenarien wie bisher
        for kind, T, eps0, paw in rows_paw:
            tstr  = f"{T:>5.2f}"    if T is not None else "  -- "
            eostr = f"{eps0:>6.2f}" if eps0 is not None else "  --  "
            print(f"{kind:<16} {tstr} {eostr} {paw:>8.3f}")

    # --- Save VG results to JSON with simple scenario IDs ---
    vg_out = []
    scenario_id = 1

    for kind, T, eps0, ths, a, nn, thr, Ks_val, BD_val in rows_vg:
        # ID vergeben:
        if kind == "scenario":
            id_ = scenario_id      # Szenarien: 1..6
            scenario_id += 1
        elif kind == "initial VG":
            id_ = 7                # gemessener Datensatz
            kind = "measured VG"   # schönerer Name
        else:
            id_ = 0                # z.B. "T=0 PSD→VG" -> brauchst du später nicht

        vg_out.append({
            "id": id_,
            "kind": kind,
            "T": T,
            "eps0": eps0,
            "theta_r": thr,
            "theta_s": ths,
            "alpha": a,
            "n": nn,
            "Ks": Ks_val,
            "BD": BD_val,
        })

    with open("VG_params_af3_from_scenario_analysis.json", "w") as f:
        json.dump(vg_out, f, indent=2)

    print("\nJSON file written: VG_parameters_from_scenario_analysis.json\n")


def scale_bulk_density_physically_corrected(theta_s_initial, theta_s_new, BD_initial, rho_assumed=1.4):
    """
    Estimate scenario-specific Bulk Density (BD) from hydraulic porosity shrinkage.
    
    IMPORTANT:
    - This BD scaling is NOT used inside the Leij model itself.
    - It is only a post-processing helper: we use the BD_new values
      as input for the *root-growth* model (penetration resistance).
    - Leij computes hydraulic compaction (θs changes), and we convert
      that hydraulic shrinkage into a reasonable structural BD estimate
      for the other script (Soil3_Kompost_Baseline.py).

    Reason:
    Using BD_initial together with θs_initial directly would imply an
    unrealistically high particle density. We instead anchor the structural
    density to a realistic compost particle density (≈1.4 g/cm³).
    """

    # 1. Structural porosity at initial state (from measured BD)
    phi_struct_initial = 1.0 - (BD_initial / rho_assumed)

    # 2. Shrinkage factor from the hydraulic model (Leij PSD collapse)
    shrinkage = theta_s_new / theta_s_initial

    # 3. New structural porosity after compaction
    phi_struct_new = phi_struct_initial * shrinkage

    # 4. Convert back to BD
    return (1.0 - phi_struct_new) * rho_assumed



if __name__ == "__main__":
    # --- measured compost+sand (t0) ---
    theta_s, theta_r = 0.707, 0.001
    alpha_cm_inv, n = 0.0452, 1.41
    Ks_measured = 538.0
    BD_initial = 0.84        # measured compost BD
    p0 = lm2.SoilParams.from_van_genuchten(theta_s, theta_r, alpha_cm_inv, n)
    print("p0:", p0)

    # --- fixed particle density ---
    rho_s = 2.65  # g/cm³

    # --- five representative (T, eps0) pairs ---
    pairs = [
        (0.3, 0.08),  # Setzung: gering; Heterogenisierung: gering
        (0.3, 0.25),  # Setzung: gering; Heterogenisierung: stark
        (1.0, 0.08),  # Setzung: gering; Heterogenisierung: gering
        (1.0, 0.25),  # Setzung: stark; Heterogenisierung: stark
        (2.0, 0.08),  # Setzung: stark; Heterogenisierung: gering
        (2.0, 0.25),  # Setzung: stark; Heterogenisierung: stark
    ]
    # --- common grids ---
    r_min_um_psd = 0.001
    r_max_um_all = 3000.0
    N_plot_psd   = 900
    N_plot_vg    = 1000

    # ---------- single BD-anchored case: from 0.84 -> 1.64 ----------
    BD_env = 1.64           # umliegender Sand
    phi0   = p0.phi_0
    phi_env = 1.0 - BD_env / rho_s
    phi_drop_anchor = max(1e-6, min(1.0, phi_env / phi0))
    print(f"Using fixed rho_s={rho_s:.2f} g/cm³, BD={BD_env:.2f} g/cm³ -> φ_end/φ0 ≈ {phi_drop_anchor:.2f}")

    plot_psd_and_vg_side_by_side(
        p0, pairs,
        rho_s=rho_s,
        phi_end_frac=phi_drop_anchor,
        theta_r=theta_r,
        theta_s_initial=theta_s,
        alpha0=alpha_cm_inv,
        n0=n,
        save_png="combined_psd_vg_bd_1_64.png",
    )
    plot_vg_switched_axes(
        p0, pairs,
        theta_r=theta_r,
        theta_s_initial=theta_s,
        alpha0=alpha_cm_inv,
        n0=n,
        phi_end_frac=phi_drop_anchor,
        title="Water retention (VG fits) — current run",
        save_png="combined_psd_vg_current.png",
        rho_s=rho_s,
        h_fc=63.0,
        r_min_um=0.001,     # wie im Plot
        r_max_um=3000.0,    # wie im Plot
        N_plot=1000,
        Ks0_ref=Ks_measured,
        BD_initial=BD_initial,
    )