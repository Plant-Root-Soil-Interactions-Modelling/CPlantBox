"""
Toy model ONLY for testing the pipeline. Biologically meaningless.

Reads the same XML as the real simulation and produces a travelling wave over the rings of the
Petri dish in the format of model_hook.simulate (raw.npz): filament density as a sigmoid, tips as a
bell around the arrival time, anastomoses ~ C rho tips, with Poisson noise. Front speed,
amplitude, width and BAS share depend arbitrarily on the factors.
`b` of the runners is deliberately NOT used (b_prob = 0) -> control factor.
"""
import xml.etree.ElementTree as ET

import numpy as np

import config

FWHM_BELL = 2 * np.log(3 + 2 * np.sqrt(2))


def _read(xml_path):
    root = ET.parse(xml_path).getroot()
    out = {}
    for hy in root.iter("hyphae"):
        out[int(hy.get("subType"))] = {p.get("name"): (float(p.get("value")), float(p.get("dev") or 0))
                                       for p in hy.iter("parameter") if p.get("value") is not None}
    return out


def ring_geometry():
    """Outer ring radii [cm] and ring areas [mm²] as amfs.ringLocations / amfs.ringVolume (without plantbox)."""
    radius, barrier_thickness = config.RING_RADIUS, 0.16   # rings up to RING_RADIUS around the crossing point (here: 0)
    outer = np.array([radius * np.sqrt(i / config.N_RINGS) for i in range(1, config.N_RINGS + 1)])
    c = barrier_thickness / 2
    R = np.concatenate(([0.], outer))
    c_ = np.minimum(c, R)
    area = np.pi / 2 * R * R + c_ * np.sqrt(R * R - c_ * c_) + R * R * np.arcsin(np.divide(c_, R, out=np.zeros_like(R), where=R > 0))
    return outer, np.diff(area) * 100


def influx(rng, time_star, hours, roots_on, scale=1.):
    """Supply diagnostics as model_hook.influx_diagnostics: rises until t*, then decays with the roots turned off."""
    edges = np.arange(0., time_star + hours / 24 + 1 / 24 + 1e-9, 1 / 24)
    h = (edges[:-1] - time_star) * 24
    rate = 20 * scale * np.clip(1 + h / 48, 0, 1) * np.where((h > 0) & ~roots_on, np.exp(-np.maximum(h, 0) / 40), 1.)
    cross = rng.poisson(rate).astype(float)
    return {"diag_edges": edges, "diag_cross_in": cross, "diag_cross_out": rng.poisson(0.1 * rate).astype(float),
            "diag_cross_in_runner": rng.binomial(cross.astype(int), 0.7).astype(float),
            "diag_cross_in_bas": np.zeros_like(cross), "diag_root_hyphae_new": rng.poisson(2 * rate).astype(float),
            "diag_root_length_new": np.where((h > 0) & ~roots_on, 0., 0.1 * (1 + rng.normal(0, 0.1, len(h))))}


def simulate(xml_path, seed, outdir, hours=None, roots_on=False, root_days=None, sim_args=None, layer=None):
    hours = hours or config.HOURS_HYPHAE
    root_days = root_days or config.ROOT_DAYS
    sim_args = sim_args or {}
    emerge = sim_args.get("hyphal_emergence_density", 4.0) / 4.0
    dist = (sim_args.get("seed_pos_x", -4.2) + 0.26) / -4.46      # 1 = 5 mm from the edge, small = close to the barrier
    rng = np.random.default_rng(seed)
    P = _read(xml_path)
    R, B = P[config.RUNNER], P[config.BAS]
    val = lambda prm, k: prm[k][0]
    cv = R["v"][1] / R["v"][0]
    dx = R["dx"][0]
    noise = lambda s: np.exp(s * rng.normal())

    # front speed [mm/h], amplitude of the tip wave [mm^-2], width [h], saturation [mm mm^-2]
    vf = 0.6 * (val(R, "v") / 0.576) ** 0.7 * (val(B, "v") / 0.576) ** 0.2 * (val(R, "hlt") / 8) ** 0.3 \
        * (val(R, "lmax") / 3) ** 0.1 * np.cos(val(R, "theta") - 1.047) * noise(0.05)
    if vf < 0.42:
        from model_hook import NoTstar
        crossed = 10 * max(vf - 0.38, 0) / 0.04          # very slow: ≤ 1 cm -> excluded, otherwise failed
        raise NoTstar(f"status not_crossed (toy model: front too slow, {crossed:.2f} cm)", crossed)
    K2 = 0.15 * (val(B, "b") / 0.96) ** 0.5 * (0.1 / val(R, "ln")) ** 0.4 * noise(0.05)
    width = 10 * (1 + cv) * (val(B, "hlt") / 2) ** 0.2 * noise(0.05)
    K1 = 1.4 * (val(B, "hlt") / 2) ** 0.5 * (val(R, "lmax") / 3) ** 0.3 * noise(0.05)
    frac_bas = np.clip(0.3 * (val(B, "hlt") / 2) ** 0.4 * (val(B, "b") / 0.96) ** 0.3
                       * (val(B, "v") / val(R, "v")) ** 0.2 * (1 + 0.3 * np.sin(val(B, "theta"))) * noise(0.05), 0, 1)
    C = 2e-5 * (val(R, "distTH") / 0.01) ** 2 * val(R, "ana") * (1 + 20 * (dx - 0.01)) * (config.HYPHAE_LAYER / (layer or config.HYPHAE_LAYER))

    outer, area = ring_geometry()
    inner = np.concatenate(([0.], outer[:-1]))
    r = np.sqrt((inner * inner + outer * outer) / 2) * 10           # ring centres [mm]
    wait_days = 2 / val(R, "v") * dist ** 0.8 * (5 / root_days) ** 0.3 / emerge ** 0.2 * noise(0.1)
    time_star = root_days + wait_days
    t_h = np.arange(0, hours + 1)[1:].astype(float)   # columns as getParasPerRing: times[1:]
    t_arr = (r[:, None] - r[0]) / vf                                 # arrival per ring [h]
    lam = FWHM_BELL / width
    e = np.exp(-np.abs(lam * (t_h[None, :] - t_arr)))
    tip_dens = 4 * K2 * e / (1 + e) ** 2
    rho = K1 / (1 + np.exp(np.minimum(lam * (t_arr - t_h[None, :]), 700.)))
    tips = rng.poisson(tip_dens * area[:, None]).astype(float)
    tips_bas = rng.binomial(tips.astype(int), frac_bas).astype(float)
    rho = rho * noise(0.02)
    ana_rate = C * rho * 1000 * tips                                 # per hour
    ana_bas = rng.poisson(ana_rate * frac_bas).astype(float)
    ana_runner = rng.poisson(ana_rate * (1 - frac_bas)).astype(float)
    branch_rate = (val(B, "b") * frac_bas + 0.576 / val(R, "ln") * 0.05 * (1 - frac_bas)) / 24   # per tip and hour
    raw = {"seed": seed, "time_star": time_star, "location": outer, "area": area,
           "times": time_star + t_h / 24, "tip_dens": tips / area[:, None], "active_tips": tips,
           "branch_new": rng.poisson(branch_rate * tips).astype(float),
           "branch_new_runner": rng.poisson(0.576 / val(R, "ln") * 0.25 / 24 * (tips - tips_bas)).astype(float),
           "angle_rh_mean": float(np.degrees(val(R, "theta")) * 1.15 * noise(0.03)),
           "angle_bas_mean": float(np.degrees(val(B, "theta")) * 1.2 * noise(0.03)),
           "persistence_length_median": float(15 * (val(R, "tropismN") / 4) ** 0.8 * noise(0.1)),
           # tip speeds per hour, as model_hook.tip_speed_diagnostics [µm/h]
           "diag_speed_hours": np.arange(1, hours + 1, dtype=float),
           "diag_front_tip_speed": 240 * val(R, "v") / 0.576 * (1 + 1.2 * cv) * noise(0.03) * np.ones(int(hours)),
           "diag_tip_speed_all": 240 * val(R, "v") / 0.576 * noise(0.03) * np.ones(int(hours)),
           # share of the tips beyond the rings: the front (r ≈ vf t) passes RING_RADIUS
           "diag_frac_tips_outside_rings": np.clip((r[0] + vf * np.arange(1, hours + 1) - config.RING_RADIUS * 10) / 10, 0, 1),
           # image graph at IMAGE_HOURS and the end, as model_hook
           **{f"{k}_{h:g}h": v for h in sorted({x for x in config.IMAGE_HOURS + config.EXTRA_HOURS if x <= hours} | {hours})
              for k, v in {"bas_image_share": float(np.clip(frac_bas * 1.2 * (0.1 / val(R, "ln")) ** 0.1 * (0.8 + 0.2 * h / 60), 0, 1)),
                           "bas_image_share_subtype": float(frac_bas),
                           "bas_image_length_xy_um": float((rho[:, min(int(h), rho.shape[1]) - 1] * area).sum() * 1000),
                           "root_eff_E_rel": 0.5 * val(R, "ana") * (val(R, "distTH") / 0.01) ** 0.3 * noise(0.05),
                           "root_eff_E_g_rel": 0.4 * val(R, "ana") * (h / 60) ** 0.5 * noise(0.05),
                           "root_eff_C_rel": 0.1 * (0.1 / val(R, "ln")) ** 0.3 * noise(0.05)}.items()},
           "root_eff_E_rel": 0.5 * val(R, "ana") * (val(R, "distTH") / 0.01) ** 0.3 * (val(B, "a") / 0.0025) ** 0.2 * noise(0.05),
           "ana_new": ana_bas + ana_runner, "length_densities": rho, "wait_days": wait_days,
           "length_runner": rho * (1 - frac_bas), "length_bas": rho * frac_bas,
           "active_tips_runner": tips - tips_bas, "active_tips_bas": tips_bas,
           "ana_new_runner": ana_runner, "ana_new_bas": ana_bas, "roots_on": roots_on,
           "ring_center_y": 0., "n_crossings_tstar": int(rng.integers(1, 9)), "ring_radius": config.RING_RADIUS,
           **influx(rng, time_star, hours, roots_on, emerge)}
    path = outdir / "raw.npz"
    np.savez_compressed(path, **raw)
    return [path]
