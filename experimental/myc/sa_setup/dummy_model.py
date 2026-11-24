"""
Spielzeugmodell NUR zum Testen der Pipeline. Biologisch bedeutungslos.

Liest dieselbe XML wie die echte Simulation und erzeugt eine laufende Welle über die Ringe der
Petrischale im Format von model_hook.simulate (raw.npz): Filamentdichte als Sigmoid, Spitzen als
Glocke um die Ankunftszeit, Anastomosen ~ C rho tips, mit Poisson-Rauschen. Frontgeschwindigkeit,
Amplitude, Breite und BAS-Anteil hängen willkürlich von den Faktoren ab.
`b` der Runner wird bewusst NICHT verwendet (b_prob = 0) -> Kontrollfaktor.
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
    """Äußere Ringradien [cm] und Ringflächen [mm²] wie amfs.ringLocations / amfs.ringVolume (ohne plantbox)."""
    radius, barrier_thickness = 9.4 / 2, 0.16
    outer = np.array([radius * np.sqrt(i / config.N_RINGS) for i in range(1, config.N_RINGS + 1)])
    c = barrier_thickness / 2
    R = np.concatenate(([0.], outer))
    c_ = np.minimum(c, R)
    area = np.pi / 2 * R * R + c_ * np.sqrt(R * R - c_ * c_) + R * R * np.arcsin(np.divide(c_, R, out=np.zeros_like(R), where=R > 0))
    return outer, np.diff(area) * 100


def simulate(xml_path, seed, outdir, hours=None):
    hours = hours or config.HOURS_HYPHAE
    rng = np.random.default_rng(seed)
    P = _read(xml_path)
    R, B = P[config.RUNNER], P[config.BAS]
    val = lambda prm, k: prm[k][0]
    cv = R["v"][1] / R["v"][0]
    dx = R["dx"][0]
    noise = lambda s: np.exp(s * rng.normal())

    # Frontgeschwindigkeit [mm/h], Amplitude der Spitzenwelle [mm^-2], Breite [h], Sättigung [mm mm^-2]
    vf = 0.6 * (val(R, "v") / 0.576) ** 0.7 * (val(B, "v") / 0.576) ** 0.2 * (val(R, "hlt") / 8) ** 0.3 \
        * (val(R, "lmax") / 3) ** 0.1 * np.cos(val(R, "theta") - 1.047) * noise(0.05)
    if vf < 0.42:
        raise RuntimeError("status not_crossed (Spielzeugmodell: Front zu langsam)")
    K2 = 0.15 * (val(B, "b") / 0.96) ** 0.5 * (0.1 / val(R, "ln")) ** 0.4 * noise(0.05)
    width = 10 * (1 + cv) * (val(B, "hlt") / 2) ** 0.2 * noise(0.05)
    K1 = 1.4 * (val(B, "hlt") / 2) ** 0.5 * (val(R, "lmax") / 3) ** 0.3 * noise(0.05)
    frac_bas = np.clip(0.3 * (val(B, "hlt") / 2) ** 0.4 * (val(B, "b") / 0.96) ** 0.3
                       * (val(B, "v") / val(R, "v")) ** 0.2 * (1 + 0.3 * np.sin(val(B, "theta"))) * noise(0.05), 0, 1)
    C = 2e-5 * (val(R, "distTH") / 0.01) ** 2 * val(R, "ana") * (1 + 20 * (dx - 0.01))

    outer, area = ring_geometry()
    inner = np.concatenate(([0.], outer[:-1]))
    r = np.sqrt((inner * inner + outer * outer) / 2) * 10           # Ringmitten [mm]
    wait_days = 2 / val(R, "v") * noise(0.1)
    time_star = config.ROOT_DAYS + wait_days
    t_h = np.arange(0, hours + 1)[1:].astype(float)   # Spalten wie getParasPerRing: times[1:]
    t_arr = (r[:, None] - r[0]) / vf                                 # Ankunft je Ring [h]
    lam = FWHM_BELL / width
    e = np.exp(-np.abs(lam * (t_h[None, :] - t_arr)))
    tip_dens = 4 * K2 * e / (1 + e) ** 2
    rho = K1 / (1 + np.exp(np.minimum(lam * (t_arr - t_h[None, :]), 700.)))
    tips = rng.poisson(tip_dens * area[:, None]).astype(float)
    tips_bas = rng.binomial(tips.astype(int), frac_bas).astype(float)
    rho = rho * noise(0.02)
    ana_rate = C * rho * 1000 * tips                                 # pro Stunde
    ana_bas = rng.poisson(ana_rate * frac_bas).astype(float)
    ana_runner = rng.poisson(ana_rate * (1 - frac_bas)).astype(float)
    branch_rate = (val(B, "b") * frac_bas + 0.576 / val(R, "ln") * 0.05 * (1 - frac_bas)) / 24   # pro Spitze und Stunde
    raw = {"seed": seed, "time_star": time_star, "location": outer, "area": area,
           "times": time_star + t_h / 24, "tip_dens": tips / area[:, None], "active_tips": tips,
           "branch_new": rng.poisson(branch_rate * tips).astype(float),
           "branch_new_runner": rng.poisson(0.576 / val(R, "ln") * 0.25 / 24 * (tips - tips_bas)).astype(float),
           "root_eff_E_rel": 0.5 * val(R, "ana") * (val(R, "distTH") / 0.01) ** 0.3 * (val(B, "a") / 0.0025) ** 0.2 * noise(0.05),
           "ana_new": ana_bas + ana_runner, "length_densities": rho, "wait_days": wait_days,
           "length_runner": rho * (1 - frac_bas), "length_bas": rho * frac_bas,
           "active_tips_runner": tips - tips_bas, "active_tips_bas": tips_bas,
           "ana_new_runner": ana_runner, "ana_new_bas": ana_bas}
    path = outdir / "raw.npz"
    np.savez_compressed(path, **raw)
    return [path]
