"""
Morris-Design (Elementary Effects) mit Trajektorien-Auswahl nach Campolongo et al. (2007)
und Auswertung (mu*, mu, sigma) inkl. Bootstrap-Konfidenzintervallen.

Bewusst ohne SALib-Abhängigkeit implementiert, damit Design und Auswertung exakt
zusammenpassen (gemeinsame Seeds pro Trajektorie, Mittelung über Replikate).
Die Ergebnisse entsprechen SALib.analyze.morris (mu_star = mean |EE|).
"""
import numpy as np


def _one_trajectory(k, p, rng):
    delta = p / (2 * (p - 1))
    grid = np.arange(p) / (p - 1)
    base_levels = grid[grid <= 1 - delta + 1e-12]
    x = rng.choice(base_levels, size=k)
    sign = rng.choice([-1, 1], size=k)
    x = np.where(sign < 0, x + delta, x)          # Start oben, falls Richtung negativ
    order = rng.permutation(k)
    traj = [x.copy()]
    for i in order:
        x = x.copy()
        x[i] += sign[i] * delta
        traj.append(x)
    return np.array(traj), order


def _traj_distance(a, b):
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    return d.sum()


def sample(k, r, p=4, n_candidates=300, seed=0):
    """Liefert (X, order) mit X: (r, k+1, k) in [0,1] und order: (r, k)."""
    rng = np.random.default_rng(seed)
    cands = [_one_trajectory(k, p, rng) for _ in range(n_candidates)]
    T = [c[0] for c in cands]
    D = np.zeros((n_candidates, n_candidates))
    for a in range(n_candidates):
        for b in range(a + 1, n_candidates):
            D[a, b] = D[b, a] = _traj_distance(T[a], T[b])
    # Greedy-Auswahl: zuerst das Paar mit größtem Abstand, dann jeweils maximaler Zuwachs
    a, b = np.unravel_index(D.argmax(), D.shape)
    chosen = [a, b]
    while len(chosen) < r:
        rest = [i for i in range(n_candidates) if i not in chosen]
        gain = [np.sqrt((D[i, chosen] ** 2).sum()) for i in rest]
        chosen.append(rest[int(np.argmax(gain))])
    X = np.stack([T[i] for i in chosen])
    order = np.stack([cands[i][1] for i in chosen])
    return X, order


def elementary_effects(X, Y, order):
    """X: (r, k+1, k) normiert, Y: (r, k+1) Output. Gibt EE: (r, k) zurück.
    EE ist auf den normierten Faktorbereich bezogen (Änderung des Outputs pro
    Durchlaufen von Δ des jeweiligen Bereichs) -> zwischen Faktoren vergleichbar."""
    r, kp1, k = X.shape
    EE = np.full((r, k), np.nan)
    for t in range(r):
        for step, i in enumerate(order[t]):
            dx = X[t, step + 1, i] - X[t, step, i]
            EE[t, i] = (Y[t, step + 1] - Y[t, step]) / dx
    return EE


def analyze(EE, n_boot=1000, ci=0.95, seed=0):
    """mu*, mu, sigma und Bootstrap-KI für mu* über die Trajektorien."""
    rng = np.random.default_rng(seed)
    r, k = EE.shape
    out = {
        "mu_star": np.nanmean(np.abs(EE), 0),
        "mu": np.nanmean(EE, 0),
        "sigma": np.nanstd(EE, 0, ddof=1),
        "n_valid": np.sum(~np.isnan(EE), 0),
    }
    boots = np.array([np.nanmean(np.abs(EE[rng.integers(0, r, r)]), 0) for _ in range(n_boot)])
    lo, hi = np.nanpercentile(boots, [(1 - ci) / 2 * 100, (1 + ci) / 2 * 100], axis=0)
    out["mu_star_lo"], out["mu_star_hi"] = lo, hi
    return out
