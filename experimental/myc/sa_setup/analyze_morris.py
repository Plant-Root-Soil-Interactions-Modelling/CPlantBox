"""
Phase 1: Morris-Auswertung für alle Outputs, Teilbereiche und Zeitpunkte.

    python analyze_morris.py
      -> results/morris_indices.csv     mu*, KI, mu, sigma je (time, region, output, factor)
      -> results/plots/morris_ranking_t<T>.png    mu* mit 95%-KI, gesamter Zielbereich
      -> results/plots/morris_musigma_t<T>.png    mu* vs. sigma (Nichtlinearität/Interaktion)
      -> results/plots/morris_regions_t<T>.png    mu* je Teilbereich (zeilenweise normiert)

mu* ist in Output-Einheiten pro Durchlaufen des gesamten Faktorbereichs angegeben
(normierte Faktoren) und damit zwischen Faktoren, nicht aber zwischen Outputs vergleichbar.
Für den Vergleich dient mu*_rel = mu* / max(mu*) innerhalb eines Outputs.
mu*, mu, sigma und das Bootstrap-KI kommen aus SALib.analyze.morris; das Rauschniveau
aus den Replikaten wird hier berechnet.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from SALib.analyze import morris as morris_analyze

import config
from design import morris_problem

KEY = ["time", "region", "output"]
NAMES = [f.name for f in config.FACTORS]
ROLE = {f.name: f.role for f in config.FACTORS}
SUBTYPE = {f.name: f.subtype for f in config.FACTORS}

# Farben (validiert: Kategorien Puller/Runner, Sequenz für Heatmap)
C_PULLER, C_RUNNER, C_CTRL = "#2a78d6", "#eb6834", "#8a8a85"
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
SEQ = ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def _color(name):
    if ROLE[name] == "control":
        return C_CTRL
    return C_PULLER if SUBTYPE[name] == 1 else C_RUNNER


def _elementary_effects(X, Y):
    """EEs je Trajektorie (nur fürs Rauschniveau, SALib gibt sie nicht heraus).
    X: (r, k+1, k) normiert, Y: (..., r, k+1). Gibt EE: (..., r, k) zurück, normiert wie
    SALib: Output-Änderung pro Schritt Δ im normierten Faktorbereich."""
    dX = np.diff(X, axis=1)                              # (r, k, k): je Schritt ändert sich ein Faktor
    i = np.abs(dX).argmax(-1)                            # (r, k): welcher Faktor in welchem Schritt
    dx = np.take_along_axis(dX, i[..., None], -1)[..., 0]
    ee_step = np.diff(Y, axis=-1) / dx                   # (..., r, k) in Schrittreihenfolge
    EE = np.empty_like(ee_step)
    np.put_along_axis(EE, np.broadcast_to(i, ee_step.shape), ee_step, -1)
    return EE


def compute():
    design = pd.read_csv(config.RESULTS_DIR / "design_morris.csv")
    problem = morris_problem()
    k = problem["num_vars"]
    X = np.load(config.RESULTS_DIR / "morris_X.npy").reshape(-1, k + 1, k)   # (r, k+1, k)
    out = pd.read_csv(config.RESULTS_DIR / "outputs_morris.csv").merge(
        design[["run_id", "traj", "step", "rep"]], on="run_id")
    r, kp1, _ = X.shape
    R = int(design["rep"].max()) + 1
    vals = out.set_index(KEY + ["rep", "traj", "step"])["value"]
    rows = []
    for key, s in vals.groupby(level=[0, 1, 2]):
        Y = np.full((R, r, kp1), np.nan)
        for (_, _, _, rep, t, st), v in s.items():
            Y[rep, t, st] = v
        # Jedes Replikat ist eine in sich konsistente Trajektorie (gemeinsamer Seed) -> zählt nur
        # vollständig. Gültig sind Trajektorien mit mind. einem vollständigen Replikat.
        complete = np.isfinite(Y).all(-1)                  # (R, r)
        Y = np.where(complete[..., None], Y, np.nan)
        valid = complete.any(0)
        n_valid = int(valid.sum())
        res = {m: np.full(k, np.nan) for m in ("mu_star", "mu", "sigma")}
        res["n_valid"] = np.full(k, n_valid)
        res["mu_star_lo"], res["mu_star_hi"] = np.full(k, np.nan), np.full(k, np.nan)
        if n_valid:
            # EE ist linear in Y: Replikate vorher mitteln == EEs über Replikate mitteln
            Y_mean = np.nanmean(Y[:, valid], 0)
            Si = morris_analyze.analyze(problem, X[valid].reshape(-1, k), Y_mean.reshape(-1),
                                        num_resamples=1000, conf_level=0.95,
                                        num_levels=config.MORRIS_LEVELS, seed=config.BASE_SEED)
            for m in ("mu_star", "mu", "sigma"):
                res[m] = np.asarray(Si[m], dtype=float)
            # SALib liefert das KI als Halbbreite (Normalapproximation des Bootstraps)
            res["mu_star_lo"] = np.maximum(res["mu_star"] - Si["mu_star_conf"], 0)
            res["mu_star_hi"] = res["mu_star"] + Si["mu_star_conf"]
        # Rauschniveau je Faktor: mu*, das ein wirkungsloser Faktor allein durch Seed-Rauschen
        # hätte. SE des Replikat-Mittels je EE, E|N(0,se)| = se*sqrt(2/pi), gemittelt über Trajektorien.
        if R > 1:
            EE_rep = _elementary_effects(X, Y)             # (R, r, k)
            se = np.nanstd(EE_rep, 0, ddof=1) / np.sqrt(R)
            res["noise_level"] = np.nanmean(se, 0) * np.sqrt(2 / np.pi)
        else:
            res["noise_level"] = np.full(k, np.nan)
        for j, n in enumerate(NAMES):
            rows.append(dict(zip(KEY, key), factor=n, role=ROLE[n],
                             **{m: float(res[m][j]) for m in res}))
    ind = pd.DataFrame(rows)
    mx = ind.groupby(KEY)["mu_star"].transform("max").replace(0, np.nan)
    ind["mu_star_rel"] = ind["mu_star"] / mx
    # "im Rauschen": mu* nicht deutlich (Faktor 2) über dem, was Seed-Rauschen allein erzeugt.
    # Hinweis: Der Kontrollfaktor liefert mit gemeinsamen Seeds exakt 0, wenn er wirklich
    # wirkungslos ist — er prüft den b_prob-Schalter, nicht das Rauschniveau.
    ind["in_noise"] = ind["mu_star"] <= 2 * ind["noise_level"]
    ind.to_csv(config.RESULTS_DIR / "morris_indices.csv", index=False)
    return ind


def plot_ranking(ind, t):
    outs = list(ind.output.unique())
    sub = ind[(ind.time == t) & (ind.region == "all")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 9), sharey=False)
    for ax, o in zip(axes.flat, outs):
        d = sub[sub.output == o].sort_values("mu_star")
        y = np.arange(len(d))
        ax.barh(y, d.mu_star, color=[_color(n) for n in d.factor], height=0.6, edgecolor="white", lw=1)
        ax.errorbar(d.mu_star, y, xerr=[d.mu_star - d.mu_star_lo, d.mu_star_hi - d.mu_star],
                    fmt="none", ecolor=INK2, elinewidth=0.8, capsize=2)
        if d.noise_level.notna().any():   # Rauschniveau je Faktor als kurze Markierung
            ax.scatter(2 * d.noise_level, y, marker="|", s=90, color=INK, lw=1.5, zorder=4)
        ax.set_yticks(y, d.factor, fontsize=8, color=INK)
        ax.set_title(o, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("μ*", fontsize=8, color=INK2)
        _style(ax)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig.legend(handles=[Patch(color=C_PULLER, label="Puller"), Patch(color=C_RUNNER, label="Runner"),
                        Patch(color=C_CTRL, label="Kontrolle"),
                        Line2D([], [], marker="|", ls="", ms=10, mew=1.5, color=INK,
                               label="2 × Rauschniveau (aus Replikaten)")],
               loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle(f"Morris μ* mit 95%-Bootstrap-KI — gesamter Zielbereich, t = {t:g} d",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(config.RESULTS_DIR / "plots" / f"morris_ranking_t{t:g}.png", dpi=150)
    plt.close(fig)


def plot_musigma(ind, t, n_label=6):
    outs = list(ind.output.unique())
    sub = ind[(ind.time == t) & (ind.region == "all")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))
    for ax, o in zip(axes.flat, outs):
        d = sub[sub.output == o]
        for n, row in d.set_index("factor").iterrows():
            mk = "X" if ROLE[n] == "control" else ("o" if SUBTYPE[n] == 1 else "s")
            ax.scatter(row.mu_star, row.sigma, s=40, marker=mk, color=_color(n),
                       edgecolor="white", lw=1, zorder=3)
        lim = max(d.mu_star.max(), d.sigma.max()) * 1.1
        ax.plot([0, lim], [0, lim], color=GRID, lw=1, zorder=1)        # σ = μ*
        for _, row in d.nlargest(n_label, "mu_star").iterrows():
            ax.annotate(row.factor, (row.mu_star, row.sigma), xytext=(4, 3),
                        textcoords="offset points", fontsize=7, color=INK)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_title(o, fontsize=10, color=INK, loc="left")
        ax.set_xlabel("μ*", fontsize=8, color=INK2); ax.set_ylabel("σ", fontsize=8, color=INK2)
        _style(ax); ax.grid(axis="y", color=GRID, lw=0.6)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", color=C_PULLER, label="Puller"),
                        Line2D([], [], marker="s", ls="", color=C_RUNNER, label="Runner"),
                        Line2D([], [], marker="X", ls="", color=C_CTRL, label="Kontrolle"),
                        Line2D([], [], color=GRID, label="σ = μ*")],
               loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle(f"Morris μ* vs. σ (σ groß: nichtlinear oder Interaktionen) — t = {t:g} d",
                 x=0.01, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(config.RESULTS_DIR / "plots" / f"morris_musigma_t{t:g}.png", dpi=150)
    plt.close(fig)


def plot_regions(ind, t):
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("seq", SEQ)
    outs = list(ind.output.unique())
    regions = list(config.REGIONS)
    sub = ind[(ind.time == t) & (ind.region.isin(regions))]
    order = (ind[(ind.time == t) & (ind.region == "all")]
             .groupby("factor").mu_star_rel.mean().sort_values(ascending=False).index)
    fig, axes = plt.subplots(1, len(outs), figsize=(2.1 * len(outs) + 1.5, 7), sharey=True)
    for ax, o in zip(np.atleast_1d(axes), outs):
        M = (sub[sub.output == o].pivot(index="factor", columns="region", values="mu_star_rel")
             .reindex(index=order, columns=regions))
        im = ax.imshow(M.values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(regions)), regions, fontsize=8, color=INK2)
        ax.set_yticks(range(len(order)), order, fontsize=8, color=INK)
        ax.set_title(o, fontsize=9, color=INK, loc="left")
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=axes, fraction=0.02, pad=0.02)
    cb.set_label("μ* / max μ* (je Output und Teilbereich)", fontsize=8, color=INK2)
    cb.outline.set_visible(False)
    fig.suptitle(f"Relative Sensitivität je Teilbereich — t = {t:g} d", x=0.01, ha="left",
                 fontsize=12, color=INK)
    fig.savefig(config.RESULTS_DIR / "plots" / f"morris_regions_t{t:g}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main():
    (config.RESULTS_DIR / "plots").mkdir(parents=True, exist_ok=True)
    ind = compute()
    for t in config.SNAPSHOT_TIMES:
        plot_ranking(ind, t); plot_musigma(ind, t); plot_regions(ind, t)
    t = max(config.SNAPSHOT_TIMES)
    sub = ind[(ind.time == t) & (ind.region == "all")]
    top = (sub[~sub.in_noise].groupby("factor").mu_star_rel.max().sort_values(ascending=False))
    print(f"\nFaktoren über dem Rauschniveau in mind. einem Output (t = {t:g}, gesamter Bereich),"
          f"\nWert = größtes relatives μ* über die Outputs:")
    print(top.round(2).to_string())
    rest = sorted(set(NAMES) - set(top.index))
    print(f"\nIn keinem Output über dem Rauschniveau: {', '.join(rest) or '-'}")
    ctrl = sub[sub.role == "control"]
    if len(ctrl) and ctrl.mu_star.max() > 0:
        print("\nAchtung: Kontrollfaktor hat μ* > 0 -> wirkt der b_prob-Schalter wie erwartet?")
    nv = ind.n_valid.min()
    if nv < config.MORRIS_TRAJECTORIES:
        print(f"\nAchtung: wegen fehlgeschlagener Läufe teils nur {nv} gültige Trajektorien.")


if __name__ == "__main__":
    main()
