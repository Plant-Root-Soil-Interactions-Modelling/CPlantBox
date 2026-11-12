import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import AMFDensities as amfd
import AMFIO as amfio
import AMFSimulationSetup as amfs

# Is the tip density a travelling wave? Velocity, amplitude and width (wavelength) of the tip wave and their deviations,
# from the per-seed results of petridishsimauto.py (run $PETRIDISH_RUN, else the newest run).
# The tip wave is a single pulse (a band of growing tips at the front, Oyarte Galvez et al. 2025 Fig. 2b), not a periodic wave,
# so its "wavelength" is the width of the band: full width at half maximum (FWHM) of the tip density.
# Two views of the same wave, for a travelling wave n(r, t) = f(r - v t) they have to agree:
#   in time:  per ring, the bell fit of n(t) (the fits of SpatialTemporalProfile.png, amfd.collectArrivalFits): peak time t_tip,
#             amplitude K2, duration (FWHM in h); with the velocity, duration * v is the width in space
#   in space: per time, the bell fit of n(r) over the rings (amfd.fitSpatialFronts): front position r_f, amplitude A, width (FWHM in mm)
# Velocity three ways: from the arrival of the filament density t_n(r), from the tip peak t_tip(r) and from the front position r_f(t).
# Every quantity is reduced to one value per seed (median, robust slope), the seeds are the replicates: mean, sd, cv and
# 95% bootstrap confidence interval over the seeds. Constant shape: trend of amplitude and width over time (relative change per day, t-test).
out = amfio.petridishOutputDirs(latest=True)
print("run " + out["run"])
R2_min_length = 0.9  # as in petridishstats.py
R2_min_tips = 0.5
R2_min_space = 0.5   # minimal R2 of the bell fit of n(r) at one time
min_points = 3       # minimal number of rings (or times) per seed for a velocity or trend
v_tip = 0.64 * 10    # tip elongation rate [mm/day], mean of PullerHyphae and RunnerHyphae, TwoHyphaePlusBAS.xml

simulations = amfio.load_petridish_results(out["densities"] + "seed_*.npz")
if not simulations:
    sys.exit("no simulations with the current geometry found")
location = amfs.ringCenters(simulations[0]["location"]) # [cm]
r_ring = location * 10                                  # [mm]
sim_by_seed = {int(sim["seed"]): sim for sim in simulations}

## fits in time (per ring) and in space (per time)
rows = []
for f in amfd.collectArrivalFits(simulations):
    K2, lam_tip, t_tip = f["tips"] if f["tips"] is not None else (np.nan, np.nan, np.nan)
    rows.append({"seed": f["seed"], "ring": f["ring"], "r": r_ring[f["ring"]], "t_n": f["length"][2], "R2_length": f["R2_length"],
                 "K2": K2, "t_tip": t_tip, "duration": amfd.FWHM_BELL / lam_tip, "R2_tips": f["R2_tips"], "t_first": f["t"][0], "t_last": f["t"][-1]})
ring_fits = pd.DataFrame(rows, columns=["seed", "ring", "r", "t_n", "R2_length", "K2", "t_tip", "duration", "R2_tips", "t_first", "t_last"])
good_len = ring_fits[ring_fits["R2_length"] >= R2_min_length]
# tip peak and both half maximum points within the simulated time, else amplitude and duration are extrapolated
inside = (good_len["t_tip"] - good_len["duration"] / 2 >= good_len["t_first"]) & (good_len["t_tip"] + good_len["duration"] / 2 <= good_len["t_last"])
good_tip = good_len[(good_len["R2_tips"] >= R2_min_tips) & inside]
space = pd.DataFrame(amfd.fitSpatialFronts(simulations, location, R2_min_space), columns=["seed", "j", "t", "r_f", "A", "k", "width", "R2"])
space.to_csv(amfio.plotFile(out, "TipWaveFitsInSpace.csv"), index=False)
print(f"{len(good_tip)} rings with a good fit in time (of {len(ring_fits)}), {len(space)} times with a good fit in space, {len(simulations)} seeds")

## one value per seed
per_seed = []
for seed in sim_by_seed:
    gl, gt, sp = good_len[good_len["seed"] == seed], good_tip[good_tip["seed"] == seed], space[space["seed"] == seed]
    row = {"seed": seed}
    # velocities [mm/h]: robust slope of the front position against time
    row["v_length"] = amfd.theilSlope(gl["t_n"], gl["r"])[0] if len(gl) >= min_points else np.nan
    row["v_tips"] = amfd.theilSlope(gt["t_tip"], gt["r"])[0] if len(gt) >= min_points else np.nan
    row["v_space"] = amfd.theilSlope(sp["t"], sp["r_f"])[0] if len(sp) >= min_points else np.nan
    # amplitude [mm^-2] and width [mm]
    row["A_space"] = sp["A"].median() if len(sp) else np.nan
    row["K2_time"] = gt["K2"].median() if len(gt) else np.nan
    row["width_space"] = sp["width"].median() if len(sp) else np.nan
    row["duration"] = gt["duration"].median() if len(gt) else np.nan
    v = row["v_tips"] if np.isfinite(row["v_tips"]) else row["v_space"]
    row["width_time"] = row["duration"] * v # duration of the tip peak in a ring times the velocity
    # trends over time (constant shape: 0), relative change per day
    if len(sp) >= min_points:
        row["A_trend"] = amfd.theilSlope(sp["t"], sp["A"])[0] * 24 / sp["A"].median()
        row["width_trend"] = amfd.theilSlope(sp["t"], sp["width"])[0] * 24 / sp["width"].median()
    else:
        row["A_trend"] = row["width_trend"] = np.nan
    # damping or amplification: exponential rate of A(t), exponent of A against the front position (A ~ r_f^alpha, alpha = -1: the same
    # number of tips spread over a longer front, geometric dilution), and rate of the number of tips in the band N = A * 4/k * pi r_f
    # (integral of the bell curve across the band times the length of the half circle at r_f)
    if len(sp) >= min_points:
        row["A_rate"] = amfd.theilSlope(sp["t"], np.log(sp["A"]))[0] * 24
        row["A_exponent"] = amfd.theilSlope(np.log(sp["r_f"]), np.log(sp["A"]))[0]
        row["N_rate"] = amfd.theilSlope(sp["t"], np.log(sp["A"] * 4 / sp["k"] * np.pi * sp["r_f"]))[0] * 24
    else:
        row["A_rate"] = row["A_exponent"] = row["N_rate"] = np.nan
    per_seed.append(row)
per_seed = pd.DataFrame(per_seed)
per_seed.to_csv(amfio.plotFile(out, "TipWavePerSeed.csv"), index=False)

print("\n1. Velocity of the front (mm/day)")
res = {}
for key, name in [("v_length", "from the arrival of the filament density t_n(r)"), ("v_tips", "from the tip peak t_tip(r)"),
                  ("v_space", "from the front position r_f(t) (fit in space)")]:
    res[key] = amfd.seedSummary(name, per_seed[key] * 24, "mm/day")
print(f"  tips grow with {v_tip:.1f} mm/day")
print("\n2. Amplitude of the tip wave (peak tip density, mm^-2)")
res["A_space"] = amfd.seedSummary("A, fit in space (median over the times)", per_seed["A_space"], "mm^-2")
res["K2_time"] = amfd.seedSummary("K2, fit in time (median over the rings)", per_seed["K2_time"], "mm^-2")
res["A_trend"] = amfd.seedSummary("relative change of A per day (constant: 0)", per_seed["A_trend"], "", test_zero=True)
print("\n3. Width of the tip wave = wavelength of the pulse (FWHM, mm)")
res["width_space"] = amfd.seedSummary("width in space (median over the times)", per_seed["width_space"], "mm")
res["width_time"] = amfd.seedSummary("duration in a ring x velocity", per_seed["width_time"], "mm")
res["duration"] = amfd.seedSummary("duration of the tip peak in a ring", per_seed["duration"], "h")
res["width_trend"] = amfd.seedSummary("relative change of the width per day (constant: 0)", per_seed["width_trend"], "", test_zero=True)
print("\n4. Damping (< 0) or amplification (> 0) of the amplitude")
res["A_rate"] = amfd.seedSummary("rate gamma of A(t) ~ exp(gamma t)", per_seed["A_rate"], "1/day", test_zero=True)
res["A_exponent"] = amfd.seedSummary("exponent alpha of A ~ r_f^alpha (constant: 0, geometric dilution: -1)", per_seed["A_exponent"], "", test_zero=True)
res["A_exponent_dilution"] = amfd.seedSummary("alpha + 1 (geometric dilution: 0)", per_seed["A_exponent"] + 1, "", test_zero=True, verbose=False)
if res["A_exponent_dilution"] is not None:
    print(f"  alpha = -1 (geometric dilution): p = {res['A_exponent_dilution']['p']:.3g}")
res["N_rate"] = amfd.seedSummary("rate of the number of tips in the band N(t) ~ exp(rate t)", per_seed["N_rate"], "1/day", test_zero=True)
amfio.updatePlotInfo(out, "petridishwave", {"R2_min_length": R2_min_length, "R2_min_tips": R2_min_tips, "R2_min_space": R2_min_space,
                                             "min_points": min_points, "summary": res})

## figure
def ci(r, fmt = "{:.3g}", scale = 1.):
    return (fmt + " [" + fmt + ", " + fmt + "], cv {:.2f}").format(r["mean"] * scale, r["lo"] * scale, r["hi"] * scale, r["cv"]) if r is not None else "less than 2 seeds"

def annotate(ax, text):
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.8"))

colors = {seed: plt.get_cmap("tab20")(i % 20) for i, seed in enumerate(sorted(sim_by_seed))}
fig, axes = plt.subplots(2, 3, figsize=(17, 10.5))
# a) front position against time
ax = axes[0, 0]
for seed, d in space.groupby("seed"):
    ax.plot(d["t"], d["r_f"], "o", ms=2.5, color=colors[seed], alpha=0.6)
    s, i = amfd.theilSlope(d["t"], d["r_f"])
    if np.isfinite(s):
        ax.plot(d["t"], i + s * d["t"], color=colors[seed], lw=0.8)
for seed, d in good_tip.groupby("seed"):
    ax.plot(d["t_tip"], d["r"], "x", ms=4, color=colors[seed], alpha=0.8)
ax.set_xlabel(r"Time after $t^*$ (h)")
ax.set_ylabel("Position of the front (mm from centre)")
ax.set_title("a) How fast does the tip wave move?", fontsize=10)
annotate(ax, "dots + line: peak of n(r) at each time (fit in space)\ncrosses: peak of n(t) in each ring (fit in time)\n"
             f"velocity (mm/day, mean [95% CI]):\n  front position: {ci(res['v_space'])}\n  tip peak per ring: {ci(res['v_tips'])}\n"
             f"  filament arrival: {ci(res['v_length'])}\n  tips grow with {v_tip:.1f} mm/day")
# b) collapse: tip density against the distance to the front
ax = axes[0, 1]
norm = plt.Normalize(space["t"].min() if len(space) else 0, space["t"].max() if len(space) else 1)
for _, f in space.iterrows():
    n = sim_by_seed[f["seed"]]["tip_dens"][:, 1:][:, int(f["j"])]
    ax.plot(r_ring - f["r_f"], n, ".", ms=2, color=plt.get_cmap("viridis")(norm(f["t"])), alpha=0.3)
if len(space):
    x = np.linspace(-25, 25, 200)
    ax.plot(x, amfd.sigmoidTipDens(x, space["A"].median(), space["k"].median(), 0.), color="k", lw=2, label="median fit")
    ax.legend(frameon=False, loc="upper right")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax, label=r"Time after $t^*$ (h)")
ax.set_xlabel("Distance to the front $r - r_f(t)$ (mm)")
ax.set_ylabel(r"Tip density $n$ (mm$^{-2}$)")
ax.set_title("b) Does the wave keep its shape? (all times on top of each other)", fontsize=10)
# c) amplitude over time
ax = axes[0, 2]
for seed, d in space.groupby("seed"):
    ax.plot(d["t"], d["A"], "-", lw=0.8, color=colors[seed], alpha=0.7)
for seed, d in good_tip.groupby("seed"):
    ax.plot(d["t_tip"], d["K2"], "x", ms=4, color=colors[seed], alpha=0.8)
ax.set_xlabel(r"Time after $t^*$ (h)")
ax.set_ylabel(r"Amplitude = peak tip density (mm$^{-2}$)")
ax.set_title("c) Amplitude of the tip wave", fontsize=10)
annotate(ax, "lines: A(t), fit in space; crosses: K2 per ring at its peak time\n"
             f"A: {ci(res['A_space'])} mm$^{{-2}}$\nK2: {ci(res['K2_time'])} mm$^{{-2}}$\n"
             + (f"change per day: {res['A_trend']['mean'] * 100:+.0f}% (constant: p = {res['A_trend']['p']:.2g})" if res["A_trend"] else ""))
# d) width over time
ax = axes[1, 0]
for seed, d in space.groupby("seed"):
    ax.plot(d["t"], d["width"], "-", lw=0.8, color=colors[seed], alpha=0.7)
for seed, d in good_tip.groupby("seed"):
    vs = per_seed.loc[per_seed["seed"] == seed, "v_tips"].values[0]
    if np.isfinite(vs):
        ax.plot(d["t_tip"], d["duration"] * vs, "x", ms=4, color=colors[seed], alpha=0.8)
ax.set_xlabel(r"Time after $t^*$ (h)")
ax.set_ylabel("Width of the tip wave (FWHM, mm)")
ax.set_title("d) Width (wavelength) of the tip wave", fontsize=10)
annotate(ax, "lines: width of n(r) (fit in space)\ncrosses: duration in a ring x velocity\n"
             f"width in space: {ci(res['width_space'])} mm\nduration x velocity: {ci(res['width_time'])} mm\n"
             f"duration in a ring: {ci(res['duration'])} h\n"
             + (f"change per day: {res['width_trend']['mean'] * 100:+.0f}% (constant: p = {res['width_trend']['p']:.2g})" if res["width_trend"] else ""))
# e) width: space against time view
ax = axes[1, 1]
ax.plot(per_seed["width_time"], per_seed["width_space"], "o", color="tab:blue")
m = np.nanmax([per_seed["width_time"].max(), per_seed["width_space"].max(), 1.])
ax.plot([0, m], [0, m], "k--", lw=1)
ax.set_xlabel("Duration in a ring x velocity (mm)")
ax.set_ylabel("Width in space (mm)")
ax.set_title("e) Do both views give the same width?", fontsize=10)
annotate(ax, "one point = one seed (medians)\ndashed: 1:1, as for a travelling wave")
# f) velocities per seed
ax = axes[1, 2]
keys = [("v_space", "front position"), ("v_tips", "tip peak per ring"), ("v_length", "filament arrival")]
for i, (key, name) in enumerate(keys):
    y = per_seed[key] * 24
    ax.plot(np.full(len(y), i) + np.linspace(-0.15, 0.15, len(y)), y, "o", alpha=0.7)
    if res[key] is not None:
        ax.errorbar(i + 0.3, res[key]["mean"], yerr=[[res[key]["mean"] - res[key]["lo"]], [res[key]["hi"] - res[key]["mean"]]], fmt="ks", capsize=4)
ax.axhline(v_tip, color="0.5", ls=":", lw=1, label="tip elongation rate")
ax.set_xticks(range(len(keys)))
ax.set_xticklabels([n for _, n in keys])
ax.set_ylabel("Velocity of the front (mm/day)")
ax.set_title("f) Do the three velocities agree?", fontsize=10)
ax.legend(frameon=False, loc="lower right")
annotate(ax, "points: seeds, black: mean and 95% CI")
for ax in axes.flat:
    ax.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.005,
         "Tip wave = band of growing tips at the front. For a travelling wave $n(r, t) = f(r - vt)$ the front moves with constant velocity v and keeps its amplitude and width.\n"
         "Fit in space: at every time $n(r) = 4A\\,e^{-k|r - r_f|}/(1 + e^{-k|r - r_f|})^2$ over the rings (front position $r_f$, amplitude A, width = 3.53/k), "
         f"kept if peak and half maximum points lie within the rings and $R^2 \\geq$ {R2_min_space}.\n"
         "Fit in time: per ring the same bell curve over time (peak time $t_{tip}$, amplitude $K_2$, duration = FWHM in h), as in SpatialTemporalProfile.png, "
         "kept if peak and half maximum points lie within the simulated time.\n"
         f"Velocities and trends: robust (Theil-Sen) lines per seed with at least {min_points} points.\n"
         "Values: mean over the seeds [95% bootstrap confidence interval], cv = sd / mean over the seeds (variation between seeds). "
         f"{len(simulations)} seeds, run {out['run']}.",
         fontsize=8.5, va="bottom", ha="left")
plt.tight_layout(rect=(0, 0.08, 1, 1))
plt.savefig(amfio.plotFile(out, "TipWave.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print("\nwritten " + amfio.plotFile(out, "TipWave.png"))

## figure: damping or amplification of the amplitude
def pval(r):
    return f"p = {r['p']:.2g}" if r is not None else ""

fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
# a) A(t) relative to the first fitted time, log scale: straight line = exponential damping or amplification
ax = axes[0]
for seed, d in space.groupby("seed"):
    d = d.sort_values("t")
    ax.semilogy(d["t"], d["A"] / d["A"].iloc[0], "o", ms=2.5, color=colors[seed], alpha=0.6)
    s, i = amfd.theilSlope(d["t"], np.log(d["A"] / d["A"].iloc[0]))
    if np.isfinite(s):
        ax.semilogy(d["t"], np.exp(i + s * d["t"]), color=colors[seed], lw=1)
ax.axhline(1, color="k", ls=":", lw=1)
ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2g"))
ax.yaxis.set_minor_formatter(plt.FormatStrFormatter("%.2g"))
ax.set_xlabel(r"Time after $t^*$ (h)")
ax.set_ylabel("Amplitude relative to the first time $A(t) / A(t_0)$")
ax.set_title("a) Is the tip wave damped or amplified over time?", fontsize=10)
annotate(ax, "one colour = one seed, line: $A \\sim e^{\\gamma t}$\n"
             f"$\\gamma$ = {ci(res['A_rate'])} 1/day\n"
             f"(< 0 damped, > 0 amplified; $\\gamma$ = 0: {pval(res['A_rate'])})")
# b) amplitude against the front position, log-log: slope alpha, reference slope -1 for geometric dilution
ax = axes[1]
for seed, d in space.groupby("seed"):
    ax.loglog(d["r_f"], d["A"], "o", ms=2.5, color=colors[seed], alpha=0.6)
    s, i = amfd.theilSlope(np.log(d["r_f"]), np.log(d["A"]))
    if np.isfinite(s):
        rr = np.array([d["r_f"].min(), d["r_f"].max()])
        ax.loglog(rr, np.exp(i) * rr ** s, color=colors[seed], lw=1)
for seed, d in good_tip.groupby("seed"):
    ax.loglog(d["r"], d["K2"], "x", ms=4, color=colors[seed], alpha=0.8)
if len(space):
    rr = np.array([space["r_f"].min(), space["r_f"].max()])
    ref = space["A"].median() * space["r_f"].median()
    ax.loglog(rr, ref / rr, "k--", lw=1, label="$A \\sim 1/r$: same number of tips\non a longer front (dilution)")
    ax.legend(frameon=False, loc="lower left", fontsize=8)
ax.set_xlabel("Position of the front $r_f$ (mm)")
ax.set_ylabel(r"Amplitude $A$ (mm$^{-2}$)")
ax.set_title("b) Amplitude against the distance travelled", fontsize=10)
annotate(ax, "dots + line: $A \\sim r_f^{\\alpha}$ per seed, crosses: $K_2$ per ring\n"
             f"$\\alpha$ = {ci(res['A_exponent'])}\n"
             f"constant ($\\alpha$ = 0): {pval(res['A_exponent'])},  dilution ($\\alpha$ = -1): {pval(res['A_exponent_dilution'])}")
# c) number of tips in the band
ax = axes[2]
for seed, d in space.groupby("seed"):
    d = d.sort_values("t")
    N = d["A"] * 4 / d["k"] * np.pi * d["r_f"]
    ax.semilogy(d["t"], N, "o", ms=2.5, color=colors[seed], alpha=0.6)
    s, i = amfd.theilSlope(d["t"], np.log(N))
    if np.isfinite(s):
        ax.semilogy(d["t"], np.exp(i + s * d["t"]), color=colors[seed], lw=1)
ax.set_xlabel(r"Time after $t^*$ (h)")
ax.set_ylabel(r"Tips in the band $N = A \cdot 4/k \cdot \pi r_f$")
ax.set_title("c) Does the number of tips at the front grow?", fontsize=10)
annotate(ax, "line: $N \\sim e^{\\mu t}$ per seed\n"
             f"$\\mu$ = {ci(res['N_rate'])} 1/day ($\\mu$ = 0: {pval(res['N_rate'])})\n"
             "(branching adds tips, anastomosis and death remove them)")
for ax in axes:
    ax.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.005,
         "Amplitude A, width (3.53/k) and front position $r_f$ from the fit of the tip density over the rings at each time (as in TipWave.png). "
         "The front spreads over a half circle of length $\\pi r_f$, so the same number of tips gives\n"
         "$A \\sim 1/r_f$ (geometric dilution, b). Without dilution and growth A stays constant ($\\gamma$ = 0, $\\alpha$ = 0). "
         "N counts the tips in the band; its rate $\\mu$ separates real growth of the front from the dilution. "
         f"Rates and exponents: robust (Theil-Sen) fits per seed, mean over the seeds [95% CI], p: t-test over the seeds. {len(simulations)} seeds, run {out['run']}.",
         fontsize=8.5, va="bottom", ha="left", wrap=True)
plt.tight_layout(rect=(0, 0.1, 1, 1))
plt.savefig(amfio.plotFile(out, "TipWaveAmplitude.png"), dpi=200, bbox_inches="tight")
plt.close(fig)
print("written " + amfio.plotFile(out, "TipWaveAmplitude.png"))
