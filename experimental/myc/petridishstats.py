import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import AMFDensities as amfd
import AMFIO as amfio
import AMFSimulationSetup as amfs

# Statistics of the travelling wave in the per-seed results of petridishsimauto.py (see Oyarte Galvez et al. 2025 Fig. 2b),
# from the arrival time fits of each ring (amfd.collectArrivalFits, the same fits as in SpatialTemporalProfile.png of petridishcombine.py):
# 1. front velocity: t_n against r per seed, v = 1/slope, constant velocity tested by the quadratic term
# 2. constant wave shape: saturation K1, tip peak K2 and width of the tip peak against r, slope per seed
# 3. consistency: the length behind the front is grown by the tips, K1 ~ v_tip * integral of the tip density = v_tip * 4 K2 / lamda_tips
# The rings of one seed are not independent, so every quantity is first reduced to one value per seed (slope, median),
# the seeds are the independent replicates (t-test, bootstrap confidence interval over the seeds).
out = amfio.petridishOutputDirs(latest=True) # run $PETRIDISH_RUN, else the newest run, below $PETRIDISH_OUT (cluster) or petridish_results
print("run " + out["run"])
results_dir = out["densities"]
R2_min_length = 0.9   # minimal coefficient of determination of the sigmoid fit of the length density
R2_min_tips = 0.5     # minimal coefficient of determination of the sigmoid derivative fit of the tip density
saturation_min = 0.8  # for K1: the length density at the end has to reach this fraction of the fitted saturation (else K1 is extrapolated)
min_rings = 3         # minimal number of rings per seed for a slope
v_tip = 0.64 * 1.e4 / 24  # tip elongation rate [um/h], mean of PullerHyphae (0.68 cm/day) and RunnerHyphae (0.60 cm/day), TwoHyphaePlusBAS.xml
n_boot = 10000
rng = np.random.default_rng(1)
paper = {"K1": 1400.}  # saturation of the filament density [um mm^-2], read from Fig. 2b (dashed line); tip peak there about 0.1-0.2 mm^-2

simulations = amfio.load_petridish_results(results_dir + "seed_*.npz")
if not simulations:
    sys.exit("no simulations with the current geometry found")
r_ring = amfs.ringCenters(simulations[0]["location"]) * 10  # ring centers [mm]

## table of all fits (one row per ring and seed)
fwhm = 2 * np.log(3 + 2 * np.sqrt(2))  # full width at half maximum of sigmoidTipDens is fwhm / lamda_tips
rows = []
for f in amfd.collectArrivalFits(simulations):
    K1, lam_len, t_n = f["length"]
    K2, lam_tip, t_tip = f["tips"] if f["tips"] is not None else (np.nan, np.nan, np.nan)
    rows.append({"seed": f["seed"], "ring": f["ring"], "r": r_ring[f["ring"]], "t_n": t_n, "K1": K1, "lam_len": lam_len, "R2_length": f["R2_length"],
                 "saturation": f["rho_data"][-1] / K1, "K2": K2, "lam_tip": lam_tip, "t_tip": t_tip, "width": fwhm / lam_tip, "R2_tips": f["R2_tips"]})
df = pd.DataFrame(rows)
df.to_csv(amfio.plotFile(out, "WaveFitsPerRing.csv"), index=False)
amfio.updatePlotInfo(out, "petridishstats", {"R2_min_length": R2_min_length, "R2_min_tips": R2_min_tips, "saturation_min": saturation_min,
                                              "min_rings": min_rings, "v_tip": v_tip, "n_boot": n_boot, "seeds": sorted(df["seed"].unique().tolist()) if "seed" in df else []})
good = df[df["R2_length"] >= R2_min_length]
good_tips = good[good["R2_tips"] >= R2_min_tips]
good_K1 = good[good["saturation"] >= saturation_min]
print(f"{len(df)} ring fits of {df['seed'].nunique()} seeds, {len(good)} with R2_length >= {R2_min_length}, {len(good_tips)} of them with R2_tips >= {R2_min_tips}, "
      f"{len(good_K1)} saturated (>= {saturation_min} K1); all fits in WaveFitsPerRing.csv")

def boot_ci(x, stat=np.mean):
    # 95% bootstrap confidence interval of stat over the seeds
    x = np.asarray(x)
    b = stat(x[rng.integers(0, len(x), (n_boot, len(x)))], axis=1)
    return np.percentile(b, [2.5, 97.5])

def summary(name, x, unit = "", test_zero = False):
    # one value per seed: mean, median, sd and 95% CI over the seeds, optionally t-test against 0
    # returns mean, 95% CI, p-value (nan without test) and the number of seeds, for the annotations in the figure (None for < 2 seeds)
    x = np.asarray(x, dtype=float) # an empty column (no seed with enough rings) has dtype object
    x = x[np.isfinite(x)]
    if len(x) < 2:
        print(f"  {name}: less than 2 seeds")
        return None
    lo, hi = boot_ci(x)
    s = f"  {name}: mean {x.mean():.4g} {unit}, median {np.median(x):.4g}, sd {x.std(ddof=1):.3g}, 95% CI of mean [{lo:.4g}, {hi:.4g}], n = {len(x)} seeds"
    p = np.nan
    if test_zero:
        p = stats.ttest_1samp(x, 0).pvalue
        s += f", t-test mean = 0: p = {p:.3g}"
    print(s)
    return {"mean": x.mean(), "lo": lo, "hi": hi, "p": p, "n": len(x)}

def annotate(ax, text, lower_right = False):
    # result of the panel in its upper left (or lower right) corner
    x, y, va, ha = (0.97, 0.03, "bottom", "right") if lower_right else (0.03, 0.97, "top", "left")
    ax.text(x, y, text, transform=ax.transAxes, va=va, ha=ha, fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="0.8"))

def per_seed_slopes(data, y):
    # slope and intercept of y against r for each seed with at least min_rings rings
    res = []
    for seed, d in data.groupby("seed"):
        if len(d) >= min_rings and d["r"].nunique() >= 2:
            lr = stats.linregress(d["r"], d[y])
            res.append({"seed": seed, "slope": lr.slope, "intercept": lr.intercept, "mean": d[y].mean(), "n": len(d)})
    return pd.DataFrame(res, columns=["seed", "slope", "intercept", "mean", "n"])

## 1. front velocity
# the slope of t_n(r) [h/mm] is averaged over the seeds, not v = 1/slope: a seed with slope near 0 gives an arbitrarily large v.
# Theil-Sen slope (median of the pairwise slopes), robust against single rings with a bad arrival time fit
print("\n1. Front velocity (t_n against r per seed, Theil-Sen slope)")
vel = []
for seed, d in good.groupby("seed"):
    if len(d) < min_rings:
        continue
    slope, intercept = stats.theilslopes(d["t_n"], d["r"])[:2]
    R2 = stats.pearsonr(d["r"], d["t_n"])[0]**2
    c2 = np.polyfit(d["r"], d["t_n"], 2)[0] if len(d) >= 4 else np.nan # quadratic term [h mm^-2]: > 0 front slows down, < 0 it accelerates
    vel.append({"seed": seed, "slope": slope, "intercept": intercept, "R2": R2, "c2": c2, "n": len(d)})
vel = pd.DataFrame(vel, columns=["seed", "slope", "intercept", "R2", "c2", "n"])
res_slope = summary("slope of t_n(r)", vel["slope"], "h/mm")
v = lambda s: 24 / s if s > 0 else np.inf
if res_slope is not None:
    print(f"  velocity v = 1 / mean slope: {v(res_slope['mean']):.3g} mm/day, 95% CI [{v(res_slope['hi']):.3g}, {v(res_slope['lo']):.3g}] (from the CI of the slope)")
res_R2 = summary("R2 of t_n(r)", vel["R2"])
res_c2 = summary("quadratic term of t_n(r) (constant velocity: 0)", vel["c2"], "h/mm^2", test_zero=True)

## 2. constant wave shape
print("\n2. Wave shape against r (slope per seed, constant shape: slope 0; relative change = slope * 10 mm / mean)")
shape = {"K1": (good_K1, "um mm^-2 per mm"), "K2": (good_tips, "mm^-2 per mm"), "width": (good_tips, "h per mm")}
slopes, res_shape = {}, {}
for y, (data, unit) in shape.items():
    s = per_seed_slopes(data, y)
    slopes[y] = s
    res_shape[y] = {"median": summary(f"{y}: seed median", data.groupby("seed")[y].median(), unit.split(" per")[0]),
                    "slope": summary(f"{y}: slope against r", s["slope"], unit, test_zero=True),
                    "relative": summary(f"{y}: relative change over 10 mm", s["slope"] * 10 / s["mean"])}
print(f"  paper: K1 about {paper['K1']:.0f} um mm^-2, tip peak about 0.1-0.2 mm^-2 (Fig. 2b)")

## 3. consistency of length and tip density
print(f"\n3. Consistency: K1 against v_tip * 4 K2 / lamda_tips (v_tip = {v_tip:.0f} um/h)")
cons = good_tips[good_tips["saturation"] >= saturation_min].copy()
cons["K1_pred"] = v_tip * 4 * cons["K2"] / cons["lam_tip"]
cons["ratio"] = cons["K1"] / cons["K1_pred"]
cons["lam_ratio"] = cons["lam_len"] / cons["lam_tip"]
res_ratio = summary("K1 / predicted, seed median (1: all length grown by the front tips, > 1: length grown behind the front)", cons.groupby("seed")["ratio"].median())
res_lam = summary("lamda_length / lamda_tips, seed median (1: tips are the derivative of the length)", cons.groupby("seed")["lam_ratio"].median())

## figure
# every panel says in words what is plotted, the question it answers (title) and its result over the seeds (box),
# the fit functions behind K1, K2, t_n and the width are explained below the panels
def ci_text(res, fmt = "{:.3g}"):
    return (fmt + " [" + fmt + ", " + fmt + "]").format(res["mean"], res["lo"], res["hi"]) if res is not None else "less than 2 seeds"

fig, axes = plt.subplots(2, 3, figsize=(16, 10.5))
ax = axes[0, 0]
for seed, d in good.groupby("seed"):
    ax.plot(d["r"], d["t_n"], "o", ms=3, alpha=0.5)
    row = vel[vel["seed"] == seed]
    if len(row):
        rr = np.array([d["r"].min(), d["r"].max()])
        ax.plot(rr, row["intercept"].values[0] + row["slope"].values[0] * rr, lw=0.8, color=ax.lines[-1].get_color())
ax.set_xlabel("Distance from centre $r$ (mm)")
ax.set_ylabel("Arrival time of the network at the ring $t_n$ (h after t*)")
ax.set_title("1a. When does the network reach each ring?", fontsize=10)
annotate(ax, "one colour = one seed, line = robust fit (Theil-Sen)\n"
             + (f"$R^2$ of the line: {ci_text(res_R2, '{:.2f}')}" if res_R2 is not None else ""))
ax = axes[0, 1]
if len(vel):
    ax.hist(vel["slope"], bins=max(5, len(vel) // 3), color="#4c72b0", edgecolor="white")
    ax.axvline(vel["slope"].mean(), color="k", ls="--", lw=1)
ax.set_xlabel("Slope of the lines in 1a (h/mm) = time the front needs per mm")
ax.set_ylabel("Number of seeds")
ax.yaxis.get_major_locator().set_params(integer=True)
ax.set_title("1b. How fast does the front move?", fontsize=10)
if res_slope is not None:
    annotate(ax, f"dashed: mean slope {ci_text(res_slope)} h/mm\n"
                 f"front velocity = 1 / slope: {v(res_slope['mean']):.3g} mm/day\n"
                 f"tips grow with {v_tip * 24 / 1000:.1f} mm/day")
ax = axes[0, 2]
ax.plot(cons["K1_pred"], cons["K1"], "o", ms=3, alpha=0.5)
m = np.nanmax([cons["K1_pred"].max(), cons["K1"].max(), 1.])
ax.plot([0, m], [0, m], "k--", lw=1)
ax.set_xlabel("Filament density the passing tip wave lays down\n"
              r"$v_{tip} \cdot \int n\,dt = v_{tip} \cdot 4 K_2 / \lambda_{tips}$ ($\mu$m mm$^{-2}$)")
ax.set_ylabel("Final filament density in the ring $K_1$ ($\\mu$m mm$^{-2}$)")
ax.set_title("3. Is the filament length laid down by the tips of the front?", fontsize=10)
annotate(ax, "one point = one ring\ndashed: 1:1 (all length from the front tips)\n"
             f"ratio $K_1$ / expected: {ci_text(res_ratio, '{:.2f}')}\n"
             r"$\lambda_{length} / \lambda_{tips}$: " + ci_text(res_lam, '{:.2f}'))
labels = {"K1": "Final filament density in the ring $K_1$ ($\\mu$m mm$^{-2}$)",
          "K2": "Height of the tip wave $K_2$\n(peak tip density, mm$^{-2}$)",
          "width": "Duration of the tip wave in the ring\n(full width at half maximum, h)"}
titles = {"K1": "2a. Does the final filament density change with r?",
          "K2": "2b. Does the height of the tip wave change with r?",
          "width": "2c. Does the duration of the tip wave change with r?"}
for ax, (y, (data, unit)) in zip(axes[1], shape.items()):
    ax.plot(data["r"], data[y], "o", ms=3, alpha=0.3, color="0.5")
    for _, s in slopes[y].iterrows():
        d = data[data["seed"] == s["seed"]]
        rr = np.array([d["r"].min(), d["r"].max()])
        ax.plot(rr, s["intercept"] + s["slope"] * rr, lw=0.6, alpha=0.6)
    if y == "K1":
        ax.axhline(paper["K1"], color="k", ls="--", lw=1, label="paper (Fig. 2b)")
        ax.legend(frameon=False, loc="lower left")
    ax.set_xlabel("Distance from centre $r$ (mm)")
    ax.set_ylabel(labels[y])
    ax.set_title(titles[y], fontsize=10)
    res = res_shape[y]
    text = "grey = rings, line = linear fit per seed\n"
    if res["relative"] is not None:
        text += f"change over 10 mm: {res['relative']['mean'] * 100:+.0f}% [{res['relative']['lo'] * 100:+.0f}%, {res['relative']['hi'] * 100:+.0f}%]\n"
        text += f"slope = 0 (same shape everywhere): p = {res['slope']['p']:.2g}, {res['slope']['n']} seeds"
    else:
        text += "less than 2 seeds with enough rings"
    annotate(ax, text)
for ax in axes.flat:
    ax.spines[["top", "right"]].set_visible(False)
fig.text(0.01, 0.005,
         "One point = one ring of one simulation. For each ring the filament density $\\rho(t)$ and the density of growing tips $n(t)$ are fitted:\n"
         "   $\\rho(t) = K_1 \\,/\\, (1 + e^{\\lambda_{length}(t_n - t)})$ (sigmoid):   "
         "$K_1$ = filament density long after the front has passed,   $t_n$ = arrival time (half of $K_1$ reached),   $\\lambda_{length}$ = rate of the rise.\n"
         "   $n(t) = 4 K_2\\, e^{-\\lambda_{tips}|t - t_n|} / (1 + e^{-\\lambda_{tips}|t - t_n|})^2$ (derivative of a sigmoid, bell curve):   "
         "$K_2$ = peak tip density when the front passes,   width = $3.53 / \\lambda_{tips}$ (full width at half maximum).\n"
         f"Rings with a poor fit ($R^2$ < {R2_min_length} length, < {R2_min_tips} tips) are left out, $K_1$ only where the density reached {saturation_min:.0%} of it. "
         "Values in brackets: 95% bootstrap confidence interval of the mean over the seeds; p: t-test over the per-seed slopes.",
         fontsize=9, va="bottom", ha="left")
plt.tight_layout(rect=(0, 0.1, 1, 1))
plt.savefig(amfio.plotFile(out, "WaveStatistics.png"), dpi=200, bbox_inches="tight")
print("\nwritten " + amfio.plotFile(out, "WaveStatistics.png"))
