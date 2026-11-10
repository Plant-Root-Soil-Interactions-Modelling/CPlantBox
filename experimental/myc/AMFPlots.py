import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import matplotlib.pyplot as plt
import warnings
import AMFAnalysis as amf
import AMFDensities as amfd

# All plots of the petri dish results (single simulation heatmaps, profiles over all seeds, encounter fraction)

SUBTYPE_STYLE = {"total": ("k", "all hyphae"), "runner": ("tab:blue", "Runner hyphae"), "bas": ("tab:orange", "BAS")}

def plot_data(all_results):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    n_locations = max(
        r["location"].max() + 1
        for r in all_results.values()
    )
    cmap = plt.cm.viridis
    norm = plt.Normalize(0, n_locations - 1)
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for i, (label, data) in enumerate(all_results.items()):
        ax.scatter(
            data["rho"], data["f_enc"],
            c=data["location"], cmap=cmap, norm=norm,
            marker=markers[i % len(markers)],
            alpha=0.3, s=20, edgecolors="none",
            label=label
        )

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("location")

    ax.set_xlabel(r"realized hyphal length density $\rho_h$")
    ax.set_ylabel(r"$f_{\mathrm{enc}}$")
    ax.set_title(r"Encounter fraction vs. $\rho_h$")
    ax.grid(alpha=0.2)
    ax.legend()

    plt.tight_layout()
    plt.show()

def plot_linear_analysis(all_results):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    for label, data in all_results.items():
        linear = data["linear"]
        if not linear:
            continue
        rho = np.array([r["rho_center"] for r in linear])
        C = np.array([r["C"] for r in linear])
        R2 = np.array([r["R2"] for r in linear])
        ax1.plot(rho, C, "o-", ms=4, label=label)
        ax2.plot(rho, R2, "o-", ms=4, label=label)
    ax1.set_xlabel(r"density $\rho_h$")
    ax1.set_ylabel(r"$C$")
    ax1.set_title(r"Linear coefficient $C$")
    ax2.set_xlabel(r"density $\rho_h$")
    ax2.set_ylabel(r"$R^2$")
    ax2.set_title(r"Linearity")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.2)
        ax.legend()
    plt.tight_layout()
    plt.show()

def plot_ring_heatmap(times, location, mat, cbar_label, filename):
    plt.pcolormesh(times, location, mat, shading='auto', cmap='plasma')
    plt.colorbar(label=cbar_label)
    plt.xlabel("Time [days]")
    plt.ylabel("Distance [cm] from centre")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

def plot_radial_profiles(location, t_common, interp, indices, ylabel="Tip amount", stat="median", filename="TipDensOverAllSim.png"):
    # profiles over the simulations at the time indices given, location should be the ring centers (ringCenters) [cm], plotted in mm
    # stat = "median": median with 25%-75% quantiles (robust, the distribution over seeds is skewed), "mean": mean +- std
    if stat == "median":
        mid = np.median(interp, axis=0)
        low, high = np.percentile(interp, [25, 75], axis=0)
        band = "25%-75% quantile"
    else:
        mid = interp.mean(axis=0)
        std = interp.std(axis=0)
        low, high = mid - std, mid + std
        band = r"mean $\pm$ std"
    fig, ax = plt.subplots(figsize=(6, 5))
    colors = plt.get_cmap("plasma")(np.linspace(0.2, 0.9, len(indices)))
    for color, i in zip(colors, indices):
        ax.plot(location * 10, mid[:, i], color=color, lw=2, marker="o", ms=3, label=f"t = {t_common[i]:.2f} d")
        ax.fill_between(location * 10, low[:, i], high[:, i], color=color, alpha=0.3)
    ax.set_xlabel(r"Distance from centre $r$ (mm)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{stat} over {interp.shape[0]} simulations, band: {band}", fontsize=9)
    ax.legend(title=r"time after $t^*$")
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_spatial_temporal(simulations, location, n_profiles=4, filename="SpatialTemporalProfile.png"):
    # figure as Oyarte Galvez et al. 2025 Fig. 2b. simulations: dicts with "times" (days since t*), "length_densities" [mm mm^-2]
    # and "tip_dens" [mm^-2] (rings x times), location: ring centers (ringCenters) [cm]
    # left, spatial profile: median (band 25%-75% quantile) over the simulations at n_profiles times, coloured by time
    # right, temporal profile: each ring of each simulation shifted by its arrival time t_n (fitArrivalTimes), points are data,
    #        lines the sigmoid (length) and sigmoid derivative (tips) fits, coloured by t_n
    # the first column of length_densities only counts the first interval (see getParasPerRing), so it is left out
    sims = [{"times": sim["times"][1:], "rho": sim["length_densities"][:, 1:] * 1000, "tips": sim["tip_dens"][:, 1:]} for sim in simulations] # rho in um mm^-2
    t_common, rho_i = amf.interpolate_to_common_time(sims, "rho")
    _, tips_i = amf.interpolate_to_common_time(sims, "tips")
    indices = np.linspace(0, len(t_common) - 1, n_profiles).astype(int)
    fits = amfd.collectArrivalFits(simulations)
    t_max = max(t_common[-1], max([f["t_n"] / 24 for f in fits], default=0))
    norm = plt.Normalize(0, t_max)
    cmaps = [plt.get_cmap("viridis"), plt.get_cmap("plasma")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex="col")
    r = location * 10
    for row, (data, cmap) in enumerate(zip([rho_i, tips_i], cmaps)):
        ax = axes[row, 0]
        mid = np.median(data, axis=0)
        low, high = np.percentile(data, [25, 75], axis=0)
        for i in indices:
            c = cmap(norm(t_common[i]))
            ax.plot(r, mid[:, i], color=c, lw=2)
            ax.fill_between(r, low[:, i], high[:, i], color=c, alpha=0.3)
        ax = axes[row, 1]
        key, func = ["length", "tips"][row], [amfd.sigmoidLengthDens, amfd.sigmoidTipDens][row]
        for f in fits:
            c = cmap(norm(f["t_n"] / 24))
            ax.scatter(f["t"] - f["t_n"], f[["rho_data", "tip_data"][row]], color=c, s=4, alpha=0.4, edgecolors="none")
            if f[key] is not None:
                ts = np.linspace(f["t"][0], f["t"][-1], 200)
                ax.plot(ts - f["t_n"], func(ts, *f[key]), color=c, lw=1.2)
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axes[row, :], label="Time (days)")
    if fits: # saturation of the length density
        K1 = np.median([f["length"][0] for f in fits])
        for ax in axes[0]:
            ax.axhline(K1, color="k", ls="--", lw=1)
    axes[0, 0].set_title("Spatial profile")
    axes[0, 1].set_title(f"Temporal profile ({len(fits)} rings of {len(sims)} simulations)")
    axes[0, 0].set_ylabel(r"Filament density ($\mu$m mm$^{-2}$)")
    axes[1, 0].set_ylabel(r"Tip density (mm$^{-2}$)")
    axes[1, 0].set_xlabel(r"Distance from centre $r$ (mm)")
    axes[1, 1].set_xlabel(r"Shifted time $t - t_n$ (h)")
    axes[1, 1].set_xlim(-30, 30)
    for ax in axes.flat:
        ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return fits

def _plot_band(ax, x, data, key, label = None, marker = None):
    # median (line) and 25%-75% quantile (band) over the simulations, data: (n_sim, len(x))
    color, name = SUBTYPE_STYLE[key]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning) # rings without hyphae give NaN in all simulations
        low, mid, high = np.nanpercentile(data, [25, 50, 75], axis=0)
    lw = 3.5 if key == "total" else 2 # thicker, so that it stays visible below a sub type with the same values
    ax.plot(x, mid, color=color, lw=lw, marker=marker, ms=3, label=name if label is None else label)
    ax.fill_between(x, low, high, color=color, alpha=0.25, lw=0)

def plot_hyphal_length_composition(simulations, location, filename = "HyphalLengthComposition.png", target_fraction = 0.3):
    # length of Runner hyphae and BAS in the hyphal compartment (sum over the rings) after t*, simulations as from load_petridish_results
    # (needs "length_runner", "length_bas" [mm mm^-2] and "area" [mm^2] of petridishsimauto.py), location: ring centers [cm]
    # a: length of both sub types, b: BAS fraction of the total hyphal length, c: BAS fraction per ring at the end
    # target_fraction: BAS share of the hyphal length in the experiment (dashed line)
    sims = _with_keys(simulations, ["length_runner", "length_bas", "area"], "plot_hyphal_length_composition")
    if not sims:
        return
    t_end = min(sim["times"][-1] for sim in sims)
    t = np.linspace(max(sim["times"][1] for sim in sims), t_end, 100) # first column only counts the first interval (getParasPerRing)
    L = {}
    for key in ["runner", "bas"]:
        L[key] = np.array([np.interp(t, sim["times"][1:], (sim["length_" + key][:, 1:] * sim["area"][:, None]).sum(axis=0) / 10) for sim in sims]) # [cm]
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = L["bas"] / (L["bas"] + L["runner"])
        ring_frac = np.array([[np.interp(t_end, sim["times"], sim["length_" + k][r]) for k in ["bas", "runner"]] for sim in sims for r in range(len(location))])
        ring_frac = (ring_frac[:, 0] / ring_frac.sum(axis=1)).reshape(len(sims), len(location)) # at the end of the shortest simulation
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for key in ["runner", "bas"]:
        _plot_band(axes[0], t * 24, L[key], key)
    _plot_band(axes[0], t * 24, L["runner"] + L["bas"], "total")
    axes[0].set_ylabel("Hyphal length in the hyphal compartment (cm)")
    axes[0].set_title("a) How much of each hyphal type?")
    _plot_band(axes[1], t * 24, frac, "bas", label="BAS / (BAS + Runner)")
    axes[1].set_ylabel("BAS share of the hyphal length")
    axes[1].set_title("b) BAS share over time")
    _plot_band(axes[2], location * 10, ring_frac, "bas", label=f"BAS / (BAS + Runner) at {t_end * 24:.0f} h", marker="o")
    axes[2].set_xlabel(r"Distance from centre $r$ (mm)")
    axes[2].set_ylabel("BAS share of the hyphal length in the ring")
    axes[2].set_title("c) BAS share along the radius (end)")
    for ax in axes[:2]:
        ax.set_xlabel(r"Time after $t^*$ (h)")
    for ax in axes[1:]:
        ax.axhline(target_fraction, color="k", ls="--", lw=1, label=f"experiment: {target_fraction:.0%}")
        ax.set_ylim(0, 1)
    for ax in axes:
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Median over {len(sims)} simulations, band: 25%-75% quantile. Length = sum over the rings of the filament density times the ring area", fontsize=9)
    fig.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"t": t, "frac": frac, **L}

def plot_anastomosis_frequency(simulations, filename = "AnastomosisFrequency.png", bin_h = 6):
    # anastomoses in the hyphal compartment after t*, all hyphae and split by the sub type of the hypha whose tip fused
    # (needs "ana_new", "ana_new_runner", "ana_new_bas" of petridishsimauto.py: new anastomoses per ring and interval)
    # a: frequency = anastomoses per hour, counted in bins of bin_h hours, b: cumulative number since t*
    keys = ["total", "runner", "bas"]
    names = {"total": "ana_new", "runner": "ana_new_runner", "bas": "ana_new_bas"}
    sims = _with_keys(simulations, list(names.values()), "plot_anastomosis_frequency")
    if not sims:
        return
    edges, centers = _time_bins(sims, bin_h)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    t_cum = np.linspace(0, edges[-1], 100)
    for key in keys:
        rate = np.array([_binned_rate(sim, sim[names[key]].sum(axis=0), edges) for sim in sims])
        _plot_band(axes[0], centers, rate, key, marker="o")
        cum = np.array([np.interp(t_cum, sim["times"] * 24, np.cumsum(sim[names[key]].sum(axis=0))) for sim in sims])
        _plot_band(axes[1], t_cum, cum, key)
    axes[0].set_ylabel("Anastomoses per hour")
    axes[0].set_title(f"a) Anastomosis frequency (bins of {bin_h} h)")
    axes[1].set_ylabel(r"Anastomoses since $t^*$")
    axes[1].set_title("b) Cumulative anastomoses")
    for ax in axes:
        ax.set_xlabel(r"Time after $t^*$ (h)")
        ax.set_ylim(bottom=0)
        ax.legend(fontsize=8, title="tip that fused")
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"Hyphal compartment (all rings), median over {len(sims)} simulations, band: 25%-75% quantile", fontsize=9)
    fig.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_encounter_fraction(simulations, filename = "EncounterFraction.png", min_tips = 5, n_bins = 12):
    # fraction of the growing tips that come close enough to another hypha to fuse (distance < distTH, with ana = 1 every such tip fuses),
    # per hour, against the local filament density, one point per ring, interval and simulation (only with at least min_tips growing tips).
    # f = new anastomoses in the interval / tips growing in the interval / interval length [h], rho = filament density of the ring [um mm^-2]
    # lines: sum of the anastomoses / sum of the tips in bins of rho (equal number of points per bin), for all hyphae and per sub type
    # (tips and anastomoses of that sub type, rho always of all hyphae); dashed: f = C rho, weighted least squares through the origin
    # (needs "ana_new*", "active_tips*" of petridishsimauto.py)
    sims = _with_keys(simulations, ["ana_new", "ana_new_runner", "ana_new_bas", "active_tips", "active_tips_runner", "active_tips_bas"], "plot_encounter_fraction")
    if not sims:
        return
    suffix = {"total": "", "runner": "_runner", "bas": "_bas"}
    data = {key: {"rho": [], "ana": [], "tips": []} for key in suffix}
    for sim in sims:
        dt_h = np.diff(np.concatenate([[0.], sim["times"]]))[1:] * 24 # interval lengths [h], first column left out (see getParasPerRing)
        rho = sim["length_densities"][:, 1:] * 1000
        for key, s in suffix.items():
            ana, tips = sim["ana_new" + s][:, 1:], sim["active_tips" + s][:, 1:]
            data[key]["rho"].append(rho.ravel())
            data[key]["ana"].append((ana / dt_h).ravel())
            data[key]["tips"].append(tips.ravel())
    data = {key: {k: np.concatenate(v) for k, v in d.items()} for key, d in data.items()}
    fig, ax = plt.subplots(figsize=(8, 5.5))
    d = data["total"]
    ok = d["tips"] >= min_tips
    ax.scatter(d["rho"][ok], d["ana"][ok] / d["tips"][ok], s=4, color="0.6", alpha=0.3, edgecolors="none",
               label=f"ring, hour and simulation (all hyphae, $\\geq${min_tips} tips)")
    edges = np.unique(np.quantile(d["rho"][ok], np.linspace(0, 1, n_bins + 1))) if ok.any() else np.array([])
    result = {}
    f_max = 0.
    for key in suffix:
        d = data[key]
        ok = d["tips"] >= min_tips
        if not ok.any() or len(edges) < 2:
            continue
        idx = np.clip(np.searchsorted(edges, d["rho"][ok], side="right") - 1, 0, len(edges) - 2)
        ana = np.bincount(idx, weights=d["ana"][ok], minlength=len(edges) - 1)
        tips = np.bincount(idx, weights=d["tips"][ok], minlength=len(edges) - 1)
        rho_c = np.bincount(idx, weights=d["rho"][ok] * d["tips"][ok], minlength=len(edges) - 1)
        has = tips > 0
        color, name = SUBTYPE_STYLE[key]
        ax.plot(rho_c[has] / tips[has], ana[has] / tips[has], color=color, marker="o", lw=3.5 if key == "total" else 2, label=name)
        f_max = max(f_max, np.max(ana[has] / tips[has]))
        C = np.sum(d["ana"][ok] * d["rho"][ok]) / np.sum(d["tips"][ok] * d["rho"][ok] ** 2) # minimises sum tips (f - C rho)^2
        result[key] = C
        if key == "total":
            r = np.linspace(0, edges[-1], 50)
            ax.plot(r, C * r, color=color, ls="--", lw=1, label=f"$f = C\\rho$, $C$ = {C:.2e} mm$^2$ $\\mu$m$^{{-1}}$ h$^{{-1}}$")
    ax.set_xlabel(r"Local filament density $\rho$ ($\mu$m mm$^{-2}$)")
    ax.set_ylabel("Fraction of growing tips that fuse per hour")
    ax.set_title(f"Tips close enough to another hypha to fuse vs. local density ({len(sims)} simulations)", fontsize=10)
    ax.set_xlim(left=0)
    ax.set_ylim(0, 3 * f_max if f_max > 0 else None) # the single points scatter much more than the binned values
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return result

def _with_keys(simulations, keys, what):
    # simulations that contain all keys (written by the current petridishsimauto.py), prints how many were skipped
    sims = [sim for sim in simulations if all(k in sim for k in keys)]
    if len(sims) < len(simulations):
        print(f"{what}: {len(simulations) - len(sims)} of {len(simulations)} simulations have no {', '.join(keys)} (older files), skipped")
    return sims

def _time_bins(simulations, bin_h):
    # common bins [h since t*] for per-interval counts: from 0 to the end of the shortest simulation
    t_end = min(sim["times"][-1] for sim in simulations) * 24
    edges = np.arange(0, t_end + 1.e-9, bin_h)
    return edges, 0.5 * (edges[1:] + edges[:-1])

def _binned_rate(sim, counts, edges):
    # counts per interval (len(times)) summed into time bins, divided by the bin width -> events per hour
    t_h = sim["times"] * 24 # end of each interval
    summed, _ = np.histogram(t_h - 1.e-9, bins=edges, weights=counts)
    return summed / np.diff(edges)

