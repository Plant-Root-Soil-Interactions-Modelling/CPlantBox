import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import matplotlib.pyplot as plt
import AMFAnalysis as amf

# Run statistics of the per-seed results written by petridishsimauto.py: how long (simulated) it takes until t* (10 cm hyphae
# behind the barrier), how long each phase takes on the wall clock, and which seeds were stopped by max_wait_days.
# Used to choose max_wait_days and the SLURM time limit.
results_dir = "petridish_results/densities/"
quantile = 0.95

simulations = [sim for sim in amf.load_simulations(results_dir + "seed_*.npz") if "status" in sim]
if not simulations:
    sys.exit("no run statistics found (files written before petridishsimauto.py stored them)")
ok = [sim for sim in simulations if str(sim["status"]) == "ok"]
stopped = [sim for sim in simulations if str(sim["status"]) != "ok"]

wait_ok = np.array([float(sim["wait_days"]) for sim in ok])
phases = ["wall_root", "wall_wait", "wall_hyphae", "wall_analysis"]
wall = {p: np.array([float(sim.get(p, 0.)) / 60 for sim in simulations]) for p in phases} # [min]
wall_total = sum(wall.values())

print(str(len(ok)) + " seeds reached t*, " + str(len(stopped)) + " were stopped: " + str([int(sim["seed"]) for sim in stopped]))
if len(wait_ok):
    print("simulated days until t*: median %.2f, %d%% quantile %.2f, max %.2f" % (np.median(wait_ok), quantile * 100, np.quantile(wait_ok, quantile), wait_ok.max()))
for p in phases:
    print("%-14s median %.1f min, max %.1f min" % (p, np.median(wall[p]), wall[p].max()))
print("%-14s median %.1f min, max %.1f min" % ("total", np.median(wall_total), wall_total.max()))

c_ok, c_stopped, c_ref = "#4c72b0", "#dd8452", "0.4"
fig, axes = plt.subplots(2, 2, figsize=(11, 8))

# (a) simulated days until t*
ax = axes[0, 0]
if len(wait_ok):
    ax.hist(wait_ok, bins=20, color=c_ok, edgecolor="white")
    q = np.quantile(wait_ok, quantile)
    ax.axvline(q, color=c_ref, ls="--", lw=1)
    ax.annotate("%d%%: %.1f d" % (quantile * 100, q), (q, ax.get_ylim()[1] * 0.95), xytext=(4, 0), textcoords="offset points", va="top", fontsize=9)
ax.set_xlabel("simulated days from end of root phase to t* [d]")
ax.set_ylabel("seeds")
ax.set_title("Waiting time until t* (%d reached, %d stopped)" % (len(ok), len(stopped)), fontsize=10)

# (b) hyphal length behind the barrier while waiting
ax = axes[0, 1]
for sim in ok:
    ax.plot(sim["crossed_t"] - float(sim["t_root_end"]), sim["crossed_len"], color=c_ok, lw=0.8, alpha=0.4)
for sim in stopped:
    ax.plot(sim["crossed_t"] - float(sim["t_root_end"]), sim["crossed_len"], color=c_stopped, lw=1.2)
ax.axhline(10, color=c_ref, ls="--", lw=1)
ax.plot([], [], color=c_ok, label="reached t*")
if stopped:
    ax.plot([], [], color=c_stopped, label="stopped")
ax.legend(frameon=False, fontsize=9)
ax.set_xlabel("simulated days from end of root phase [d]")
ax.set_ylabel("hyphal length behind the barrier [cm]")
ax.set_title("Growth behind the barrier until t*", fontsize=10)

# (c) wall clock time per phase, seeds sorted by total time
ax = axes[1, 0]
order = np.argsort(wall_total)
x = np.arange(len(order))
bottom = np.zeros(len(order))
colors = ["#c7c7c7", "#937860", "#8172b3", "#55a868"] # not the colors of "reached t*" / "stopped" used in the other panels
labels = ["root phase", "waiting for t*", "hyphae after t*", "analysis"]
for p, c, l in zip(phases, colors, labels):
    ax.bar(x, wall[p][order], bottom=bottom, color=c, width=0.9, label=l)
    bottom += wall[p][order]
ax.legend(frameon=False, fontsize=9)
ax.set_xticks([])
ax.set_xlabel("seeds, sorted by total wall clock time")
ax.set_ylabel("wall clock time [min]")
ax.set_title("Wall clock time per phase", fontsize=10)

# (d) does the wall clock time of waiting follow the simulated waiting time?
ax = axes[1, 1]
wait_all = np.array([float(sim["wait_days"]) for sim in simulations])
is_ok = np.array([str(sim["status"]) == "ok" for sim in simulations])
ax.scatter(wait_all[is_ok], wall["wall_wait"][is_ok], s=12, color=c_ok, label="reached t*")
if (~is_ok).any():
    ax.scatter(wait_all[~is_ok], wall["wall_wait"][~is_ok], s=12, color=c_stopped, label="stopped")
ax.legend(frameon=False, fontsize=9)
ax.set_xlabel("simulated waiting time [d]")
ax.set_ylabel("wall clock time of waiting [min]")
ax.set_title("Simulated vs. wall clock waiting time", fontsize=10)

for ax in axes.flat:
    ax.spines[["top", "right"]].set_visible(False)
plt.tight_layout()
plt.savefig("petridish_timing.png", dpi=200, bbox_inches="tight")
print("written petridish_timing.png")
