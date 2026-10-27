import sys; sys.path.append("../.."); sys.path.append("../../src/")

import AMFAnalysis as amf

# Combines the per-seed results written by petridishsimauto.py into the radial profile over all seeds
results_dir = "petridish_results/densities/"

simulations = amf.load_simulations(results_dir + "seed_*.npz")
# seeds stopped before t* (see max_wait_days in petridishsimauto.py) have no densities, older files have no status
stopped = [int(sim["seed"]) for sim in simulations if str(sim.get("status", "ok")) != "ok"]
simulations = [sim for sim in simulations if str(sim.get("status", "ok")) == "ok"]
print(str(len(simulations)) + " simulations reached t*, " + str(len(stopped)) + " were stopped before: " + str(stopped))
location = amf.ringCenters(simulations[0]["location"])   # "location" stores the outer ring radii

# time since t* (hyphae behind the barrier reached 10 cm total length), the absolute times of different seeds do not overlap
time_stars = [float(sim["time_star"]) for sim in simulations]
print("time_star from " + str(min(time_stars)) + " to " + str(max(time_stars)) + " days")
for sim in simulations:
    sim["times"] = sim["times"] - sim["time_star"]

t_common, tip_dens_interp = amf.interpolate_to_common_time(simulations, "tip_dens")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [cm$^{-3}$]")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [cm$^{-3}$]", stat="mean", filename="TipDensOverAllSim_mean.png")
