import sys; sys.path.append("../.."); sys.path.append("../../src/")

import AMFAnalysis as amf

# Combines the per-seed results written by petridishsimauto.py into the radial profile over all seeds
results_dir = "petridish_results/densities/"

simulations = amf.load_simulations(results_dir + "seed_*.npz")
# seeds stopped before t* (status "not_crossed", see max_wait_days in petridishsimauto.py) or because the front stalled after t*
# (status "stalled", see min_growing_tips, stall_fraction) have no densities, older files have no status
stopped = [int(sim["seed"]) for sim in simulations if str(sim.get("status", "ok")) != "ok"]
simulations = [sim for sim in simulations if str(sim.get("status", "ok")) == "ok"]
# older files have densities per cm^3: with hyphae_layer (thin layer geometry) they are converted to mm^-2 (cm^-3 * layer [cm] = cm^-2 = mm^-2 * 100),
# without it (old geometry, densities per volume of the whole height) they are skipped
for sim in simulations:
    if "density_per" not in sim and "hyphae_layer" in sim:
        h = float(sim["hyphae_layer"])
        for key in ["tip_dens", "stopped_tip_dens"]:
            sim[key] = sim[key] * h / 100
        sim["length_densities"] = sim["length_densities"] * h / 10 # cm cm^-3 -> cm cm^-2 -> mm mm^-2
        sim["density_per"] = "mm^2"
old = [int(sim["seed"]) for sim in simulations if "density_per" not in sim]
if old:
    print("skipping " + str(len(old)) + " seeds of the old geometry (rerun them): " + str(old))
simulations = [sim for sim in simulations if "density_per" in sim]
if not simulations:
    sys.exit("no simulations with the current geometry found")
print(str(len(simulations)) + " simulations reached t*, " + str(len(stopped)) + " were stopped before: " + str(stopped))
location = amf.ringCenters(simulations[0]["location"])   # "location" stores the outer ring radii

# time since t* (hyphae behind the barrier reached 10 cm total length), the absolute times of different seeds do not overlap
time_stars = [float(sim["time_star"]) for sim in simulations]
print("time_star from " + str(min(time_stars)) + " to " + str(max(time_stars)) + " days")
for sim in simulations:
    sim["times"] = sim["times"] - sim["time_star"]

t_common, tip_dens_interp = amf.interpolate_to_common_time(simulations, "tip_dens")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [mm$^{-2}$]")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [mm$^{-2}$]", stat="mean", filename="TipDensOverAllSim_mean.png")
# spatial and temporal profiles of filament and tip density (Oyarte Galvez et al. 2025 Fig. 2b), arrival time t_n of each ring from a sigmoid fit
fits = amf.plot_spatial_temporal(simulations, location, filename="SpatialTemporalProfile.png")
print(str(len(fits)) + " rings with arrival time fit, t_n from " + str(round(min(f["t_n"] for f in fits), 1)) + " to " + str(round(max(f["t_n"] for f in fits), 1)) + " h after t*")
