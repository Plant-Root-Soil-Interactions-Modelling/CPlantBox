import sys; sys.path.append("../.."); sys.path.append("../../src/")

import AMFAnalysis as amf

# Combines the per-seed results written by petridishsimauto.py into the radial profile over all seeds
results_dir = "petridish_results/densities/"

simulations = amf.load_simulations(results_dir + "seed_*.npz")
location = simulations[0]["location"]

t_common, tip_dens_interp = amf.interpolate_to_common_time(simulations, "tip_dens")
amf.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density")
