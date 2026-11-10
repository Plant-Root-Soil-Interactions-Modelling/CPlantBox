import sys; sys.path.append("../.."); sys.path.append("../../src/")

import numpy as np
import AMFAnalysis as amf
import AMFIO as amfio
import AMFPlots as amfp
import AMFSimulationSetup as amfs

# Combines the per-seed results written by petridishsimauto.py into the radial profile over all seeds
out = amfio.petridishOutputDirs() # below $PETRIDISH_OUT (cluster) or petridish_results

simulations = amfio.load_petridish_results(out["densities"] + "seed_*.npz") # only seeds that can be combined, times in days since t*
if not simulations:
    sys.exit("no simulations with the current geometry found")
location = amfs.ringCenters(simulations[0]["location"])   # "location" stores the outer ring radii

time_stars = [float(sim["time_star"]) for sim in simulations]
print("time_star from " + str(min(time_stars)) + " to " + str(max(time_stars)) + " days")

t_common, tip_dens_interp = amf.interpolate_to_common_time(simulations, "tip_dens")
amfp.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [mm$^{-2}$]", filename=out["plots"] + "TipDensOverAllSim.png")
amfp.plot_radial_profiles(location, t_common, tip_dens_interp, indices=[0, 30, 60, 90], ylabel="Tip density [mm$^{-2}$]", stat="mean", filename=out["plots"] + "TipDensOverAllSim_mean.png")
# spatial and temporal profiles of filament and tip density (Oyarte Galvez et al. 2025 Fig. 2b), arrival time t_n of each ring from a sigmoid fit
fits = amfp.plot_spatial_temporal(simulations, location, filename=out["plots"] + "SpatialTemporalProfile.png")
print(str(len(fits)) + " rings with arrival time fit, t_n from " + str(round(min(f["t_n"] for f in fits), 1)) + " to " + str(round(max(f["t_n"] for f in fits), 1)) + " h after t*")
# Runner vs BAS, anastomoses and encounter fraction (only seeds written by the current petridishsimauto.py)
comp = amfp.plot_hyphal_length_composition(simulations, location, filename=out["plots"] + "HyphalLengthComposition.png", target_fraction=0.3)
if comp is not None:
    print("BAS share of the hyphal length at the end: median " + str(round(float(np.nanmedian(comp["frac"][:, -1])), 3)) + " (experiment 0.3)")
amfp.plot_anastomosis_frequency(simulations, filename=out["plots"] + "AnastomosisFrequency.png", bin_h=6)
C = amfp.plot_encounter_fraction(simulations, filename=out["plots"] + "EncounterFraction.png", min_tips=5)
if C:
    print("encounter fraction f = C rho: " + ", ".join(f"{k} C = {v:.3e} mm^2 um^-1 h^-1" for k, v in C.items()))
