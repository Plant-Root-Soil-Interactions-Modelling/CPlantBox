import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import time
from sklearn.linear_model import LinearRegression
from scipy.interpolate import interp1d

# SegmentAnalyser with the mycorrhizal data, encounter fraction vs. density, combining simulations, carbon cost
# (geometry and simulation: AMFSimulationSetup, densities per ring and arrival fits: AMFDensities, network: AMFNetwork,
#  plots: AMFPlots, saving and loading: AMFIO, old unused functions: AMFLegacy)

def getMycSegmentAnalyser(plant, write=[],filename = "write",step = 0, stdwrite = False,):
        ana = pb.SegmentAnalyser(plant)
        ana.addData("colonization", plant.getNodeColonizations(2))
        ana.addData("colonizationTime", plant.getNodeColonizationTime(2))
        ana.addData("anastomosis", plant.getAnastomosisPoints(5))
        ana.addData("nodeTips", plant.getNodeTips(5))
        # ana.crop(crop)
        if write:
            ana.write(filename + "_" + str(step) + "_" +time.strftime("%Y%m%d-%H%M%S") + ".vtp",write)
        if stdwrite:
            ana.write(filename + "_" + str(step) + "_" +time.strftime("%Y%m%d-%H%M%S") + ".vtp",["radius","creationTime","subType","organType","colonization","colonizationTime","anastomosis","nodeTips"])
        return ana

def prepare_data(lenMat, vol, anaMat, tipMat):
    rho = lenMat / np.reshape(vol, (-1, 1)) # vol: one value or one per ring/box
    # print(rho)
    f_enc = np.divide(
        anaMat, tipMat,
        out=np.zeros_like(anaMat, dtype=float),
        where= (tipMat != 0)
    )
    n_locations, n_times = rho.shape
    rho = rho.ravel()
    f_enc = f_enc.ravel()
    # Zeile = Ort, Spalte = Zeit
    location = np.repeat(np.arange(n_locations), n_times)
    time = np.tile(np.arange(n_times), n_locations)
    valid = np.isfinite(rho) & np.isfinite(f_enc)
    return rho[valid], f_enc[valid], location[valid], time[valid]


def find_linear_range(rho, f_enc, window, step, min_points=10):
    results = []
    starts = np.arange(rho.min(), rho.max() - window, step)
    for start in starts:
        end = start + window
        mask = (rho >= start) & (rho <= end)
        if mask.sum() < min_points:
            continue
        r = rho[mask]
        f = f_enc[mask]
        model = LinearRegression(fit_intercept=False).fit(r.reshape(-1, 1), f)
        results.append({
            "rho_min": start,
            "rho_max": end,
            "rho_center": (start + end) / 2,
            "C": model.coef_[0],
            "R2": model.score(r.reshape(-1, 1), f),
            "n": mask.sum()
        })
    return results

def analyze_simulations(simulations, window, step, min_points=10):
    all_results = {}
    for sim in simulations:
        rho, f_enc, location, time = prepare_data(
            sim["lenMat"], sim["vol"], sim["anaMat"], sim["tipMat"]
        )
        all_results[sim["label"]] = {
            "rho": rho,
            "f_enc": f_enc,
            "location": location,
            "time": time,
            "linear": find_linear_range(
                rho, f_enc, window, step, min_points
            )
        }
    return all_results


def interpolate_to_common_time(simulations, key, n_times=99):
    # interpolates sim[key] (rings x times) of all simulations onto the time range they share

    t_start = max(sim["times"][0] for sim in simulations)
    t_end = min(sim["times"][-1] for sim in simulations)
    if t_start >= t_end:
        raise ValueError(f"interpolate_to_common_time: time ranges of the simulations do not overlap (latest start {t_start}, earliest end {t_end})")
    t_common = np.linspace(t_start, t_end, n_times)
    interp = np.stack([interp1d(sim["times"], sim[key], axis=1)(t_common) for sim in simulations])
    return t_common, interp   # (n_sim, n_r, n_t)


def estimate_carbon_cost_per_Org(length,radius,pdry,pcarbon,dcell,CUE):
    gamma = np.pi * radius *radius*dcell*pdry*pcarbon
    gamma = gamma / CUE
    carbon_cost = gamma*length
    return carbon_cost
