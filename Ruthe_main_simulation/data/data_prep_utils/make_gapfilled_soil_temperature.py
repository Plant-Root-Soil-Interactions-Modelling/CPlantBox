"""Gap-filled daily soil temperature profiles (10-200 cm) per growing season.

Inputs
  - Ruthe climate files (Ruthe_95.xls, Ruthe_96.xls, Ruthe_97.xls), sheet 'Stundenwerte' with the hourly
    soil temperatures Bodentemp1-6 (Bodentemp1-5 = 10-50 cm, Bodentemp6 = second sensor at 10 cm).
  - A HYDRUS-1D project folder with Obs_Node.out (hourly output) and PROFILE.DAT (to get the depth of each
    observation node).

Processing
  1. Measured: hourly values outside -20...40 C (incl. the 999 error code) are discarded; a daily mean is
     only formed from at least 20 valid hours. Gaps in Bodentemp1 are filled with Bodentemp6.
  2. HYDRUS: daily mean of the hourly Obs_Node.out output. HYDRUS time t in (n-1, n] is day n,
     day 1 = 01.04.1994 (start of the simulation).
  3. Per growing season: measured value where available (10-50 cm), HYDRUS everywhere else,
     i.e. measurement gaps and all depths below 50 cm.

Output: soil_temperature_<season>.csv, one row per day, columns:
  date, 10, 20, ..., 200 (C, 2 decimals), source_10 ... source_50 (Bodentemp1-6 or 'HYDRUS').
"""
from pathlib import Path

import numpy as np
import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]  # the data folder

RAW_TEMP_DIR = BASE_DIR / "data_raw"
HYDRUS_DIR = BASE_DIR / "hydrus" / "soil_temperature_simulation"
OUT_DIR = BASE_DIR / "data_processed"


YEARS = (1995, 1996, 1997)
SIM_START = pd.Timestamp("1994-04-01")          # date of HYDRUS day 1
PERIODS = {"95": ("1994-11-06", "1995-08-09"),  # sowing - harvest
           "96": ("1995-11-02", "1996-08-27"),
           "97": ("1996-11-05", "1997-08-12")}
MEASURED = {10: "Bodentemp1", 20: "Bodentemp2", 30: "Bodentemp3", 40: "Bodentemp4", 50: "Bodentemp5"}
SENSORS = [f"Bodentemp{i}" for i in range(1, 7)]
T_VALID = (-20, 40)                             # plausible range, C (exclusive)
MIN_HOURS = 20                                  # valid hours needed for a daily mean


def read_measured(raw_dir, years):
    """Daily mean soil temperature Bodentemp1-6 from the hourly sheets of the Ruthe climate files."""
    hourly = []
    for y in years:
        f = Path(raw_dir) / f"Ruthe_{str(y)[2:]}.xls"
        sheet = next(s for s in pd.ExcelFile(f).sheet_names if s.strip() == "Stundenwerte")
        h = pd.read_excel(f, sheet_name=sheet, header=None)
        cols = [str(c).replace("\n", " ").replace(" ", "") for c in h.iloc[1]]
        h = h.iloc[2:]
        h.columns = cols
        h = h[pd.to_numeric(h["Jahr"], errors="coerce") == y]
        d = pd.DataFrame({"date": pd.Timestamp(f"{y}-01-01") + pd.to_timedelta(h["Tag"].astype(int) - 1, unit="D")})
        for s in SENSORS:
            d[s] = pd.to_numeric(h[next(c for c in cols if c.startswith(s))], errors="coerce").values
        hourly.append(d)
    hourly = pd.concat(hourly, ignore_index=True)
    valid = hourly[SENSORS].where((hourly[SENSORS] > T_VALID[0]) & (hourly[SENSORS] < T_VALID[1]))
    g = valid.groupby(hourly["date"])
    daily = g.mean().where(g.count() >= MIN_HOURS)

    # 10 cm: Bodentemp1, gaps filled with the second 10 cm sensor Bodentemp6
    daily["source_10"] = None
    daily.loc[daily.Bodentemp1.notna(), "source_10"] = "Bodentemp1"
    fill = daily.Bodentemp1.isna() & daily.Bodentemp6.notna()
    daily.loc[fill, "source_10"] = "Bodentemp6"
    daily.loc[fill, "Bodentemp1"] = daily.loc[fill, "Bodentemp6"]
    return daily


def observation_depths(profile_dat):
    """Depths (cm, positive down) of the observation nodes listed at the end of PROFILE.DAT."""
    lines = Path(profile_dat).read_text().split("\n")
    node_depth, last = {}, 0
    for i, l in enumerate(lines):
        f = l.split()
        if len(f) == 10 and f[0].isdigit():                     # node table: node, x, h, Mat, Lay, ..., Temp
            node_depth[int(f[0])] = -float(f[1])
            last = i
    tok = " ".join(lines[last + 1:]).split()
    nodes = [int(t) for t in tok[1:1 + int(tok[0])]]
    return [round(node_depth[n]) for n in nodes]


def read_hydrus(hydrus_dir):
    """Daily mean simulated soil temperature at the observation nodes; columns = depth in cm."""
    depths = observation_depths(Path(hydrus_dir) / "PROFILE.DAT")
    rows = []
    with open(Path(hydrus_dir) / "Obs_Node.out") as obs_file:
        for line in obs_file:
            f = line.split()
            if len(f) == 1 + 3 * len(depths):                  # time, then h, theta, Temp per node
                try:
                    rows.append([float(x) for x in f])
                except ValueError:
                    pass
    a = np.array(rows)
    sim = pd.DataFrame(a[:, 3::3], index=a[:, 0], columns=depths)
    day = np.ceil(sim.index - 1e-6).astype(int)                 # (n-1, n] -> day n
    daily = sim.groupby(day).mean()
    daily.index = SIM_START + pd.to_timedelta(daily.index - 1, unit="D")
    return daily


def gap_filled(measured, sim, start, end):
    """Measured daily means where available, HYDRUS elsewhere, for one growing season."""
    idx = pd.date_range(start, end, freq="D", name="date")
    s = sim.reindex(idx)
    assert s.notna().all().all(), "HYDRUS output does not cover the whole period"
    out, src = s.copy(), pd.DataFrame(index=idx)
    for z, c in MEASURED.items():
        m = measured[c].reindex(idx)
        out[z] = m.fillna(s[z])
        src[f"source_{z}"] = "HYDRUS"
        if z == 10:
            src.loc[m.notna(), "source_10"] = measured.source_10.reindex(idx)[m.notna()]
        else:
            src.loc[m.notna(), f"source_{z}"] = c
    out.columns = [str(z) for z in out.columns]
    return out.round(2).join(src)


if __name__ == "__main__":
    measured = read_measured(RAW_TEMP_DIR, YEARS)
    sim = read_hydrus(HYDRUS_DIR)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for season, (start, end) in PERIODS.items():
        df = gap_filled(measured, sim, start, end)
        fn = OUT_DIR / f"soil_temperature_{season}.csv"
        df.to_csv(fn, date_format="%Y-%m-%d")
        n_meas = {z: int((df[f"source_{z}"] != "HYDRUS").sum()) for z in MEASURED}
        print(f"{fn.name}: {start} - {end}, {len(df)} days, measured days per depth {n_meas}")
