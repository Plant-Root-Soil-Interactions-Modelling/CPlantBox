"""
Summarize maximum length and radius for main, first-order, and second-order roots.

Estimation about the radius and maximum length were restricted to
  - First 3 laterals of 1st order
  - First 50 % of 2nd order laterals emerging from these 1st order laterals

Outputs:
- results_lmax_radius.csv
- terminal summary with applied selection rules
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


# ===== PATTERNS =====

MAIN_RE = re.compile(r"^(main(_\d+)?|root(_\d+)?)$", re.IGNORECASE)
LAT1_RE = re.compile(r"^(?:lat_(\d+)|lat_(\d+)_(\d+))$", re.IGNORECASE)
LAT2_RE = re.compile(r"^(?:lat_(\d+)_lat_(\d+)|latlat_(\d+)_(\d+))$", re.IGNORECASE)


# ===== IO =====

def read_csv_clean(path):
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df["root_name"] = df["root_name"].astype(str).str.strip()
    return df


def find_global_files(data_dir, recursive=False):
    it = data_dir.rglob("*.csv") if recursive else data_dir.glob("*.csv")
    return [p for p in sorted(it) if p.is_file() and not p.name.lower().endswith("_marks.csv")]


# ===== PARSING =====

def parse_lat1_n(name):
    m = LAT1_RE.match(name)
    if not m:
        return None
    if m.group(1):
        return int(m.group(1))
    if m.group(2):
        return int(m.group(2))
    return None


def parse_lat2_nm(name):
    m = LAT2_RE.match(name)
    if not m:
        return None
    if m.group(1) and m.group(2):
        return int(m.group(1)), int(m.group(2))
    if m.group(3) and m.group(4):
        return int(m.group(3)), int(m.group(4))
    return None


# ===== FILTER =====

def first_half_lat2(df):
    df = df.copy()
    nm = df["root_name"].map(parse_lat2_nm)
    df["n"] = nm.map(lambda t: t[0] if t else np.nan)
    df["m"] = nm.map(lambda t: t[1] if t else np.nan)

    if "insertion_position" in df.columns:
        df["_ins"] = pd.to_numeric(df["insertion_position"], errors="coerce")
    else:
        df["_ins"] = np.nan

    parts = []
    for _, g in df.groupby("n"):
        if g["_ins"].notna().any():
            g = g.sort_values("_ins")
        else:
            g = g.sort_values("m")

        take = int(np.ceil(0.5 * len(g)))
        parts.append(g.head(take))

    if parts:
        return pd.concat(parts, ignore_index=True)
    return df.iloc[0:0]


# ===== STATS =====

def summarize(x):
    x = pd.to_numeric(x, errors="coerce").dropna()
    n = len(x)
    mean = x.mean() if n else np.nan
    std = x.std(ddof=1) if n >= 2 else np.nan
    return n, mean, std


# ===== MAIN =====

def main():
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=script_dir / "data" / "UptSmartRoot")
    parser.add_argument("--out", type=Path, default=Path("results_lmax_radius.csv"))
    parser.add_argument("--recursive", action="store_true")
    args = parser.parse_args()

    files = find_global_files(args.data_dir, args.recursive)

    lmax = {"main": [], "lat1": [], "lat2": []}
    radius = {"main": [], "lat1": [], "lat2": []}

    for path in files:
        df = read_csv_clean(path)

        df["length"] = pd.to_numeric(df["length"], errors="coerce")
        df["diameter"] = pd.to_numeric(df["diameter"], errors="coerce")
        df["radius"] = df["diameter"] / 2.0

        main_df = df[df["root_name"].str.match(MAIN_RE)]
        lmax["main"].extend(main_df["length"].dropna().tolist())
        radius["main"].extend(main_df["radius"].dropna().tolist())

        lat1_df = df[df["root_name"].map(lambda s: parse_lat1_n(s) in [1, 2, 3])]
        lmax["lat1"].extend(lat1_df["length"].dropna().tolist())
        radius["lat1"].extend(lat1_df["radius"].dropna().tolist())

        lat2_df = df[df["root_name"].map(lambda s: (parse_lat2_nm(s) or (None, None))[0] in [1, 2, 3])]
        lat2_half = first_half_lat2(lat2_df)
        lmax["lat2"].extend(lat2_half["length"].dropna().tolist())
        radius["lat2"].extend(lat2_half["radius"].dropna().tolist())

    rows = []

    for group in ["main", "lat1", "lat2"]:
        n_l, mean_l, std_l = summarize(pd.Series(lmax[group]))
        rows.append({"group": group, "metric": "lmax", "n": n_l, "mean": mean_l, "std": std_l})

        n_r, mean_r, std_r = summarize(pd.Series(radius[group]))
        rows.append({"group": group, "metric": "radius", "n": n_r, "mean": mean_r, "std": std_r})

    out = pd.DataFrame(rows).sort_values(["group", "metric"])
    out.to_csv(args.out, index=False)

    print("\n===== RESULTS =====")
    print(out.to_string(index=False))

    print("\n===== OUTPUT =====")
    print(f"Results saved to: {args.out.resolve()}")
    print("Applied rules:")
    print(" - Estimation about the radius and maximum length were restricted to")
    print("   First 3 laterals of 1st order")
    print("   First 50 % of 2nd order laterals emerging from these 1st order laterals")
    print("Finished.")


if __name__ == "__main__":
    main()