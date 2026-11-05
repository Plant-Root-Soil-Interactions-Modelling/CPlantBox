"""
Compute lateral position relative to the matched main root and create
two scatter plots separated by E and F groups.

Outputs:
- results_lat_origin_nearest_main_pos.csv
- Distance_from_base_vs_Length_by_EF.png
- Distance_from_tip_vs_Length_by_EF.png
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# ===== SETTINGS =====

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 11,
    "figure.titlesize": 20
})

MAIN_RE = re.compile(r"^(main(_\d+)?|root(_\d+)?)$", re.IGNORECASE)
LAT1_RE = re.compile(r"^(?:lat_(\d+)|lat_(\d+)_(\d+))$", re.IGNORECASE)


# ===== IO =====

def read_csv_clean(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]

    for col in [
        "root_name", "mark_type", "value", "image", "source",
        "parent_name", "root", "root_ontology"
    ]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip()

    return df


def find_marks_files(data_dir: Path, recursive: bool = False) -> list[Path]:
    if not data_dir.exists():
        raise FileNotFoundError(f"Data dir does not exist: {data_dir}")

    it = data_dir.rglob("*_marks.csv") if recursive else data_dir.glob("*_marks.csv")
    marks_files = [p for p in sorted(it) if p.is_file()]

    print(f"[INFO] Reading from: {data_dir.resolve()}")
    print(f"[INFO] recursive={recursive}")
    print(f"[INFO] Found {len(marks_files)} marks CSV files.")
    for p in marks_files[:20]:
        print(f"       marks: {p.relative_to(data_dir)}")
    if len(marks_files) > 20:
        print(f"       ... ({len(marks_files) - 20} more)")

    return marks_files


def get_image_name(marks_df: pd.DataFrame) -> str | None:
    for col in ["image", "source"]:
        if col in marks_df.columns:
            s = marks_df[col].dropna()
            if len(s) > 0:
                value = str(s.iloc[0]).strip()
                return value if value else None
    return None


def get_base_file_from_marks(marks_path: Path) -> Path:
    name = marks_path.name
    if name.lower().endswith("_marks.csv"):
        return marks_path.with_name(name[:-10] + ".csv")
    return marks_path


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


# ===== GROUPS =====

def classify_image_group(image_name: str) -> str:
    if pd.isna(image_name):
        return "Other"

    s = str(image_name).strip()
    if re.search(r"\bE\.\d+\b", s, flags=re.IGNORECASE):
        return "E"
    if re.search(r"\bF\.\d+\b", s, flags=re.IGNORECASE):
        return "F"
    return "Other"


# ===== DIAMETER =====

def build_diameter_map(base_df: pd.DataFrame) -> dict:
    if "diameter" not in base_df.columns:
        return {}

    root_col = None
    for candidate in ["root_name", "root"]:
        if candidate in base_df.columns:
            root_col = candidate
            break

    if root_col is None:
        return {}

    tmp = base_df[[root_col, "diameter"]].copy()
    tmp[root_col] = tmp[root_col].astype(str).str.strip()
    tmp["diameter"] = to_num(tmp["diameter"])
    tmp = tmp.dropna(subset=[root_col, "diameter"])

    if len(tmp) == 0:
        return {}

    return tmp.groupby(root_col)["diameter"].first().to_dict()


# ===== ANALYSIS =====

def compute_lat_origin_to_main_position(marks_path: Path, main_mark_types: list[str]) -> pd.DataFrame:
    marks = read_csv_clean(marks_path)
    image_name = get_image_name(marks)

    base_path = get_base_file_from_marks(marks_path)
    diameter_map = {}

    if base_path.exists():
        try:
            base_df = read_csv_clean(base_path)
            diameter_map = build_diameter_map(base_df)
        except Exception as e:
            print(f"[WARN] Could not read base file {base_path.name}: {e}")
    else:
        print(f"[WARN] No matching base file found for {marks_path.name}")

    required = {"root_name", "mark_type", "x", "y", "position_from_base"}
    if not required.issubset(marks.columns):
        return pd.DataFrame(columns=[
            "image", "main", "lat_root",
            "position_from_base", "position_from_tip",
            "length", "diameter"
        ])

    main_len = marks[
        marks["root_name"].apply(lambda x: bool(MAIN_RE.match(str(x))))
        & (marks["mark_type"].astype(str).str.lower() == "length")
    ].copy()

    main_len["position_from_base"] = to_num(main_len["position_from_base"])
    main_len = main_len.dropna(subset=["position_from_base"])

    main_len_map = (
        main_len.groupby(main_len["root_name"].astype(str))["position_from_base"].max().to_dict()
        if len(main_len) > 0 else {}
    )

    mm_types = {t.strip().lower() for t in main_mark_types}
    main_pts = marks[
        marks["root_name"].apply(lambda x: bool(MAIN_RE.match(str(x))))
        & marks["mark_type"].astype(str).str.lower().isin(mm_types)
    ].copy()

    main_pts["x"] = to_num(main_pts["x"])
    main_pts["y"] = to_num(main_pts["y"])
    main_pts["position_from_base"] = to_num(main_pts["position_from_base"])
    main_pts = main_pts.dropna(subset=["x", "y", "position_from_base"])

    lat_len = marks[
        marks["root_name"].apply(lambda x: bool(LAT1_RE.match(str(x))))
        & (marks["mark_type"].astype(str).str.lower() == "length")
    ].copy()

    lat_len["position_from_base"] = to_num(lat_len["position_from_base"])
    lat_len = lat_len.dropna(subset=["position_from_base"])

    lat_len_map = (
        lat_len.groupby(lat_len["root_name"].astype(str))["position_from_base"].first().to_dict()
        if len(lat_len) > 0 else {}
    )

    lat_org = marks[
        marks["root_name"].apply(lambda x: bool(LAT1_RE.match(str(x))))
        & (marks["mark_type"].astype(str).str.lower() == "origin")
    ].copy()

    lat_org["x"] = to_num(lat_org["x"])
    lat_org["y"] = to_num(lat_org["y"])
    lat_org = lat_org.dropna(subset=["x", "y"])

    if len(lat_org) == 0 or len(main_pts) == 0:
        return pd.DataFrame(columns=[
            "image", "main", "lat_root",
            "position_from_base", "position_from_tip",
            "length", "diameter"
        ])

    main_xy = main_pts[["x", "y"]].to_numpy(dtype=float)

    rows = []
    for _, lr in lat_org.iterrows():
        lx, ly = float(lr["x"]), float(lr["y"])
        d2 = (main_xy[:, 0] - lx) ** 2 + (main_xy[:, 1] - ly) ** 2
        j = int(np.argmin(d2))
        best = main_pts.iloc[j]

        main_root = str(best["root_name"])
        main_length = main_len_map.get(main_root, np.nan)

        lat_root = str(lr["root_name"])
        pos_base = float(best["position_from_base"])
        pos_tip = float(main_length) - pos_base if pd.notna(main_length) else np.nan
        length_val = lat_len_map.get(lat_root, np.nan)
        diameter_val = diameter_map.get(lat_root, np.nan)

        rows.append({
            "image": image_name,
            "main": main_root,
            "lat_root": lat_root,
            "position_from_base": pos_base,
            "position_from_tip": pos_tip,
            "length": float(length_val) if pd.notna(length_val) else np.nan,
            "diameter": float(diameter_val) if pd.notna(diameter_val) else np.nan,
        })

    return pd.DataFrame(rows)


# ===== PLOTS =====

def save_plot(fig, plot_dir: Path, filename: str):
    plot_dir.mkdir(exist_ok=True)
    plot_path = plot_dir / filename
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"Plot saved: {plot_path.name}")


def plot_by_ef(df: pd.DataFrame, plot_dir: Path, x_col: str):
    df = df.copy()
    df["EF_group"] = df["image"].apply(classify_image_group)

    styles = {
        "E": {"color": "blue", "marker": "o"},
        "F": {"color": "green", "marker": "^"},
    }

    if x_col == "position_from_tip":
        xlabel = "Position from Tip [cm]"
        filename = "Distance_from_tip_vs_Length_by_EF.png"
    elif x_col == "position_from_base":
        xlabel = "Position from Base [cm]"
        filename = "Distance_from_base_vs_Length_by_EF.png"
    else:
        raise ValueError(f"Unknown x column: {x_col}")

    fig, ax = plt.subplots(figsize=(9, 6))

    for grp in ["E", "F"]:
        sub = df[df["EF_group"] == grp]
        if len(sub) == 0:
            continue
        ax.scatter(
            sub[x_col],
            sub["length"],
            label=grp,
            color=styles[grp]["color"],
            marker=styles[grp]["marker"],
            s=70,
            alpha=0.85
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Length [cm]")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), borderaxespad=0.0)

    fig.tight_layout()
    save_plot(fig, plot_dir, filename)
    plt.close(fig)


# ===== MAIN =====

def main():
    script_dir = Path(__file__).resolve().parent
    default_data_dir = script_dir / "data" / "UptSmartRoot"
    default_out = script_dir / "results_lat_origin_nearest_main_pos.csv"

    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=default_data_dir)
    ap.add_argument("--out", type=Path, default=default_out)
    ap.add_argument("--recursive", action="store_true")
    ap.add_argument(
        "--main-mark-types",
        nargs="+",
        default=["number", "measure", "free text"]
    )
    args = ap.parse_args()

    marks_files = find_marks_files(args.data_dir, recursive=args.recursive)

    all_rows = []
    for marks_path in marks_files:
        try:
            all_rows.append(compute_lat_origin_to_main_position(marks_path, args.main_mark_types))
        except Exception as e:
            print(f"[ERROR] {marks_path.name}: {e}")

    if not all_rows:
        print("[WARN] No results produced.")
        return

    out = pd.concat(all_rows, ignore_index=True)
    out = out.sort_values(["image", "main", "lat_root"], na_position="last").reset_index(drop=True)
    out = out[[
        "image", "main", "lat_root",
        "position_from_base", "position_from_tip",
        "length", "diameter"
    ]]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)

    n_E = 0
    n_F = 0

    df = out.dropna(subset=["position_from_tip", "position_from_base", "length"]).copy()

    if len(df) > 0:
        df["EF_group"] = df["image"].apply(classify_image_group)
        n_E = int((df["EF_group"] == "E").sum())
        n_F = int((df["EF_group"] == "F").sum())

        plot_dir = script_dir / "plots"
        plot_by_ef(df, plot_dir, "position_from_base")
        plot_by_ef(df, plot_dir, "position_from_tip")

    print("\n===== OUTPUT =====")
    print(f"Results saved to: {args.out.resolve()}")
    print(f"n(E) = {n_E}")
    print(f"n(F) = {n_F}")
    print("Finished.")


if __name__ == "__main__":
    main()