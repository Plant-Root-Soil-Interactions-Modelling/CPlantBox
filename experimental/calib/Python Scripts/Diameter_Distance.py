"""
Plot radius vs distance separated by E and F.

Outputs:
- two plots (base / tip)
- sample sizes
- mean and std of radius (total, E, F)
"""

from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


# ===== SETTINGS =====

CSV_PATH = Path("results_lat_origin_nearest_main_pos.csv")
PLOTS_DIR = Path("plots_radius")
PLOTS_DIR.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 11
})


# ===== GROUP =====

def classify_ef(image):
    s = str(image)
    if "Run4 E." in s:
        return "E"
    if "Run4 F." in s:
        return "F"
    return "Other"


# ===== PLOT =====

def save_plot(fig, name):
    fig.savefig(PLOTS_DIR / name, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Plot saved: {name}")


def plot_ef(df, xcol, name, xlabel, mean_E, mean_F):

    styles = {
        "E": {"color": "blue", "marker": "o"},
        "F": {"color": "green", "marker": "^"}
    }

    fig, ax = plt.subplots(figsize=(9, 6))

    for g in ["E", "F"]:
        sub = df[df["EF"] == g]
        if len(sub) == 0:
            continue

        ax.scatter(
            sub[xcol],
            sub["radius"],
            color=styles[g]["color"],
            marker=styles[g]["marker"],
            s=70,
            alpha=0.85,
            label=g
        )

    ax.axhline(mean_E, color="blue", lw=1.8, linestyle="--", label="Mean E")
    ax.axhline(mean_F, color="green", lw=1.8, linestyle="--", label="Mean F")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Radius [cm]")

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()

    save_plot(fig, name)


# ===== MAIN =====

def main():

    df = pd.read_csv(CSV_PATH)

    for c in ["position_from_base", "position_from_tip", "diameter"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df.dropna(subset=["position_from_base", "position_from_tip", "diameter"]).copy()

    df["radius"] = df["diameter"] / 2.0
    df["EF"] = df["image"].apply(classify_ef)
    df = df[df["EF"].isin(["E", "F"])]

    n_total = len(df)
    n_E = int((df["EF"] == "E").sum())
    n_F = int((df["EF"] == "F").sum())

    mean_all = df["radius"].mean()
    std_all = df["radius"].std()

    mean_E = df[df["EF"] == "E"]["radius"].mean()
    std_E = df[df["EF"] == "E"]["radius"].std()

    mean_F = df[df["EF"] == "F"]["radius"].mean()
    std_F = df[df["EF"] == "F"]["radius"].std()

    plot_ef(
        df,
        "position_from_base",
        "Radius_vs_base_EF.png",
        "Position from Root Base [cm]",
        mean_E,
        mean_F
    )

    plot_ef(
        df,
        "position_from_tip",
        "Radius_vs_tip_EF.png",
        "Position from Root Tip [cm]",
        mean_E,
        mean_F
    )

    print("\n===== OUTPUT =====")
    print(f"n(total) = {n_total}")
    print(f"n(E) = {n_E}")
    print(f"n(F) = {n_F}")
    print(f"\nMean radius total = {mean_all:.4f} ± {std_all:.4f}")
    print(f"Mean radius E = {mean_E:.4f} ± {std_E:.4f}")
    print(f"Mean radius F = {mean_F:.4f} ± {std_F:.4f}")
    print("Finished.")


if __name__ == "__main__":
    main()