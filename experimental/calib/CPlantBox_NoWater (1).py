import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import plantbox as pb

plt.rcParams.update({
    'font.size': 14,          # Standard Schrift
    'axes.titlesize': 16,     # Titel
    'axes.labelsize': 14,     # Achsenbeschriftung
    'xtick.labelsize': 13,    # x Tick Labels
    'ytick.labelsize': 13,    # y Tick Labels
    'legend.fontsize': 12,    # Legende
    'figure.titlesize': 20    # Figure Titel
})

plant = pb.Plant()

param_path = Path("../modelparameter/structural/rootsystem")
param_name = "wheat"
xml_file = param_path / f"{param_name}.xml"
plant.readParameters(str(xml_file),fromFile = True, verbose = False)

def get_or_create_root_param(plant: pb.Plant, sub_type: int) -> pb.RootRandomParameter:
    plist = plant.getOrganRandomParameter(pb.root)  # LISTE
    for p in plist:
        if hasattr(p, "subType") and int(p.subType) == int(sub_type):
            return p
    p_new = pb.RootRandomParameter(plant)
    p_new.subType = sub_type
    plant.setOrganRandomParameter(p_new)
    return p_new

p0 = get_or_create_root_param(plant, 1)  # taproot
p1 = get_or_create_root_param(plant, 2)  # lateral 1st
p2 = get_or_create_root_param(plant, 3)  # lateral 2nd

# --- Taproot (subType 1)
p0.name = "taproot"
p0.a = 0.023612          # [cm] 
p0.a_s = 0.0132
p0.lb = 0.342549        # [cm]
p0.lbs = 0.359216
p0.la = 3.075370        # [cm]
p0.las = 3.747081
p0.ln = 0.351840        # [cm]
p0.lns = 0.197011
p0.lmax = 112#83.0           # [cm]
p0.lmaxs = 112 / 10. #83.0           # [cm]
## Linear Growth
p0.r = 2.197             # [cm/day]
p0.gf = pb.GrowthFunctionType.linear
## Sigmoid Growth
#p0.r = 6.349
#p0.gf = pb.GrowthFunctionType.sigmoid
print(p0)
raise Exception
p0.successor = [[2]]
p0.successorP = [[1.0]]

# --- Lateral 1. Ordnung (subType 2)
p1.name = "lateral_1st"
p1.a = 0.0097
p1.a_s = 0.0023
p1.lb = 0.7765
p1.lbs = 0.7049
p1.la = 1.261800
p1.las = 1.382121
p1.ln = 0.422848
p1.lns = 0.248059
p1.lmax = 2.809707
p1.lmaxs = 3.101439
p1.r = 0.244
p1.gf = pb.GrowthFunctionType.linear
p1.successor = [[3]]
p1.successorP = [[1.0]]

# --- Lateral 2. Ordnung (subType 3)
p2.name = "lateral_2nd"
p2.a = 0.007442
p2.a_s = 0.0015
p2.r = 0.224
p2.lmax = 0.334268
p2.lmaxs = 0.376370
p2.gf = pb.GrowthFunctionType.linear

srp = pb.SeedRandomParameter(plant)
srp.seedPos = pb.Vector3d(0.0, 0.0, -3.0)  # [cm]

# Basal Root Parameter
srp.maxB   = 16.33
srp.firstB = 11.46
srp.delayB = 2.09
plant.setOrganRandomParameter(srp)

# Experiment Setup
simtime_days = 48.0
cyl_d_cm = 11.0
cyl_r_cm = cyl_d_cm / 2.0
cyl_depth_cm = 80.0

cylinder = pb.SDF_PlantContainer(cyl_r_cm, cyl_r_cm, cyl_depth_cm, False)
plant.setGeometry(cylinder)

# Simulate
plant.initialize()
plant.simulate(simtime_days)

print("Done. Roots:", len(plant.getPolylines()))


# -----------------------------
# Root-System Plot (2D)
# -----------------------------
def plot_root_system_2d(plant, stride=1, lw=0.8, depth=80):
    polylines = plant.getPolylines()

    # -------- Seitenansicht (X-Z) --------
    fig_xz = plt.figure(figsize=(6, 10))
    ax_xz = fig_xz.add_subplot(111)

    for i, pl in enumerate(polylines):
        if stride > 1 and (i % stride != 0):
            continue

        x = np.array([p.x for p in pl])
        z = np.array([p.z for p in pl])

        depth_vals = -z                # 0..80 (positiv nach unten)
        ax_xz.plot(x, depth_vals, linewidth=lw)

    ax_xz.set_title("Seitenansicht (X–Z)")
    ax_xz.set_xlabel("X [cm]")
    ax_xz.set_ylabel("Depth [cm]")

    # Harter Flip: 0 oben, 80 unten
    ax_xz.set_ylim(depth + 1, 0)

    ax_xz.set_aspect("equal", adjustable="box")
    ax_xz.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()


    # -------- Draufsicht (X-Y) --------
    fig_xy = plt.figure(figsize=(8,8))
    ax_xy = fig_xy.add_subplot(111)

    for i, pl in enumerate(polylines):
        if stride > 1 and (i % stride != 0):
            continue

        x = np.array([p.x for p in pl])
        y = np.array([p.y for p in pl])

        ax_xy.plot(x, y, linewidth=lw)

    ax_xy.set_title("Draufsicht (X–Y)")
    ax_xy.set_xlabel("X [cm]")
    ax_xy.set_ylabel("Y [cm]")

    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

plot_root_system_2d(plant, stride=1, lw=0.8)

# -----------------------------
# RLD (Root Length Density) Plot über Tiefe, getrennt nach Zeitpunkten
# -----------------------------
def plot_rld_over_time(
    plant: pb.Plant,
    times_days,
    depth_cm: float,
    radius_cm: float,
    layers: int = 32,
    out_path: Path = Path("plots/rld_over_time.png"),
):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Volumen pro Tiefenlayer im Zylinder
    layer_thickness = depth_cm / layers
    layer_volume = layer_thickness * (np.pi * radius_cm * radius_cm)  # [cm^3]

    z_ = np.linspace(0.0, -depth_cm, layers) 

    base_ana = pb.SegmentAnalyser(plant) 

    rld_all = []
    for t in times_days:
        ana = pb.SegmentAnalyser(base_ana)     
        ana.filter("creationTime", 0.0, float(t)) 
        ana.pack()

        distrib = ana.distribution("length", 0.0, -depth_cm, layers, True)  # Länge pro Layer
        rld = np.array(distrib) / layer_volume  # [cm/cm^3]
        rld_all.append(rld)

    rld_all = np.asarray(rld_all)
    x_max = max(1e-6, float(np.nanmax(rld_all)))  

    # 10 Plots: 2 Reihen x 5 Spalten
    n = len(times_days)
    nrows, ncols = 2, 5
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(14, 6))
    axes = axes.flatten()

    for i, t in enumerate(times_days):
        ax = axes[i]
        ax.plot(rld_all[i], z_)
        ax.set_title(f"t = {t:g} d")
        ax.set_xlim(0, x_max * 1.05)

    for j in range(n, nrows * ncols):
        axes[j].axis("off")

    plt.setp(axes[-ncols:], xlabel="RLD (cm/cm$^3$)")
    for ax in axes[::ncols]:
        ax.set_ylabel("Depth (cm)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()


def plot_rld_tap_basal_over_time(
    plant: pb.Plant,
    times_days,
    depth_cm: float,
    radius_cm: float,
    layers: int = 32,
    out_path: Path = Path("plots/rld_tap_basal_over_time.png"),
):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    layer_thickness = depth_cm / layers
    layer_volume = layer_thickness * (np.pi * radius_cm * radius_cm)

    z_ = np.linspace(0.0, -depth_cm, layers)

    base_ana = pb.SegmentAnalyser(plant)

    rld_all = []
    for t in times_days:
        # --- Taproot
        ana_tap = pb.SegmentAnalyser(base_ana)
        ana_tap.filter("creationTime", 0.0, float(t))
        ana_tap.filter("subType", 1, 1)
        ana_tap.pack()

        distrib_tap = np.array(
            ana_tap.distribution("length", 0.0, -depth_cm, layers, True)
        )

        # --- Basal roots
        ana_basal = pb.SegmentAnalyser(base_ana)
        ana_basal.filter("creationTime", 0.0, float(t))
        ana_basal.filter("subType", 4, 4)
        ana_basal.pack()

        distrib_basal = np.array(
            ana_basal.distribution("length", 0.0, -depth_cm, layers, True)
        )

        # --- Summe Tap + Basal
        distrib = distrib_tap + distrib_basal
        rld = distrib / layer_volume
        rld_all.append(rld)

    rld_all = np.asarray(rld_all)
    x_max = max(1e-6, float(np.nanmax(rld_all)))

    n = len(times_days)
    nrows, ncols = 2, 5
    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, sharex=True, sharey=True, figsize=(14, 6)
    )
    axes = axes.flatten()

    for i, t in enumerate(times_days):
        ax = axes[i]
        ax.plot(rld_all[i], z_)
        ax.set_title(f"t = {t:g} d")
        ax.set_xlim(0, x_max * 1.05)

    for j in range(n, nrows * ncols):
        axes[j].axis("off")

    plt.setp(axes[-ncols:], xlabel="RLD (cm/cm$^3$)")
    for ax in axes[::ncols]:
        ax.set_ylabel("Depth (cm)")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()


times_days = [5, 10, 15, 20, 25, 30, 35, 40, 45, 48]
plot_rld_over_time(
    plant,
    times_days=times_days,
    depth_cm=cyl_depth_cm,
    radius_cm=cyl_r_cm,
    layers=32,
    out_path=Path("plots/rld_over_time.png"),
)

plot_rld_tap_basal_over_time(
    plant,
    times_days=times_days,
    depth_cm=cyl_depth_cm,
    radius_cm=cyl_r_cm,
    layers=32,
    out_path=Path("plots/rld_tap_basal_over_time.png"),
)

def plot_root_surface_density_by_diameter(
    plant: pb.Plant,
    depth_slices_cm=((0, 7), (7, 73), (73, 80)),
    cyl_radius_cm: float = 5.5,
    time_day: float | None = None,
    diameter_bins_mm=np.arange(0.0, 0.75, 0.05),
    out_path: Path = Path("plots/root_surface_density_by_diameter.png"),
):
    out_path.parent.mkdir(parents=True, exist_ok=True)

    nodes = plant.getNodes()
    segments = plant.getSegments()

    ana = pb.SegmentAnalyser(plant)
    if time_day is not None:
        ana.filter("creationTime", 0.0, float(time_day))
        ana.pack()

    lengths_cm = np.asarray(ana.getParameter("length"), dtype=float)
    radii_cm   = np.asarray(ana.getParameter("radius"), dtype=float)

    creation_times = np.asarray(ana.getParameter("creationTime"), dtype=float)

    n_seg = min(len(segments), len(lengths_cm), len(radii_cm), len(creation_times))

    seg_surface_cm2 = []
    seg_depth_cm = []
    seg_diam_mm = []

    for i in range(n_seg):
        if time_day is not None and creation_times[i] > float(time_day):
            continue

        s = segments[i]
        n1 = nodes[s.x]
        n2 = nodes[s.y]

        zmid_cm = -0.5 * (n1.z + n2.z)

        r_cm = radii_cm[i]
        L_cm = lengths_cm[i]
        d_mm = 2.0 * r_cm * 10.0
        A_cm2 = 2.0 * np.pi * r_cm * L_cm

        seg_depth_cm.append(zmid_cm)
        seg_diam_mm.append(d_mm)
        seg_surface_cm2.append(A_cm2)

    seg_depth_cm = np.asarray(seg_depth_cm, dtype=float)
    seg_diam_mm = np.asarray(seg_diam_mm, dtype=float)
    seg_surface_cm2 = np.asarray(seg_surface_cm2, dtype=float)

    bin_left = diameter_bins_mm[:-1]
    bin_right = diameter_bins_mm[1:]
    class_labels = [f"{l:.2f}–{r:.2f}" for l, r in zip(bin_left, bin_right)]

    rsd_all = []
    ymax = 0.0

    for z_top, z_bottom in depth_slices_cm:
        mask_depth = (seg_depth_cm >= z_top) & (seg_depth_cm < z_bottom)

        slice_height_cm = z_bottom - z_top
        slice_volume_cm3 = np.pi * cyl_radius_cm**2 * slice_height_cm

        rsd_vals = []
        for d_lo, d_hi in zip(bin_left, bin_right):
            mask_diam = (seg_diam_mm >= d_lo) & (seg_diam_mm < d_hi)
            surface_sum = np.sum(seg_surface_cm2[mask_depth & mask_diam])
            rsd = surface_sum / slice_volume_cm3
            rsd_vals.append(rsd)

        rsd_vals = np.asarray(rsd_vals, dtype=float)
        rsd_all.append(rsd_vals)
        if len(rsd_vals) > 0:
            ymax = max(ymax, float(np.nanmax(rsd_vals)))

    # ---- Plot
    fig, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(14, 10),
        sharex=True,
        constrained_layout=True
    )

    x = np.arange(len(class_labels))

    for ax, rsd_vals, (z_top, z_bottom) in zip(axes, rsd_all, depth_slices_cm):
        ax.bar(x, rsd_vals, width=0.8)
        ax.set_ylabel("Root surface density\n(cm$^2$/cm$^3$)")
        ax.set_title(f"Depth slice {z_top}–{z_bottom} cm")
        ax.grid(True, axis="y", alpha=0.3)
        ax.set_ylim(0, ymax * 1.08 if ymax > 0 else 1.0)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(class_labels, rotation=35, ha="right")
    axes[-1].set_xlabel("Root diameter classes (mm)")

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.show()

    # ---- Output
    print("\nRoot surface density by diameter class")
    if time_day is not None:
        print(f"Time filtered up to: {time_day:.1f} d")

    for (z_top, z_bottom), rsd_vals in zip(depth_slices_cm, rsd_all):
        print(f"\nDepth slice {z_top}-{z_bottom} cm:")
        for lbl, val in zip(class_labels, rsd_vals):
            print(f"  {lbl} mm: {val:.6f} cm²/cm³")

plot_root_surface_density_by_diameter(
    plant,
    depth_slices_cm=((0, 7), (7, 73), (73, 80)),
    cyl_radius_cm=cyl_r_cm,
    time_day=48.0,
    diameter_bins_mm=np.arange(0.0, 0.75, 0.05),
    out_path=Path("plots/root_surface_density_by_diameter.png"),
)



