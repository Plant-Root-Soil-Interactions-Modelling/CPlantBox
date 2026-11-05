import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import plantbox as pb

import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import plantbox.visualisation.vtk_plot as vp


def check_data():
    # Load the CSV file
    df = pd.read_csv("DRFW.csv")
    df = df[(df["IRR"] == 'WW') &  (df["depth"] < -100)& (df["depth"] > df["depth"].min())&  (df["DaS"] == 20)]
    irr_values = ['WW']#,'WD']#sorted(df["IRR"].unique())
    var_values = ['D','S','M']
        
       
    df["depth_cm"] = ( - df["depth"] / 10.).astype(int)
    # Get unique values
    irr_values = ['WW']#,'WD']#sorted(df["IRR"].unique())
    var_values = ['D','S','M']#sorted(df["VAR"].unique())

    # Create subplot grid
    fig, axes = plt.subplots(
        len(irr_values),
        len(var_values),
        figsize=(4*len(var_values), 4*len(irr_values)),
        squeeze=False
    )

    # Loop over panels
    for i, irr in enumerate(irr_values):
        for j, var in enumerate(var_values):
            ax = axes[i, j]
            
            subset = df[(df["IRR"] == irr) & (df["VAR"] == var)& (df["depth"] < -100 )]#& (df["DaS"] >= 42)]# ]
            
            if subset.empty:
                ax.axis("off")
                continue
            drfw = subset["DRFW"]/sum(subset["DRFW"]) * 100.
            depth_cm = -subset["depth"]/10
            # Correct scatter plot
            sc = ax.scatter(
                drfw,      # x = RLD/DRFW
                depth_cm,     # y = depth
                c=subset["DaS"],     # color by DaS
                cmap="coolwarm"
            )
            
            ax.set_title(f"IRR={irr}, VAR={var}")
            ax.set_xlabel("DRFW")
            ax.set_ylabel("Depth")
            
            # Invert depth axis (important!)
            ax.invert_yaxis()

    # One shared colorbar
    fig.colorbar(sc, ax=axes, label="DaS")

    plt.tight_layout()
    plt.show()
check_data()
raise Exception

def get_rld(r, k,ks, sim_time):
    #print('r, k, sim_time',r, k, sim_time)
    
    allRS = []
    for numplant in range(1):
        plant = pb.MappedPlant()

        param_path = Path("../../modelparameter/structural/rootsystem")
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
            
        def get_or_create_org_param(plant: pb.Plant, sub_type: int, ot = pb.root) -> pb.RootRandomParameter:
            plist = plant.getOrganRandomParameter(ot)  # LISTE
            for p in plist:
                if hasattr(p, "subType") and int(p.subType) == int(sub_type):
                    return p
            if ot == pb.root:
                p_new = pb.RootRandomParameter(plant)
                p_new.subType = sub_type
                plant.setOrganRandomParameter(p_new)
            if ot == pb.stem:
                p_new = pb.StemRandomParameter(plant)
                p_new.subType = sub_type
                plant.setOrganRandomParameter(p_new)
            return p_new

        s1 = get_or_create_org_param(plant, 1, pb.stem)  # stem
        
        p0 = get_or_create_org_param(plant, 1)  # taproot
        p1 = get_or_create_org_param(plant, 2)  # lateral 1st
        p2 = get_or_create_org_param(plant, 3)  # lateral 2nd
        p4 = get_or_create_org_param(plant, 4)  # seminal roots
        p5 = get_or_create_org_param(plant, 5)  # crown roots

        # --- stem (subType 1)
        
        s1.name = "mainstem"
        s1.a = 0.023612          # [cm] 
        s1.a_s = 0.
        s1.lb = 0.5        # [cm]
        s1.lbs = 0.
        s1.la = 0.5        # [cm]
        s1.las = 0.
        s1.ln = 1       # [cm]
        s1.lns = 0.
        s1.lmax = 5 #83.0           # [cm]
        s1.lmaxs = 0#83.0           # [cm]
        ## Linear Growth
        s1.r = s1.lb/21. # first_crown            # [cm/day]
        s1.ldelay = 0.5
        #p0.rs = r / 10.#2.197             # [cm/day]
        s1.gf = pb.GrowthFunctionType.linear
        s1.successorOT = [[pb.root]]
        s1.successorST = [[5]]
        s1.successorP = [[1.0]]
        s1.successorWhere = [[0]]
        s1.successorNo = [11]
        
        
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
        p0.lmax = 0 #83.0           # [cm]
        p0.lmaxs = 0#83.0           # [cm]
        ## Linear Growth
        p0.r = r #2.197             # [cm/day]
        #p0.rs = r / 10.#2.197             # [cm/day]
        p0.gf = pb.GrowthFunctionType.linear
        ## Sigmoid Growth
        #p0.r = 6.349
        #p0.gf = pb.GrowthFunctionType.sigmoid
        #print(p0)
        #p0.successor = [[2]]
        #p0.successorP = [[1.0]]

        # --- seminal (subType 1)
        p4.name = "seminal"
        p4.a = 0.023612          # [cm] 
        p4.a_s = 0.0132
        p4.lb = 0.342549        # [cm]
        p4.lbs = 0.359216
        p4.la = 3.075370        # [cm]
        p4.las = 3.747081
        p4.ln = 0.351840        # [cm]
        p4.lns = 0.197011
        p4.lmax = k #83.0           # [cm]
        p4.lmaxs = ks#83.0           # [cm]
        ## Linear Growth
        p4.r = r #2.197             # [cm/day]
        #p0.rs = r / 10.#2.197             # [cm/day]
        p4.gf = pb.GrowthFunctionType.linear
        p4.successor = [[2]]
        p4.successorP = [[1.0]]
        ## Sigmoid Growth
        #p0.r = 6.349
        #p0.gf = pb.GrowthFunctionType.sigmoid
        #print(p0)
        #p0.successor = [[2]]
        #p0.successorP = [[1.0]]
        
        # --- crown (subType 1)
        p5.name = "crown"
        p5.a = 0.023612          # [cm] 
        p5.a_s = 0.0132
        p5.lb = 0.342549        # [cm]
        p5.lbs = 0.359216
        p5.la = 3.075370        # [cm]
        p5.las = 3.747081
        p5.ln = 0.351840        # [cm]
        p5.lns = 0.197011
        p5.lmax = k #83.0           # [cm]
        p5.lmaxs = ks#83.0           # [cm]
        ## Linear Growth
        p5.r = r #2.197             # [cm/day]
        #p0.rs = r / 10.#2.197             # [cm/day]
        p5.gf = pb.GrowthFunctionType.linear
        p5.successor = [[2]]
        p5.successorP = [[1.0]]
        ## Sigmoid Growth
        #p0.r = 6.349
        #p0.gf = pb.GrowthFunctionType.sigmoid
        #print(p0)
        #p0.successor = [[2]]
        #p0.successorP = [[1.0]]
        
        
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
        srp.maxB   = 5
        srp.firstB = 0
        srp.delayB = 10
        srp.numB = 5
        srp.delayDefinitionShoot = 1
        plant.setOrganRandomParameter(srp)

        # Experiment Setup
        simtime_days = sim_time #48.0
        cyl_d_cm = 11.0
        cyl_r_cm = cyl_d_cm / 2.0
        cyl_depth_cm = 80.0

        cylinder = pb.SDF_PlantContainer(cyl_r_cm, cyl_r_cm, cyl_depth_cm, False)
        plant.setGeometry(cylinder)

        # Simulate
        plant.initialize()
        plant.simulate(simtime_days)
        
        allRS.append(plant) #need to save the organism or the weak_ptr in SegmentAnalyser will expire
        if numplant == 0:
            allAna = pb.SegmentAnalyser(plant) #need to do that first otherwise get errors
        else:
            allAna.addSegments(plant)  
        #vp.plot_plant(plant, "subType")

    #ana = pb.SegmentAnalyser(base_ana)  
    ana = pb.SegmentAnalyser(allAna)     
    distrib = ana.distribution("length", -10., -cyl_depth_cm, int(cyl_depth_cm - 10), True)  # Länge pro Layer
    # vp.plot_plant(plant, "subType")
    # print("Done. Roots:",k, r, simtime_days, distrib, plant.getNumberOfOrgans())
    rld = np.array(distrib)  # [cm]  / layer_volume
    rld /= sum(rld) # [-]
    rld *= 100. # [%]
    return rld
    
# Load the CSV file
df = pd.read_csv("DRFW.csv")
df = df[(df["IRR"] == 'WW') &  (df["depth"] < -100)&  (df["DaS"] > 42)]
irr_values = ['WW']#,'WD']#sorted(df["IRR"].unique())
var_values = ['D','S','M']
    
   
df["depth_cm"] = ( - df["depth"] / 10.).astype(int)


def target(r, k,ks, dt_temp, VAR = 'S'):
    """ target function for optimization for fit_taproot_r, fit_taproot_rk """
    times = sorted(dt_temp["DaS"].unique())
    sum_t = 0.
    for tt in times:
        subset = dt_temp[dt_temp["DaS"] == tt]
        drfw_array = subset["DRFW"].to_numpy()

        sim_rld = get_rld(r, k,ks, tt)
        #print('lens', len(sim_rld), len(sim_rld))
        for i in range(len(sim_rld)):
            sum_t += (sim_rld[i] - drfw_array[i]/2) ** 2 # half real drfw as we only simulate one plant 

    return np.sqrt(sum_t)


def fit_taproot_rk(df, VAR = 'S'):
    """ fits lenght k """
    subset = df[df["VAR"] == VAR]
    subset_sum = (
        subset.groupby(["depth_cm", "DaS"], as_index=False)  # keep columns
        .agg({"DRFW": "sum"})
    )

    subset_sum["DRFW"] = subset_sum["DRFW"]/sum(subset_sum["DRFW"]) * 100.
    
    f = lambda x0: target(x0[0], x0[1], 0., subset_sum)
    x0 = [2.197  , 100]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    print(res.x[0], res.x[1], f(res.x))
    return res.x[0], res.x[1], f(res.x)
    
def fit_taproot_kks(df, VAR = 'S'):
    """ fits lenght k """
    subset = df[df["VAR"] == VAR]
    subset_sum = (
        subset.groupby(["depth_cm", "DaS"], as_index=False)  # keep columns
        .agg({"DRFW": "sum"})
    )

    subset_sum["DRFW"] = subset_sum["DRFW"]/sum(subset_sum["DRFW"]) * 100.
    
    f = lambda x0: target(2.197, x0[0], x0[1],  subset_sum)
    x0 = [100, 100]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    print(res.x[0], res.x[1], f(res.x))
    return res.x[0], res.x[1], f(res.x)
    
def fit_taproot_k(df, VAR = 'S'):
    """ fits lenght k """
    subset = df[df["VAR"] == VAR]
    subset_sum = (
        subset.groupby(["depth_cm", "DaS"], as_index=False)  # keep columns
        .agg({"DRFW": "sum"})
    )

    subset_sum["DRFW"] = subset_sum["DRFW"]/sum(subset_sum["DRFW"]) * 100.
    
    f = lambda x0: target(2.197, x0[0],0.,  subset_sum)
    x0 = [100]
    res = minimize(f, x0, method = 'Nelder-Mead', tol = 1e-6)  # bounds and constraints are possible, but method dependent
    print(res.x[0], f(res.x))
    return res.x[0], f(res.x)
    
# VAR == S 
# print(fit_taproot_rk(df))
# S: (np.float64(3.185649999999998), np.float64(105.0), np.float64(5.9933491276946445))
# print(fit_taproot_k(df))
# np.float64(104.06189441680905), np.float64(13.980100786982652)

# VAR == D 
# print(fit_taproot_rk(df))
# D: (np.float64(2.2752213593632495), np.float64(100.26222229344657), np.float64(11.67038230567868))
# print(fit_taproot_k(df, VAR = 'D'))
# (np.float64(112.46063172817003), np.float64(21.367354071808577))

# print(fit_taproot_kks(df, VAR = 'S'))
# print(fit_taproot_kks(df, VAR = 'D'))
# (np.float64(100.00092029571533), np.float64(102.49928712844849), np.float64(19.982049025712087))


plt.rcParams.update({
    'font.size': 14,          # Standard Schrift
    'axes.titlesize': 16,     # Titel
    'axes.labelsize': 14,     # Achsenbeschriftung
    'xtick.labelsize': 13,    # x Tick Labels
    'ytick.labelsize': 13,    # y Tick Labels
    'legend.fontsize': 12,    # Legende
    'figure.titlesize': 20    # Figure Titel
})




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

# plot_root_system_2d(plant, stride=1, lw=0.8)

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


# times_days = [5, 10, 15, 20, 25, 30, 35, 40, 45, 48]
# plot_rld_over_time(
    # plant,
    # times_days=times_days,
    # depth_cm=cyl_depth_cm,
    # radius_cm=cyl_r_cm,
    # layers=32,
    # out_path=Path("plots/rld_over_time.png"),
# )

# plot_rld_tap_basal_over_time(
    # plant,
    # times_days=times_days,
    # depth_cm=cyl_depth_cm,
    # radius_cm=cyl_r_cm,
    # layers=32,
    # out_path=Path("plots/rld_tap_basal_over_time.png"),
# )

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



