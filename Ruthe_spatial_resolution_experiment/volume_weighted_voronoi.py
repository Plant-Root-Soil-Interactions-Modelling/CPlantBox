"""
Script to simulate a single plant in a periodic domain and analyze perirhizal radii using 3D Voronoi polygons.

Features:
- Analyzes perirhizal space using 3D Voronoi tessellation with periodic boundary conditions (8 neighbors).
- Computes perirhizal outer radius for each root segment.
- Plots radii histograms with a fitted lognormal distribution for different depth layers.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plantbox as pb
from ruthe_global_variables import (  # noqa
    RUTHE_1994_95,
    RUTHE_1995_96,
    RUTHE_1996_97,
    RutheConfig,
)
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull, QhullError, Voronoi
from scipy.stats import lognorm
from simulation_utils import (
    compute_temperature_scaling,
    create_scale_elongation_grid,
    initialize_soil_temperature,
    load_soil_temperatures,
)

# Set parameters for the different vegetation periods
# 1995: RUTHE_1994_95, 227
# 1996: RUTHE_1995_96, 238
# 1997: RUTHE_1996_97, 224
DATASET = RUTHE_1994_95
SOIL_CORE_DAY = 227
year = 1997

SEED = 56789

CV_RLD_A: float = -42.96712137374045
CV_RLD_B: float = 33.359417519684904

def cv_from_rld(rld: float) -> float:
    return CV_RLD_A * np.log10(rld) + CV_RLD_B

def pdf_perirhizal_mm_analytic(
    r_mm: np.ndarray, mu_ln: float, sigma_ln: float
) -> np.ndarray:
    """Closed-form PDF for radius in mm induced by a lognormal RLD.

    Assumptions:
        - X = RLD is lognormal with ln(X) ~ N(mu_ln, sigma_ln^2)
        - Radius mapping: r_mm = 10 * sqrt(1 / (pi * X))

    Then r_mm is also lognormal with
        ln(r_mm) ~ N(mu_r, sigma_r^2)
        mu_r    = ln(10) - 0.5*ln(pi) - 0.5*mu_ln
        sigma_r = 0.5*sigma_ln
    """
    mu_r = np.log(10.0) - 0.5 * np.log(np.pi) - 0.5 * mu_ln
    sigma_r = 0.5 * sigma_ln
    return lognorm.pdf(r_mm, s=sigma_r, scale=np.exp(mu_r))

class PeriodicDomainAnalyser:
    def __init__(
        self,
        root_system: pb.RootSystem,
        domain_size: tuple[float, float],
        depth_layers: list[tuple[float, float]] = None,
    ):
        """
        Args:
            root_system: The simulated root system.
            domain_size: (Lx, Ly) dimensions of the periodic domain.
            depth_layers: List of (z_top, z_bottom) tuples for analysis.
        """
        self.rs = root_system
        self.Lx, self.Ly = domain_size

        if depth_layers is None:
            self.depth_layers = [
                (0, -30),
                (-30, -60),
                (-60, -90),
                (-90, -120),
            ]
        else:
            self.depth_layers = depth_layers

        self.segments_df = pd.DataFrame()
        self.layer_rld: dict[int, float] = {}

    def run_analysis(self):
        self.extract_nodes_and_segments()
        self.compute_periodic_bounded_voronoi()
        self.distribute_volumes_to_segments()
        self.compute_perirhizal_radius()
        self.assign_depth_layers()
        self.compute_root_length_densities()

        print("Analysis complete.")

    def extract_nodes_and_segments(self):
        # Extract nodes, segments and their radii from the root system
        sa = pb.SegmentAnalyser(self.rs)

        nodes = sa.nodes
        segments = sa.segments
        radii = sa.getParameter("radius")

        # Build clean list of nodes with their index and 3D coordinates
        node_data = []
        for i, n in enumerate(nodes):
            node_data.append({"node_id": i, "x": n.x, "y": n.y, "z": n.z})

        # Build clean list of segments with their index, node indices, length, radius, and depth
        segment_data = []
        for i, seg in enumerate(segments):
            n1 = seg.x if hasattr(seg, "x") else seg[0]  # index of first node
            n2 = seg.y if hasattr(seg, "y") else seg[1]  # index of second node

            p1 = np.array(
                [nodes[n1].x, nodes[n1].y, nodes[n1].z]
            )  # coordinates of first node
            p2 = np.array(
                [nodes[n2].x, nodes[n2].y, nodes[n2].z]
            )  # coordinates of second node

            segment_data.append(
                {
                    "seg_id": i,
                    "n1": n1,
                    "n2": n2,
                    "length": np.linalg.norm(
                        p2 - p1
                    ),  # Euclidean length of the segment
                    "root_radius": radii[i],
                    "z": 0.5 * (p1[2] + p2[2]),  # midpoint depth
                }
            )

        self.nodes_df = pd.DataFrame(node_data)
        self.segments_df = pd.DataFrame(segment_data)

        print(f"Extracted {len(self.segments_df)} segments.")
        if not self.segments_df.empty:
            print(
                f"Depth range: {self.segments_df.z.max():.2f} to {self.segments_df.z.min():.2f} cm"
            )
        else:
            print("Warning: No segments extracted!")

    def compute_periodic_bounded_voronoi(self, z_top=0.0, z_bottom=-120.0):
        # Override bottom boundary based on deepest node
        deepest_z = self.nodes_df.z.min()
        z_bottom = (
            deepest_z - 0.5  # arbitrary value that might need to be changed!!!!!
        )
        print(
            f"Adjusting Voronoi bottom boundary to {z_bottom:.2f} cm (deepest node: {deepest_z:.2f} cm)"
        )

        # Original coordinates
        original_points = self.nodes_df[["x", "y", "z"]].values
        n = len(original_points)

        # Wrap coordinates into fundamental domain [-Lx/2, Lx/2] x [-Ly/2, Ly/2]
        # Periodic wrap ensures all point of the original plant are within the domain
        # [-Lx/2, Lx/2] x [-Ly/2, Ly/2] ensures centering around the seed point (0,0)
        points_wrapped = original_points.copy()
        points_wrapped[:, 0] = (
            points_wrapped[:, 0] + self.Lx / 2
        ) % self.Lx - self.Lx / 2
        points_wrapped[:, 1] = (
            points_wrapped[:, 1] + self.Ly / 2
        ) % self.Ly - self.Ly / 2

        # Central points (wrapped) - FIRST so indices 0..n-1 match the original nodes
        all_points = [points_wrapped]

        # Add 8 periodic neighbors in XY plane using the WRAPPED points
        # Voronoi polygons arent periodic; this mimics periodicity by adding neighbors
        # which act as stand-ins for neighboring roots on the other side of the periodic boundary.
        shifts_x = [-self.Lx, 0, self.Lx]
        shifts_y = [-self.Ly, 0, self.Ly]

        for dx in shifts_x:
            for dy in shifts_y:
                if dx == 0 and dy == 0:
                    continue  # Center (the fundamental domain) is skipped
                shifted = points_wrapped + np.array([dx, dy, 0.0])
                all_points.append(shifted)

        # Z boundaries defined through mirroring
        # Ensures flat boundaries at z_top and z_bottom
        top = points_wrapped.copy()
        top[:, 2] = 2 * z_top - top[:, 2]

        bottom = points_wrapped.copy()
        bottom[:, 2] = 2 * z_bottom - bottom[:, 2]

        all_points.extend([top, bottom])
        all_points = np.vstack(all_points)

        # Compute Voronoi
        try:
            vor = Voronoi(all_points)
        except QhullError as e:
            print(f"Voronoi computation failed: {e}")
            self.nodes_df["voronoi_volume"] = np.nan
            return

        # Prepare array for volumes of the ORIGINAL points only (= first n points)
        volumes = np.zeros(n)

        # Loop over original points only
        for i in range(n):
            # Get the Voronoi region for point i
            region_idx = vor.point_region[i]
            region = vor.regions[region_idx]

            if -1 in region or len(region) == 0:  # detects invalid & unbounded regions
                volumes[i] = np.nan
                print("Warning: Unbounded Voronoi region detected for node", i)
                continue

            try:
                hull = ConvexHull(
                    vor.vertices[region]
                )  # assembles convex hull of the region
                volumes[i] = hull.volume  # computes volume of the convex hull
            except QhullError:
                volumes[i] = np.nan

        self.nodes_df["voronoi_volume"] = volumes

        # Diagnostic check - total volume vs theoretical
        total_vol = np.nansum(volumes)
        z_range = self.nodes_df.z.max() - z_bottom
        theoretical_vol = self.Lx * self.Ly * z_range
        print(f"Total Voronoi Volume: {total_vol:.2f} cm3")
        print(
            f"Theoretical Domain Volume ({self.Lx}x{self.Ly}x{z_range:.2f}): {theoretical_vol:.2f} cm3"
        )
        if total_vol > 1.2 * theoretical_vol:
            print(
                "WARNING: Total volume exceeds theoretical by >20%! Bounding might be insufficient."
            )

    def distribute_volumes_to_segments(
        self,
    ):
        node_vol = self.nodes_df["voronoi_volume"].to_dict()

        seg_volumes = []
        for _, seg in self.segments_df.iterrows():
            # Assign each segment the full Voronoi volume of its apical node.
            v = node_vol.get(seg.n2, np.nan)
            seg_volumes.append(v)

        self.segments_df["voronoi_volume"] = seg_volumes

    def compute_perirhizal_radius(self):
        df = self.segments_df

        # Filter very short segments to avoid numerical artifacts
        # here: length < 0.01 cm (100 microns)
        valid = (df.voronoi_volume > 0) & (df.length > 0.01)

        # Compute perirhizal outer radius from volume and length
        df.loc[valid, "perirhizal_radius"] = np.sqrt(
            df.loc[valid, "voronoi_volume"]
            / (np.pi * df.loc[valid, "length"])  # R = sqrt(V / (pi * L)
        )

        # Report discarded
        discarded_count = len(df) - valid.sum()
        print(
            f"Computed radii for {valid.sum()} segments. Discarded {discarded_count} (short or invalid vol)."
        )

        # Diagnostic check: max radius and suspiciously large ones
        if valid.sum() > 0:
            max_r = df.loc[valid, "perirhizal_radius"].max()
            print(f"Max computed radius: {max_r:.2f} cm")

            # Print problematic ones if any
            suspicious = df[(df.perirhizal_radius > 10.0)]
            if not suspicious.empty:
                print(f"Found {len(suspicious)} segments with Radius > 10cm. Sample:")
                print(
                    suspicious[["length", "voronoi_volume", "perirhizal_radius"]].head()
                )

    def assign_depth_layers(self):
        self.segments_df["layer_id"] = -1  # initilize everything as invalid

        # Loop over depth layers and assign layer IDs
        for i, (zt, zb) in enumerate(self.depth_layers):
            upper, lower = max(zt, zb), min(zt, zb)
            mask = (self.segments_df.z <= upper) & (self.segments_df.z > lower)
            self.segments_df.loc[mask, "layer_id"] = (
                i  # assign layer ID based on segment midpoint depth
            )

    def compute_root_length_densities(self):
        """
        Computes and prints the Root Length Density (RLD) for each depth layer.
        RLD = (Total Root Length in Layer) / (Soil Volume of Layer)
        """
        print("\nRoot Length Densities (RLD):")
        self.layer_rld = {}
        for i, (zt, zb) in enumerate(self.depth_layers):
            # Calculate layer volume
            depth_range = abs(zt - zb)
            layer_volume = self.Lx * self.Ly * depth_range

            # Sum segment lengths for this layer
            layer_segments = self.segments_df[self.segments_df["layer_id"] == i]
            total_length = layer_segments["length"].sum()

            if layer_volume > 0:
                rld = total_length / layer_volume
            else:
                rld = 0.0

            self.layer_rld[i] = rld

            print(f"Layer {i} ({zt} to {zb} cm): RLD = {rld:.6f} cm/cm3")

    def plot_layer_histograms_with_lognormal_fit(
        self,
        output_prefix: str
    ):
        """Histogram of perirhizal radii with fitted lognormal PDF.

        The histogram is a *volume-weighted* PDF (units 1/mm):
            height(bin) = (sum(voronoi_volume in bin) / sum(voronoi_volume in layer)) / bin_width

        Lognormal and RLD-derived curves are plotted on the same axis (also 1/mm).
        """

        data = self.segments_df.dropna(subset=["perirhizal_radius", "layer_id"])
        data = data[data["layer_id"] != -1]

        sorted_layers = sorted(
            range(len(self.depth_layers)),
            key=lambda k: (self.depth_layers[k][0] + self.depth_layers[k][1]) / 2,
            reverse=True,
        )

        simulation_data = {}

        for idx in sorted_layers:
            z_top, z_bot = self.depth_layers[idx]

            layer_data = data[data["layer_id"] == idx]
            radii = (layer_data["perirhizal_radius"].values * 10.0)  # mm
            volumes = layer_data["voronoi_volume"].values  # cm3 (segment-associated)

            valid = (
                np.isfinite(radii)
                & np.isfinite(volumes)
                & (radii > 0)
                & (volumes > 0)
            )
            radii = radii[valid]
            volumes = volumes[valid]

            if len(radii) < 20:
                print(f"Skipping layer {idx} (not enough data)")
                continue

            total_volume = volumes.sum()
            if total_volume <= 0:
                print(f"Skipping layer {idx} (non-positive total volume)")
                continue

            # --- compute histogram and (optional) fits first, then decide whether we need a broken y-axis ---
            vol_in_bin, bins = np.histogram(radii, bins=50, weights=volumes)
            bin_centers = 0.5 * (bins[:-1] + bins[1:])
            dx = np.diff(bins)

            # Volume-weighted probability density (1/mm)
            density = np.divide(
                vol_in_bin,
                total_volume * dx,
                out=np.zeros_like(vol_in_bin, dtype=float),
                where=dx > 0,
            )

            # Lognormal fit (to volume-weighted density)
            lognorm_fit = None  # (x_fit, pdf_fit, median)

            try:
                popt, _ = curve_fit(
                    lambda x, s, scale: lognorm.pdf(x, s=s, scale=scale),
                    bin_centers,
                    density,
                    p0=[1.0, float(np.mean(radii))],
                    maxfev=10000,
                )

                s_fit, scale_fit = popt
                x_fit = np.linspace(float(radii.min()), float(radii.max()), 300)
                pdf_fit = lognorm.pdf(x_fit, s=s_fit, scale=scale_fit)
                mode = scale_fit / np.exp(s_fit**2)
                lognorm_fit = (x_fit, pdf_fit, float(mode))

                simulation_data[idx] = {"s": float(s_fit), "scale": float(scale_fit)}
            except Exception as e:
                print(f"Fit failed in layer {idx}: {e}")

            # RLD-derived PDF (density)
            rld_curve = None
            layer_rld = self.layer_rld.get(idx, np.nan)

            if np.isfinite(layer_rld) and layer_rld > 0:
                cv = max(cv_from_rld(layer_rld), 0.0)  # ensure non-negative CV
                sigma = cv * layer_rld / 100.0
                
                if sigma >= 0:
                    sigma_ln = np.sqrt(np.log(1 + (sigma / layer_rld) ** 2))
                    mu_ln = np.log(layer_rld) - 0.5 * sigma_ln**2

                    x_rld = np.logspace(-3, 1, 500)
                    r_mm = np.sqrt(1/(np.pi * x_rld)) * 10.0  # mm
                    pdf_r_vals = pdf_perirhizal_mm_analytic(r_mm, mu_ln, sigma_ln)

                    order = np.argsort(r_mm)
                    r_mm_sorted = r_mm[order]
                    pdf_sorted = pdf_r_vals[order]

                    mode_idx = int(np.nanargmax(pdf_sorted))
                    r_mode = float(r_mm_sorted[mode_idx])
                    peak_height = float(pdf_sorted[mode_idx])

                    if np.isfinite(peak_height) and peak_height > 0:
                        rld_curve = (
                            r_mm_sorted,
                            pdf_sorted,
                            r_mode,peak_height
                        )
                else:
                    print(
                        f"Skipped RLD-derived PDF in layer {idx}: "
                        f"non-positive sigma ({sigma:.4g}%)."
                    )

            else:
                print(
                    f"Skipped RLD-derived PDF in layer {idx}: "
                    f"invalid RLD ({layer_rld})."
                )

            # If the RLD-derived curve peak (at its mode) is above 0.5, use a broken y-axis so we can
            # still keep the main plot limited to (0, 0.5) while showing the peak.
            use_broken_y = (rld_curve is not None and rld_curve[3] > 0.5)

            # Deterministic layout (avoid tight_layout differences between 1-axes and broken-axes plots)
            # Use fixed margins so layout stays identical across plot types.
            fig_left = 0.12
            fig_right = 0.78
            fig_bottom = 0.14
            fig_top = 0.86

            if use_broken_y:
                fig, (ax_top, ax_bot) = plt.subplots(
                    2,
                    1,
                    sharex=True,
                    figsize=(6, 4),
                    gridspec_kw={"height_ratios": [1, 3], "hspace": 0.05},
                )

                # Ensure the top axes (and its legend) draws above the bottom axes.
                # Without this, the bottom axes patch can cover parts of artists in the top axes
                # in the small overlap region, which can make the legend frame look "cut".
                ax_top.set_zorder(2)
                ax_bot.set_zorder(1)
                ax_top.patch.set_visible(False)

                # Break markers
                ax_top.spines["bottom"].set_visible(False)
                ax_bot.spines["top"].set_visible(False)
                ax_top.tick_params(labeltop=False)
                ax_top.xaxis.set_ticks_position("none")

                d = 0.015
                kwargs = dict(color="k", clip_on=False, linewidth=1)
                ax_top.plot((-d, +d), (-d, +d), transform=ax_top.transAxes, **kwargs)
                ax_top.plot((1 - d, 1 + d), (-d, +d), transform=ax_top.transAxes, **kwargs)
                ax_bot.plot((-d, +d), (1 - d, 1 + d), transform=ax_bot.transAxes, **kwargs)
                ax_bot.plot((1 - d, 1 + d), (1 - d, 1 + d), transform=ax_bot.transAxes, **kwargs)

                axes_main = ax_bot
                axes_peak = ax_top
            else:
                fig, ax = plt.subplots(figsize=(6, 4))
                axes_main = ax
                axes_peak = None

            # --- plot histogram on main axis ---
            axes_main.bar(
                bin_centers,
                density,
                width=dx,
                alpha=0.5,
                label="Histogram (volume-weighted PDF)",
            )

            # --- plot lognormal fit and empirical weighted median on main axis ---
            if lognorm_fit is not None:
                x_fit, pdf_fit, _fitted_median = lognorm_fit
                axes_main.plot(x_fit, pdf_fit, color="orange", lw=2, label="Lognormal fit")
                axes_main.axvline(
                    mode,
                    color="orange",
                    linestyle=":",
                    lw=2,
                    label=f"Mode {mode:.2f} mm",
                )
            # --- plot RLD-derived curve (main axis always clipped to 0..0.5) ---
            if rld_curve is not None:

                axes_main.plot(
                    r_mm_sorted,
                    pdf_sorted,
                    color="green",
                    lw=2,
                    label="RLD-derived PDF",
                )

                axes_main.axvline(
                    r_mode,
                    color="green",
                    linestyle="--",
                    lw=2,
                    label=f"RLD median {r_mode:.2f} mm",
                )

                # If peak is above 0.5, show it in the top axis with a truncated/broken y-range.
                if axes_peak is not None:
                    y_top_max = peak_height * 1.05
                    # Show only the upper part around the peak; this defines the "skipped" region.
                    y_top_min = max(0.55, peak_height * 0.8)
                    if y_top_min >= y_top_max:
                        y_top_min = 0.55

                    axes_peak.plot(r_mm_sorted, pdf_sorted, color="green", lw=2)

                    axes_peak.axvline(
                        r_mode,
                        color="green",
                        linestyle="--",
                        lw=2,
                    )
                    
                    axes_peak.set_ylim(y_top_min, y_top_max)

                    # Reduce clutter on the peak axis
                    axes_peak.set_ylabel("")
                    axes_peak.tick_params(labelbottom=False)

                    print(
                        f"Layer {idx}: RLD-derived peak at mode is {peak_height:.3g} 1/mm (> 0.5); using broken y-axis."
                    )

            # --- formatting ---
            title = f"Depth {max(z_top, z_bot)} to {min(z_top, z_bot)} cm"
            fig.suptitle(title, y=0.95)

            xlabel = "Perirhizal radius (mm)"
            ylabel = "Volume-weighted probability density (1/mm)"

            axes_main.set_xlabel(xlabel)
            if use_broken_y:
                # Center the y-label across both panels for broken-axis plots
                axes_main.set_ylabel("")
                fig.text(0.04, 0.5, ylabel, ha="left", va="center", rotation="vertical")
            else:
                axes_main.set_ylabel(ylabel)

            axes_main.set_xlim(0, 30)
            axes_main.set_ylim(0, 0.5)

            # Fixed margins for consistent title/label spacing across both plot types
            if use_broken_y:
                fig.subplots_adjust(
                    left=fig_left,
                    right=fig_right,
                    bottom=fig_bottom,
                    top=fig_top,
                    hspace=0.05,
                )
            else:
                fig.subplots_adjust(
                    left=fig_left,
                    right=fig_right,
                    bottom=fig_bottom,
                    top=fig_top,
                )

            # Legend: place inside the plot, top-right.
            handles, labels = axes_main.get_legend_handles_labels()
            if handles:
                legend_ax = axes_peak if use_broken_y and (axes_peak is not None) else axes_main
                leg = legend_ax.legend(
                    handles,
                    labels,
                    loc="upper right",
                    fontsize="small",
                    frameon=True,
                )
                leg.set_zorder(10)

            # --- save per layer ---
            path = f"{output_prefix}_layer_{idx}_lognormal.png"
            plt.savefig(path)
            plt.close()

            print(f"Saved {path}")

        print("Simulation parameters:", simulation_data)


def run_daily_growth_simulation(
    root_system: pb.RootSystem,
    df: pd.DataFrame,
    scale_elongation: pb.EquidistantGrid1D,
    dataset: RutheConfig,
    analysis_days: list[int] = None,
    callback=None,
) -> None:
    """
    Run daily growth simulation for the root system.

    Args:
        root_system (pb.RootSystem): The root system to simulate.
        df (pd.DataFrame): DataFrame containing soil temperature data.
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        dataset: Module containing global variables for a certain experiment.
        analysis_days: List of days to run analysis on.
        callback: Function to call on analysis days. Signature: callback(root_system, day)
    """

    sim_time = dataset.SIMULATION_TIME
    dt = dataset.TIME_STEP

    for time_step in range(round(sim_time / dt)):
        # Check boundary
        if time_step >= len(df["10"]):
            break

        # Update soil temperature and scaling
        soil_temp = initialize_soil_temperature(df, time_step)
        scales = compute_temperature_scaling(soil_temp, dataset)
        scale_elongation.data = scales

        root_system.simulate(dt, False)

        # Callback analysis; allows for perirhizal radii analysis on specified days not just at the end
        current_day = time_step + 1
        if analysis_days and callback and current_day in analysis_days:
            callback(root_system, current_day)


def initilize_sim_root_system(
    scale_elongation: pb.EquidistantGrid1D,
    dataset: RutheConfig,
    seed: int,
) -> pb.RootSystem:
    """
    Initialize a root system with specified parameters.

    Args:
        scale_elongation (pb.EquidistantGrid1D): Scale elongation grid.
        dataset (RutheConfig): Configuration for the experiment.
        seed (int): Random seed for reproducibility.

    Returns:
        pb.RootSystem: Initialized root system.
    """

    root_system = pb.RootSystem()
    root_system.readParameters(str(dataset.PATH_TO_PLANT_PARAMETERS))

    # Set seed for reproducibility
    root_system.setSeed(seed)

    # Assign elongation scaling
    for root_type in root_system.getRootRandomParameter():
        root_type.f_se = scale_elongation

    # Define unbounded soil space
    # PeriodicDomainAnalyser will handle periodicity and actual domain size
    soil_space = pb.SDF_PlantContainer(
        dataset.SOIL_SPACE_PARAMS[0],
        dataset.SOIL_SPACE_PARAMS[1],
        dataset.SOIL_SPACE_PARAMS[2],
        dataset.SOIL_SPACE_PARAMS[3],
    )

    root_system.setGeometry(soil_space)
    root_system.initializeDB(4, 5, False)

    return root_system


def run_single_plant_simulation():
    dataset = DATASET
    output_dir = dataset.OUTPUT_DIRECTORY

    # Load and initialize soil temperature data
    soil_temp_df = load_soil_temperatures(
        Path(dataset.PATH_TO_SOIL_TEMPERATURE_DATA),
        dataset.START_DATE,
        dataset.END_DATE,
    )
    soil_temp = initialize_soil_temperature(soil_temp_df, 0)

    # Compute temperature scaling
    temp_scales = compute_temperature_scaling(soil_temp, dataset)
    scale_elongation = create_scale_elongation_grid(temp_scales)

    # Initialize root system
    rs = initilize_sim_root_system(scale_elongation, dataset, seed=SEED)
    # Forcing the seed position to be at (0,0) in XY plane
    rs.getRootSystemParameter().seedPos = pb.Vector3d(0, 0, dataset.SOWING_DEPTH)
    # Define periodic domain size; will be centered around (0,0) by PeriodicDomainAnalyser
    domain_size = (dataset.INTER_PLANT_DISTANCE, dataset.INTER_ROW_DISTANCE)

    # Analysis callback, allows for perirhizal radii analysis on fixed set of days
    def analysis_callback(root_sys, day):
        """
        Analysis callback to compute and plot perirhizal radii on specified days.
        """
        print(f"\n--- Analyzing Day {day} ---")

        # Create a new analyser instance for the current state
        analyser = PeriodicDomainAnalyser(root_sys, domain_size)
        analyser.run_analysis()

        base_name = f"{output_dir}/simulated_perirhizal_radii/volume_weighted_voronoi_with_confidence_{year}"
        analyser.plot_layer_histograms_with_lognormal_fit(base_name)

    analysis_days = [SOIL_CORE_DAY]  # days to analyze (here: day of soil core sampling)
    analysis_days = sorted(
        list(set([d for d in analysis_days if d <= dataset.SIMULATION_TIME]))
    )

    print("Running growth simulation...")
    run_daily_growth_simulation(
        rs,
        soil_temp_df,
        scale_elongation,
        dataset,
        analysis_days=analysis_days,
        callback=analysis_callback,
    )
    print(f"Simulation done. Nodes: {len(rs.getNodes())}")


if __name__ == "__main__":
    run_single_plant_simulation()
