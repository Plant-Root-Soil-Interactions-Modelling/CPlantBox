"""
Random 3D hyphal / hyphal-tip network inside a unit cube.

Geometry model
--------------
- Hyphae are straight line segments, hyphal tips are the apical (growing) end
  of a segment.
- Hyphal tips are placed uniformly at random inside the cube to realize a
  target hyphal tip density 'n' (tips per unit volume).
- Each tip grows a straight hypha in a random (isotropic) direction; the
  segment extends from the tip to the cube boundary, i.e. the maximal
  length available in that direction.
- If the length contributed by these tip-anchored hyphae is not enough to
  reach the target hyphal length density 'rho_h' (length per unit volume),
  extra hyphae are added as random straight chords that cross the cube
  entirely (both endpoints on the cube surface, no hyphal tip inside the
  domain) until the target length is reached.
"""

import numpy as np
import matplotlib.pyplot as plt

def _random_directions(n, rng):
    """n random directions, uniformly distributed on the unit sphere."""
    v = rng.normal(size=(n, 3))
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / norms

def _random_directions_2d(n, rng):
    theta = rng.uniform(0.0, 2.0 * np.pi, size=n)
    return np.column_stack([
        np.cos(theta),
        np.sin(theta),
    ])

def _ray_square_exit(points, directions,box):
    t = np.full(directions.shape, np.inf)
    pos, neg = directions > 0, directions < 0
    t[pos] = (box - points[pos]) / directions[pos]
    t[neg] = (0.0 - points[neg]) / directions[neg]
    return t.min(axis=1)


def _ray_box_exit(points, directions, box):
    """distance from each point (inside [0, box]^3) to the box boundary along directions"""
    t = np.full(directions.shape, np.inf)
    pos, neg = directions > 0, directions < 0
    t[pos] = (box - points[pos]) / directions[pos]
    t[neg] = (0.0 - points[neg]) / directions[neg]
    return t.min(axis=1)

def create_hyphal_network_2d(n, rho_h, box=1.0, seed=None):
    rng = np.random.default_rng(seed)
    area = box**2
    n_tips = int(round(n * area))
    tips = rng.uniform(0.0,box,size=(n_tips, 2),)

    segments, is_tip_segment = [], []
    if n_tips > 0:
        directions = _random_directions_2d(n_tips,rng,)
        lengths = _ray_square_exit(tips,directions,box,)
        ends = tips + directions * lengths[:, np.newaxis]
        for p0, p1 in zip(tips, ends):
            segments.append((p0.copy(), p1.copy()))
            is_tip_segment.append(True)
    tip_hyphal_length = float(sum(np.linalg.norm(p1 - p0) for p0, p1 in segments))
    total_length = tip_hyphal_length
    target_length = rho_h * area

    while total_length < target_length:
        p = rng.uniform(0.0,box,size=2,)
        d = _random_directions_2d(1,rng,)[0]
        t_fwd = _ray_square_exit(p[np.newaxis, :],d[np.newaxis, :],box,)[0]
        t_bwd = _ray_square_exit(p[np.newaxis, :],-d[np.newaxis, :],box,)[0]
        p0 = p - d * t_bwd
        p1 = p + d * t_fwd
        segments.append((p0,p1,))
        is_tip_segment.append(False)
        total_length += t_fwd + t_bwd
    # ------------------------------------------------------------------
    # 4. Realized densities
    # ------------------------------------------------------------------

    n_realized = n_tips / area
    rho_h_realized = total_length / area

    extra_hyphal_length = total_length - tip_hyphal_length

    return {
        "tips": tips,
        "segments": segments,
        "is_tip_segment": np.array(is_tip_segment, dtype=bool),
        "tip_segment_index": np.arange(n_tips),  # segments[:n_tips] are the tip-anchored ones, in tip order
        "total_length": total_length,
        "volume": area,
        "tip_hyphal_length": tip_hyphal_length,
        "extra_hyphal_length": extra_hyphal_length,
        "n_realized": n_realized,
        "rho_h_realized": rho_h_realized,
        "n_target": n,
        "rho_h_target": rho_h,
        "rho_h_error": rho_h_realized - rho_h,
        "rho_h_relative_error": (rho_h_realized - rho_h) / rho_h if rho_h > 0 else np.nan,
    }


def create_hyphal_network(n, rho_h, box=1.0, seed=None):
    """
    Creates a random 3D hyphal network inside a cube of side length 'box' that
    realizes a hyphal tip density n and a hyphal length density rho_h.

    Parameters:
    n: hyphal tip density [number of apical tips per unit volume]
    rho_h: hyphal length density [total hyphal length per unit volume]
    box: side length of the cubic domain (default 1, i.e. unit cube)
    seed: optional random seed for reproducibility

    Returns a dict with:
    'tips': (n, 3) array of hyphal tip coordinates
    'segments': list of (p0, p1) endpoint pairs, one per hyphal segment
    'is_tip_segment': bool array, True if the segment starts at a hyphal tip
    'tip_segment_index': index of the segment growing from tip i, i.e. segments[tip_segment_index[i]]
    'total_length': realized total hyphal length
    'volume': box**3
    """
    rng = np.random.default_rng(seed)
    volume = box**3

    # 1) place hyphal tips uniformly to realize n
    n_tips = int(round(n * volume))
    tips = rng.uniform(0.0, box, size=(n_tips, 3))

    segments, is_tip_segment = [], []

    # 2) grow one straight hypha per tip, random direction, up to the cube boundary
    if n_tips > 0:
        directions = _random_directions(n_tips, rng)
        lengths = _ray_box_exit(tips, directions, box)
        ends = tips + directions * lengths[:, np.newaxis]
        for p0, p1 in zip(tips, ends):
            segments.append((p0.copy(), p1.copy()))
            is_tip_segment.append(True)
    tip_hyphal_length = float(sum(np.linalg.norm(p1 - p0) for p0, p1 in segments))
    total_length = tip_hyphal_length
    target_length = rho_h * volume

    # 3) if more length is needed, add random hyphae crossing the cube (no tip inside)
    while total_length < target_length:
        p = rng.uniform(0.0, box, size=3)
        d = _random_directions(1, rng)[0]
        t_fwd = _ray_box_exit(p[np.newaxis, :], d[np.newaxis, :], box)[0]
        t_bwd = _ray_box_exit(p[np.newaxis, :], -d[np.newaxis, :], box)[0]
        p0, p1 = p - d * t_bwd, p + d * t_fwd
        segments.append((p0, p1))
        is_tip_segment.append(False)
        total_length += t_fwd + t_bwd
    extra_hyphal_length = total_length - tip_hyphal_length
    n_realized = n_tips / volume
    rho_h_realized = total_length / volume

    return {
        "tips": tips,
        "segments": segments,
        "is_tip_segment": np.array(is_tip_segment, dtype=bool),
        "tip_segment_index": np.arange(n_tips),  # segments[:n_tips] are the tip-anchored ones, in tip order
        "total_length": total_length,
        "volume": volume,
        "tip_hyphal_length": tip_hyphal_length,
        "extra_hyphal_length": extra_hyphal_length,
        "n_realized": n_realized,
        "rho_h_realized": rho_h_realized,
        "n_target": n,
        "rho_h_target": rho_h,
        "rho_h_error": rho_h_realized - rho_h,
        "rho_h_relative_error": (rho_h_realized - rho_h) / rho_h if rho_h > 0 else np.nan,
    }


def _point_segment_distances(points, seg_p0, seg_p1):
    """distance of each point (n, 3) to each segment (seg_p0, seg_p1, both (m, 3)) -> (n, m) matrix"""
    d = seg_p1 - seg_p0
    seg_len2 = np.einsum("ij,ij->i", d, d)
    seg_len2 = np.where(seg_len2 > 0.0, seg_len2, 1.0)  # guard against zero-length segments
    w = points[:, np.newaxis, :] - seg_p0[np.newaxis, :, :]
    t = np.clip(np.einsum("nmj,mj->nm", w, d) / seg_len2[np.newaxis, :], 0.0, 1.0)
    closest = seg_p0[np.newaxis, :, :] + t[..., np.newaxis] * d[np.newaxis, :, :]
    diff = points[:, np.newaxis, :] - closest
    return np.sqrt(np.einsum("nmj,nmj->nm", diff, diff))


def compute_tip_neighbor_distances(network):
    """
    Precomputes, for every hyphal tip, the minimal distance to any other hypha
    (excluding the one segment that grows from that tip itself). The result is
    cached in network['tip_min_distances'] and, sorted, in
    network['tip_min_distances_sorted'], so that neighbor counts for many
    radii R can later be obtained in O(log n) each via count_tips_near_hyphae,
    instead of recomputing the full distance matrix every time.
    """
    tips = network["tips"]
    segments = network["segments"]
    n_tips = len(tips)

    if n_tips == 0 or len(segments) == 0:
        network["tip_min_distances"] = np.full(n_tips, np.inf)
    else:
        seg_p0 = np.array([s[0] for s in segments])
        seg_p1 = np.array([s[1] for s in segments])
        dist = _point_segment_distances(tips, seg_p0, seg_p1)  # (n_tips, n_segments)
        dist[np.arange(n_tips), network["tip_segment_index"]] = np.inf  # exclude own segment
        network["tip_min_distances"] = dist.min(axis=1)

    network["tip_min_distances_sorted"] = np.sort(network["tip_min_distances"])
    return network["tip_min_distances"]


def count_tips_near_hyphae(network,R,hyphal_radius=0.0):
    """
    Count tips whose distance to another hypha surface
    is <= R.
    The geometric distance is calculated to the hyphal
    centerline. Therefore: distance_to_centerline <= R + r
    corresponds to distance_to_hyphal_surface <= R.
    Here:
        R = tip-to-hyphal-surface distance
        r = hyphal radius
    """

    if "tip_min_distances_sorted" not in network:
        compute_tip_neighbor_distances(network)
    R = np.asarray(R, dtype=float)

    # Convert surface distance to centerline distance.
    R_effective = R + hyphal_radius

    counts = np.searchsorted(network["tip_min_distances_sorted"], R_effective, side="right")
    if R.ndim == 0:
        return int(counts)

    return counts

def measure_encounters(network, R,hyphal_radius=0.0):
    n_tips = len(network["tips"])
    if n_tips == 0:
        return 0,0,np.nan
    n_enc = count_tips_near_hyphae(network,R,hyphal_radius)
    f_enc = n_enc / n_tips
    return n_enc, n_tips, f_enc


def compute_theoretical_C(R, hyphal_radius=0.0, mean_segment_length=None):
    """
    Simple geometric estimate for C(R).
    Leading order:
        C(R) = pi (R + r)^2

    With finite segment-length correction:
        C(R) = pi (R+r)^2 [1 + 2(R+r)/(3*l)]
    """
    R_eff = np.asarray(R) + hyphal_radius
    if mean_segment_length is None:
        return np.pi * R_eff**2
    return (np.pi* R_eff**2* (1.0 + 2.0 * R_eff/ (3.0 * mean_segment_length)))

def run_proxy_study(n_values,rho_h_values,R_values,n_reps=50,box=1.0,hyphal_radius=0.0,base_seed=12345):
    """
    Run the complete parameter study.
    Both n and rho_h are varied.
    For every combination
        (n, rho_h, R)
    multiple random network realizations are generated.
    Returns all individual measurements as a list of dictionaries.
    """
    results = []
    seed_counter = 0
    for n_target in n_values:
        for rho_h_target in rho_h_values:
            for rep in range(n_reps):
                seed = base_seed + seed_counter
                seed_counter += 1
                net = create_hyphal_network_2d(n=n_target,rho_h=rho_h_target,box=box,seed=seed)
                
                compute_tip_neighbor_distances(net)

                for R in R_values:
                    n_enc, n_tips, f_enc = measure_encounters(net, R, hyphal_radius)

                    results.append({
                        "n_target": n_target,
                        "rho_h_target": rho_h_target,
                        "n": net["n_realized"],
                        "rho_h": net["rho_h_realized"],
                        "rho_h_error": net["rho_h_error"],
                        "rho_h_relative_error": net["rho_h_relative_error"],
                        "rho_h_tip": net["tip_hyphal_length"] / net["volume"],
                        "rho_h_extra": net["extra_hyphal_length"] / net["volume"],
                        "n_enc": n_enc,
                        "n_tips": n_tips,
                        "R": R,
                        "rep": rep,
                        "f_enc": f_enc,
                    })

    return results

def summarize_proxy_results(results):
    """
    Summarize encounter fractions for each target parameter combination.
    """
    summary = {}
    for row in results:
        key = (
            row["n_target"],
            row["rho_h_target"],
            row["R"],
        )
        summary.setdefault(
            key,
            [],
        ).append(row["f_enc"])

    rows = []
    for (n_target, rho_h_target, R,), values in summary.items():
        values = np.asarray(
            values,
            dtype=float,
        )
        valid = np.isfinite(values)
        rows.append({
            "n_target":n_target,
            "rho_h_target":rho_h_target,
            "R":R,
            "f_enc_mean":np.mean(values[valid]),
            "f_enc_std":np.std(values[valid],ddof=1,),
            "n":np.sum(valid),
        })
    return rows
def fit_linear_proxy_for_R(results,R, n=None):
    """
    Test the hypothesis
        f_enc = C(R) * rho_h
    using a fit through the origin.
    If n is given, only that target tip-density
    subset is used.
    """
    data = [row for row in results if np.isclose(row["R"], R)]
    if n is not None:
        data = [row for row in data if np.isclose(row["n_target"],n)]
    rho_h = np.array([row["rho_h"] for row in data])
    f_enc = np.array([row["f_enc"] for row in data])
    valid = np.isfinite(rho_h) & np.isfinite(f_enc)
    rho_h = rho_h[valid]
    f_enc = f_enc[valid]
    if len(rho_h) == 0:
        return np.nan, np.nan
    # Fit through origin.
    C = np.sum(rho_h * f_enc) / np.sum(rho_h**2)
    prediction = C * rho_h
    ss_res = np.sum((f_enc - prediction)**2)
    ss_tot = np.sum((f_enc - np.mean(f_enc))**2)
    R2 = (1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan)
    return C, R2

def estimate_C_by_tip_density(results,R,n_values):
    """
    For one R, estimate C(R) for each n.
    Returns a dict mapping n -> (C, R2)
    """
    rows = []
    for n in n_values:
        C, R2 = fit_linear_proxy_for_R(results,R,n=n)
        rows.append({R: R, "n": n, "C": C, "R2": R2})
    return rows

def compute_normalized_encounter(results, R):
    data = [row for row in results if np.isclose(row["R"], R)]
    rows = []
    for row in data:
        rows.append({"R": R, "n": row["n"], "rho_h": row["rho_h"], "f_enc": row["f_enc"], "f_enc_over_rho_h": row["f_enc"] / row["rho_h"]})
    return rows


def compute_required_p(results, alpha, dt):
    rows = []
    for row in results:
        rho_h = row["rho_h"]
        f_enc = row["f_enc"]
        P_cont = 1.0 - np.exp(-alpha * rho_h * dt)
        p_required = P_cont / f_enc if f_enc > 0 else np.inf
        rows.append({**row, "P_cont": P_cont, "p_required": p_required, "p_valid": 0.0 <= p_required <= 1.0})
    return rows


def plot_f_enc_vs_rho_h(results,R,n_values):
    """
    For one R, plot f_enc vs rho_h.
    Different n values are shown separately.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = plt.cm.viridis(np.linspace(0,1,len(n_values)))
    for n, color in zip(n_values,colors):
        data = [row for row in results if np.isclose(row["R"], R) and np.isclose(row["n_target"],n)]
        rho = np.array([row["rho_h"] for row in data])
        f_enc = np.array([row["f_enc"] for row in data])
        ax.scatter(rho,f_enc,color=color,alpha=0.35,s=20)
        C, R2 = fit_linear_proxy_for_R(results, R, n=n)
        rho_fit = np.linspace(rho.min(),rho.max(),200)
        ax.plot(rho_fit,C * rho_fit,color=color,label=(rf"$n={n:g}$, " rf"$C={C:.3e}$, " rf"$R^2={R2:.3f}$"))
    ax.set_xlabel(r"realized hyphal length density $\rho_h$")
    ax.set_ylabel(r"$f_{\mathrm{enc}}$")
    ax.set_title(rf"Encounter fraction vs. $\rho_h$ " rf"for $R={R}$")
    ax.legend()
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()


def plot_normalized_encounter(results, R):
    data = [row for row in results if np.isclose(row["R"], R)]
    n = np.array([row["n"] for row in data])
    normalized = np.array([row["f_enc"] / row["rho_h"]for row in data ])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(n,normalized,alpha=0.25)
    ax.set_xlabel(r"tip density $n$")
    ax.set_ylabel(r"$f_{\mathrm{enc}}/\rho_h$")
    ax.set_title(rf"Test of density-independent $C(R)$, "rf"$R={R}$")
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    n_values = [1.5, 2.0, 2.5]
    rho_h_values = [2.0,2.25, 2.5]
    R_values = [0.01, 0.05, 0.1, 0.15, 0.2]

    n_reps = 50

    # Hyphal radius.
    #
    # R is interpreted as the distance from the tip
    # to the hyphal SURFACE.
    #
    # The actual geometric calculation therefore uses
    # R + r as the centerline distance.
    hyphal_radius = 0.0

    results = run_proxy_study(n_values=n_values, rho_h_values=rho_h_values, R_values=R_values, n_reps=n_reps, box=5.0, hyphal_radius=hyphal_radius)

    print("\n" + "=" * 70)
    print("REALIZED DENSITIES")
    print("=" * 70)

    for n in n_values:
        data = [row for row in results if np.isclose(row["n_target"], n)]
        errors = np.array([row["rho_h_relative_error"] for row in data])
        print(f"n = {n:5.1f}: mean relative rho_h error = {100 * np.mean(errors):+.2f}%")

    print("\n" + "=" * 70)
    print("LINEAR ENCOUNTER-DENSITY RELATION")
    print("=" * 70)

    C_results = []
    for R in R_values:
        print(f"\nR = {R}")
        rows = estimate_C_by_tip_density(results, R, n_values)
        C_results.extend(rows)
        for row in rows:
            print(f"  n = {row['n']:5.1f} : C = {row['C']:.4e}, R^2 = {row['R2']:.4f}")

    print("\n" + "=" * 70)
    print("NORMALIZED ENCOUNTER")
    print("=" * 70)

    for R in R_values:
        data = compute_normalized_encounter(results, R)
        print(f"\nR = {R}")
        for n in n_values:
            values = np.array([row["f_enc_over_rho_h"] for row in data if np.isclose(row["n"], n)])
            if len(values) > 0:
                print(f"  n = {n:5.1f}: mean C = {np.mean(values):.4e}, std = {np.std(values, ddof=1):.4e}")

    # -------------------------------------------------------------------------
    # 5. Optional continuum-to-discrete parameter bridge
    #
    #    Insert the actual alpha and dt used in the continuum model.
    # -------------------------------------------------------------------------

    alpha = 0.04
    dt = 1.0
    
    results_with_p = compute_required_p(results,alpha=alpha, dt=dt, )
    
    # =========================================================================
    # 1. REALIZED HYphal LENGTH DENSITY
    # =========================================================================

    print("\n" + "=" * 80)
    print("1. REALIZED HYPHAL LENGTH DENSITY")
    print("=" * 80)

    for rho_h_target in rho_h_values:
        data = [row for row in results if np.isclose(row["rho_h_target"],rho_h_target,)]
        errors = np.array([row["rho_h_relative_error"]for row in data])
        print(f"rho_h target = {rho_h_target:6.1f} : " f"realized relative error = " f"{100 * errors.mean():+.2f}% " f"+/- {100 * errors.std(ddof=1):.2f}%")

    # =========================================================================
    # 2. MAIN RESULT: REQUIRED p FOR EACH R
    # =========================================================================

    # print("\n" + "=" * 80)
    # print("2. REQUIRED CONDITIONAL FUSION PROBABILITY p")
    # print("=" * 80)

    # for R in R_values:
    #     data = [row for row in results_with_p if np.isclose(row["R"],R,)]
    #     p_values = np.array([row["p_required"] for row in data if np.isfinite(row["p_required"])])
    #     print(f"\nR = {R:.3f}")
    #     print("-" * 80)
    #     print(f"mean p= {np.mean(p_values):.5f}")
    #     print(f"median p     = {np.median(p_values):.5f}")
    #     print(f"std p        = {np.std(p_values, ddof=1):.5f}")
    #     print(f"min p        = {np.min(p_values):.5f}")
    #     print(f"max p        = {np.max(p_values):.5f}")
    #     print(f"fraction p>1 = " f"{np.mean(p_values > 1):.3f}")
    #     print(f"fraction p<0 = "f"{np.mean(p_values < 0):.3f}")

    # =========================================================================
    # 3. p AS A FUNCTION OF rho_h
    # =========================================================================

    # print("\n" + "=" * 80)
    # print("3. DEPENDENCE OF p ON rho_h")
    # print("=" * 80)

    # for R in R_values:
    #     print(f"\nR = {R:.3f}")
    #     print("-" * 80)

    #     for rho_h_target in rho_h_values:
    #         data = [row for row in results_with_p if np.isclose(row["R"], R) and np.isclose(row["rho_h_target"], rho_h_target)]
    #         p_values = np.array([row["p_required"] for row in data if np.isfinite(row["p_required"])])

    #         if len(p_values) == 0:
    #             continue

    #         print(f"rho_h = {rho_h_target:6.1f} : p = {np.mean(p_values):.5f} +/- {np.std(p_values, ddof=1):.5f}")

    # =========================================================================
    # 4. p AS A FUNCTION OF n
    # =========================================================================

    print("\n" + "=" * 80)
    print("4. DEPENDENCE OF p ON n")
    print("=" * 80)

    for R in R_values:
        print(f"\nR = {R:.3f}")
        print("-" * 80)
        for n_target in n_values:
            data = [row for row in results_with_p if (np.isclose(row["R"],R,) and np.isclose(row["n_target"], n_target,) )]
            p_values = np.array([row["p_required"]for row in data if np.isfinite(row["p_required"])])
            if len(p_values) == 0:
                continue
            print(f"n = {n_target:6.1f} : "f"p = {np.mean(p_values):.5f} "f"+/- {np.std(p_values, ddof=1):.5f}")

    # =========================================================================
    # 5. EMPIRICAL C(R)
    #
    # This is secondary: it explains the approximate value of p.
    # =========================================================================

    # print("\n" + "=" * 80)
    # print("5. EMPIRICAL C(R)")
    # print("=" * 80)

    # for R in R_values:
    #     C_rows = estimate_C_by_tip_density(results, R, n_values)
    #     C_values = np.array([row["C"] for row in C_rows if np.isfinite(row["C"])])
    #     R2_values = np.array([row["R2"] for row in C_rows if np.isfinite(row["R2"])])

    #     print(f"\nR = {R:.3f}")
    #     print(f"mean C = {np.mean(C_values):.5e}")
    #     print(f"std C = {np.std(C_values, ddof=1):.5e}")
    #     print(f"mean R² = {np.mean(R2_values):.4f}")

    #     C_theory = compute_theoretical_C(R, hyphal_radius=hyphal_radius)
    #     print(f"C_theory = {C_theory:.5e}")


    # =========================================================================
    # 6. DIRECT COMPARISON WITH THE APPROXIMATE FORMULA
    #
    # p ≈ alpha * dt / C(R)
    #
    # This approximation is only expected to work when
    # alpha * rho_h * dt << 1 and f_enc is approximately linear.
    # =========================================================================

    # print("\n" + "=" * 80)
    # print("6. p FROM C(R) VS. DIRECTLY REQUIRED p")
    # print("=" * 80)

    # for R in R_values:
    #     C_rows = estimate_C_by_tip_density(results, R, n_values)
    #     C_values = np.array([row["C"] for row in C_rows if np.isfinite(row["C"])])
    #     C_emp = np.mean(C_values)
    #     p_from_C = alpha * dt / C_emp

    #     data = [row for row in results_with_p if np.isclose(row["R"], R)]
    #     p_required = np.array([row["p_required"] for row in data if np.isfinite(row["p_required"])])

    #     print(f"\nR = {R:.3f}")
    #     print(f"p ≈ alpha*dt/C = {p_from_C:.5f}")
    #     print(f"direct mean p = {np.mean(p_required):.5f}")
    #     print(f"direct std p = {np.std(p_required, ddof=1):.5f}")

    # for row in results_with_p[:20]:
    #     print(f"R={row['R']:.3f}, rho_h={row['rho_h']:.2f}, f_enc={row['f_enc']:.5f}, P_cont={row['P_cont']:.5f}, p={row['p_required']:.2f}")


    # -------------------------------------------------------------------------
    # 6. Optional diagnostic plots
    # -------------------------------------------------------------------------

    # Uncomment only if visual inspection is useful.

    # for R in R_values:
    #     plot_f_enc_vs_rho_h(
    #         results,
    #         R,
    #         n_values,
    #     )
    
    #     plot_normalized_encounter(
    #         results,
    #         R,
    #     )
