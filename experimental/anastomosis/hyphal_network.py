"""
Random 3D hyphal / hyphal-tip network inside a unit cube.

Geometry model
--------------
- Hyphae are straight line segments, hyphal tips are the apical (growing) end
  of a segment.
- Hyphal tips are placed uniformly at random inside the cube to realize a
  target hyphal tip density 'rho_t' (tips per unit volume).
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


def _random_directions(n, rng):
    """n random directions, uniformly distributed on the unit sphere."""
    v = rng.normal(size=(n, 3))
    norms = np.linalg.norm(v, axis=1, keepdims=True)
    return v / norms


def _ray_box_exit(points, directions, box):
    """distance from each point (inside [0, box]^3) to the box boundary along directions"""
    t = np.full(directions.shape, np.inf)
    pos, neg = directions > 0, directions < 0
    t[pos] = (box - points[pos]) / directions[pos]
    t[neg] = (0.0 - points[neg]) / directions[neg]
    return t.min(axis=1)


def create_hyphal_network(rho_t, rho_h, box=1.0, seed=None):
    """
    Creates a random 3D hyphal network inside a cube of side length 'box' that
    realizes a hyphal tip density rho_t and a hyphal length density rho_h.

    Parameters:
    rho_t: hyphal tip density [number of apical tips per unit volume]
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

    # 1) place hyphal tips uniformly to realize rho_t
    n_tips = int(round(rho_t * volume))
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

    total_length = float(sum(np.linalg.norm(p1 - p0) for p0, p1 in segments))
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

    return {
        "tips": tips,
        "segments": segments,
        "is_tip_segment": np.array(is_tip_segment, dtype=bool),
        "tip_segment_index": np.arange(n_tips),  # segments[:n_tips] are the tip-anchored ones, in tip order
        "total_length": total_length,
        "volume": volume,
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


def count_tips_near_hyphae(network, R):
    """
    Number of hyphal tips whose minimal distance to some other hypha is <= R.
    R can be a scalar or an array-like of radii; distances are precomputed
    once (cached on first call, see compute_tip_neighbor_distances) so that
    each additional radius only costs an O(log n) lookup.
    """
    if "tip_min_distances_sorted" not in network:
        compute_tip_neighbor_distances(network)
    R = np.asarray(R, dtype=float)
    counts = np.searchsorted(network["tip_min_distances_sorted"], R, side="right")
    return counts if R.ndim > 0 else int(counts)


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    rho_t, rho_h = 20.0, 80.0  # tips per unit volume, hyphal length per unit volume
    net = create_hyphal_network(rho_t, rho_h, seed=1)

    print("hyphal tip density   target {:g}, realized {:g}".format(rho_t, len(net["tips"]) / net["volume"]))
    print("hyphal length density target {:g}, realized {:g}".format(rho_h, net["total_length"] / net["volume"]))
    print("number of hyphae: {:d} ({:d} anchored at a tip, {:d} crossing hyphae)".format(len(net["segments"]), int(net["is_tip_segment"].sum()), int((~net["is_tip_segment"]).sum())))

    compute_tip_neighbor_distances(net)
    radii = [0.01, 0.02, 0.05, 0.1]
    counts = count_tips_near_hyphae(net, radii)
    for R, c in zip(radii, counts):
        print("R = {:g}: {:d} of {:d} tips are within radius R of another hypha".format(R, int(c), len(net["tips"])))

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    lines = [(p0, p1) for p0, p1 in net["segments"]]
    colors = ["C0" if is_tip else "C1" for is_tip in net["is_tip_segment"]]
    ax.add_collection3d(Line3DCollection(lines, colors=colors, linewidths=1.0))
    if len(net["tips"]) > 0:
        ax.scatter(*net["tips"].T, color="k", s=10, label="hyphal tips")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_zlim(0, 1)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.legend()

    # fraction of tips within radius R of another hypha, for a few hyphal length densities
    # (mean +/- std over repeated random realizations of the network)
    R_range = np.linspace(0.0, 0.3, 100)
    rho_h_values = [40.0, 80.0, 160.0]
    n_reps = 50
    fig2, ax2 = plt.subplots()
    for rho_h_i, color in zip(rho_h_values, ["C0", "C1", "C2"]):
        fractions = np.empty((n_reps, len(R_range)))
        for rep in range(n_reps):
            net_i = create_hyphal_network(rho_t, rho_h_i, seed=rep)
            compute_tip_neighbor_distances(net_i)
            fractions[rep] = count_tips_near_hyphae(net_i, R_range) / len(net_i["tips"])
        mean, std = fractions.mean(axis=0), fractions.std(axis=0)
        ax2.plot(R_range, mean, color=color, label=r"$\rho_h$ = {:g}".format(rho_h_i))
        ax2.fill_between(R_range, mean - std, mean + std, color=color, alpha=0.2)
    ax2.set_xlabel("R")
    ax2.set_ylabel("fraction of tips within R of another hypha")
    ax2.legend()

    plt.show()
