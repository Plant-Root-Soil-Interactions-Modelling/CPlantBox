import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import networkx as nx
from scipy.spatial import Delaunay, cKDTree
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist

import AMFAnalysis as amf

def getAnastomosisMerges(mycp):
    # {tip node id: node id it merged with}, at an anastomosis the hypha gets a new node at the position
    # of the target node (see MycorrhizalPlant::simulateAnastomosis), both are the same point in the network
    merges = {}
    for h in mycp.getOrgans(pb.hyphae):
        mergeID = int(h.getParameter("mergePointID"))
        if mergeID > -1:
            merges[h.getNodeId(h.getNumberOfNodes() - 1)] = mergeID
    for tip in merges: # a target might itself be a merged tip
        target = merges[tip]
        while target in merges and merges[target] != tip:
            target = merges[target]
        merges[tip] = target
    return merges

def getNetworkX(mycp, geometry = None, hyphae_only = True):
    # undirected networkx graph of the hyphal network (inside geometry, e.g. small_hyphae_dish)
    # nodes: node ids of the plant (cut nodes from crop get new ids), attribute "pos" (x, y, z)
    # edges: one per segment, attributes "length" (3D), "length_xy" (projected to x-y plane) and all per segment data of getMycSegmentAnalyser
    # anastomoses: the tip node is merged with the node it fused with, so the hyphal trees are connected
    ana = amf.getMycSegmentAnalyser(mycp)
    if hyphae_only:
        ana.filter("organType", 5)
    if geometry is not None:
        ana.crop(geometry) # no pack(), it would renumber the nodes and the anastomosis ids would not fit anymore
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes])
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int)
    inside = set(segs.ravel().tolist())
    merges = {tip: target for tip, target in getAnastomosisMerges(mycp).items() if target in inside} # fusions with hyphae outside geometry are not connected
    segs_merged = np.vectorize(lambda i: merges.get(i, i), otypes=[int])(segs) if (merges and len(segs) > 0) else segs
    data = {"length": np.array(ana.getParameter("length"))}
    data["length_xy"] = np.linalg.norm(nodes[segs[:, 0], :2] - nodes[segs[:, 1], :2], axis=1) if len(segs) > 0 else np.zeros(0)
    for key in ana.data.keys():
        data[key] = np.array(ana.data[key])
    G = nx.Graph(n_plant_nodes = mycp.getNumberOfNodes()) # node ids >= n_plant_nodes are cut nodes made by crop
    for i in np.unique(segs_merged):
        G.add_node(int(i), pos = tuple(nodes[i].tolist()))
    keys = list(data.keys())
    for k, (u, v) in enumerate(segs_merged):
        if u != v: # an anastomosis segment of zero length could become a self loop
            G.add_edge(int(u), int(v), **{key: float(data[key][k]) for key in keys})
    return G

def getPlanarNetwork(mycp, x_min, eps = 1.e-9, t_max = None):
    # 2D graph of the hyphal network as extracted from top view images (Oyarte Galvez et al. 2025: planar skeleton, crossings are
    # degree 4 nodes): all hyphal segments with x >= x_min (fungal side, segments across x = x_min are cut there) projected to the
    # x-y plane, every crossing of two segments becomes a node of both. Anastomoses as in getNetworkX. Node ids: plant node ids,
    # new nodes (cuts, crossings) >= graph attribute "n_plant_nodes". Edges: "length" (3D, split proportionally), "length_xy",
    # "width_len" (width × length_xy, width = 2 × radius of the segment: all segments of a hyphal organ have the same radius,
    # so the mean width of a simplified edge is width_len / length_xy) and "bas_len" (length_xy if the segment belongs to a BAS)
    # t_max [d]: the network as it was at that time (only segments created until then), from the final plant
    ana = amf.getMycSegmentAnalyser(mycp)
    ana.filter("organType", 5)
    if t_max is not None:
        ana.filter("creationTime", -1.e9, t_max + 1.e-9) # no pack(): node ids stay those of the plant
    nodes = np.array([[n.x, n.y, n.z] for n in ana.nodes])
    segs = np.array([[s.x, s.y] for s in ana.segments], dtype=int).reshape(-1, 2)
    seg_r = np.array(ana.getParameter("radius"), dtype=float)
    seg_bas = np.array(ana.getParameter("subType")).astype(int) == 2      # BAS: subType 2 of TwoHyphaePlusBAS.xml
    n_plant = max(len(nodes), mycp.getNumberOfNodes())
    merges = getAnastomosisMerges(mycp)
    if merges and len(segs) > 0:
        segs = np.vectorize(lambda i: merges.get(i, i), otypes=[int])(segs)
    nonzero = segs[:, 0] != segs[:, 1]
    segs, seg_r, seg_bas = segs[nonzero], seg_r[nonzero], seg_bas[nonzero]
    pos = [nodes]
    # cut at x = x_min
    xa, xb = nodes[segs[:, 0], 0], nodes[segs[:, 1], 0]
    keep = (xa >= x_min) & (xb >= x_min)
    cross = (xa >= x_min) != (xb >= x_min)
    cut = segs[cross]
    inside = np.where(nodes[cut[:, 0], 0] >= x_min, cut[:, 0], cut[:, 1])
    a, b = nodes[cut[:, 0]], nodes[cut[:, 1]]
    p_cut = a + ((x_min - a[:, 0]) / (b[:, 0] - a[:, 0]))[:, None] * (b - a)
    cut_ids = n_plant + np.arange(len(cut))
    pos.append(p_cut)
    segs = np.vstack([segs[keep], np.column_stack([inside, cut_ids])]) if len(cut) else segs[keep]
    seg_r = np.concatenate([seg_r[keep], seg_r[cross]])
    seg_bas = np.concatenate([seg_bas[keep], seg_bas[cross]])
    pos = np.vstack(pos)
    next_id = n_plant + len(cut)
    # crossings in 2D: candidate pairs from the midpoints (two segments of length <= L can only cross if their midpoints are closer than L)
    p, q = pos[segs[:, 0]], pos[segs[:, 1]]
    seg_len = np.linalg.norm(q - p, axis=1)
    splits = [[] for _ in range(len(segs))] # (parameter along the segment, node id)
    if len(segs) > 1:
        L = np.linalg.norm((q - p)[:, :2], axis=1).max()
        pairs = cKDTree(0.5 * (p + q)[:, :2]).query_pairs(L, output_type="ndarray")
        i, j = pairs[:, 0], pairs[:, 1]
        share = (segs[i, 0] == segs[j, 0]) | (segs[i, 0] == segs[j, 1]) | (segs[i, 1] == segs[j, 0]) | (segs[i, 1] == segs[j, 1])
        i, j = i[~share], j[~share]
        r, s = (q - p)[i, :2], (q - p)[j, :2]
        d = p[j, :2] - p[i, :2]
        den = r[:, 0] * s[:, 1] - r[:, 1] * s[:, 0]
        ok = np.abs(den) > eps # parallel segments do not cross
        t = np.where(ok, (d[:, 0] * s[:, 1] - d[:, 1] * s[:, 0]) / np.where(ok, den, 1.), -1.)
        u = np.where(ok, (d[:, 0] * r[:, 1] - d[:, 1] * r[:, 0]) / np.where(ok, den, 1.), -1.)
        hit = ok & (t > eps) & (t < 1 - eps) & (u > eps) & (u < 1 - eps)
        i, j, t, u = i[hit], j[hit], t[hit], u[hit]
        new = p[i] + t[:, None] * (q[i] - p[i]) # z of the first segment, only x-y is used in 2D
        ids = next_id + np.arange(len(i))
        pos = np.vstack([pos, new])
        for k in range(len(i)):
            splits[i[k]].append((t[k], ids[k]))
            splits[j[k]].append((u[k], ids[k]))
    G = nx.Graph(n_plant_nodes = n_plant)
    def add_edge(n1, n2, frac, k):
        if n1 == n2:
            return
        l_xy = float(np.linalg.norm(pos[n1, :2] - pos[n2, :2]))
        G.add_edge(int(n1), int(n2), length = float(frac * seg_len[k]), length_xy = l_xy,
                   width_len = 2. * float(seg_r[k]) * l_xy, bas_len = l_xy if seg_bas[k] else 0.)
    for k, (n1, n2) in enumerate(segs):
        chain = [(0., n1)] + sorted(splits[k]) + [(1., n2)]
        for (t1, m1), (t2, m2) in zip(chain[:-1], chain[1:]):
            add_edge(m1, m2, t2 - t1, k)
    for n in G.nodes:
        G.nodes[n]["pos"] = tuple(pos[n].tolist())
    return G

BAS_IMAGE_CRITERIA = {"max_len": 400., "max_len_tip": 1000., "max_width": 7., "max_width_len": 9000.}   # [µm], [µm^2]

def basShareImage(mycp, x_min, criteria = BAS_IMAGE_CRITERIA, t_max = None):
    # share of BAS in the total hyphal length as Oyarte Galvez et al. (2025, Methods) classify the edges of the image graph:
    # an edge (between junctions, crossings and ends of the planar skeleton, getPlanarNetwork + simplifyGraph) is a BAS if
    # (1) length < 400 um, or (2) length < 1000 um and one end point is a tip (degree 1, not a cut at x_min), or
    # (3) mean width < 7 um and width x length < 9000 um^2. Lengths in the x-y plane, width = 2 x radius of the segments
    # (one radius per hyphal organ in the model). Returns the share, the share by sub type in the same graph (for comparison)
    # and the numbers of edges. t_max [d]: the network at that time (getPlanarNetwork)
    G = getPlanarNetwork(mycp, x_min, t_max = t_max)
    S = simplifyGraph(G, weights = ("length", "length_xy", "width_len", "bas_len"))
    n_plant = S.graph["n_plant_nodes"]
    deg = dict(S.degree())
    total, bas, bas_subtype, n_bas, n_edges = 0., 0., 0., 0, 0
    for u, v, d in S.edges(data = True):
        L = d["length_xy"] * 1.e4 # [um]
        if L <= 0:
            continue
        w = d["width_len"] / d["length_xy"] * 1.e4 # mean width [um]
        tip = any(deg[n] == 1 and n < n_plant for n in (u, v))
        is_bas = (L < criteria["max_len"]) or (L < criteria["max_len_tip"] and tip) or \
                 (w < criteria["max_width"] and w * L < criteria["max_width_len"])
        total += L
        bas_subtype += d["bas_len"] * 1.e4
        n_edges += 1
        if is_bas:
            bas += L
            n_bas += 1
    return {"share": bas / total if total > 0 else np.nan, "share_subtype": bas_subtype / total if total > 0 else np.nan,
            "n_edges": n_edges, "n_bas_edges": n_bas, "length_xy_um": total}

def lengthWeight(dim):
    # edge attribute used as weight: 3D length or length projected to the x-y plane (dish seen from above)
    return "length" if dim == 3 else "length_xy"

def simplifyGraph(G, weights = ("length", "length_xy")):
    # multigraph with only junctions and ends of G (degree != 2), each chain of degree 2 nodes becomes one edge
    # with summed weights and the number of segments ("segments"); a multigraph keeps parallel chains and loops,
    # so total length is unchanged and shortest paths between the kept nodes are the same as in G
    keep = {n for n, d in G.degree() if d != 2}
    for comp in nx.connected_components(G): # closed loops without any junction
        if not any(n in keep for n in comp):
            keep.add(min(comp))
    S = nx.MultiGraph(**G.graph)
    for n in keep:
        S.add_node(n, **G.nodes[n])
    visited = set()
    for n in keep:
        for nb in G.neighbors(n):
            if frozenset((n, nb)) in visited:
                continue
            sums = {w: 0. for w in weights}
            segments = 0
            prev, cur = n, nb
            while True:
                visited.add(frozenset((prev, cur)))
                for w in weights:
                    sums[w] += G.edges[prev, cur][w]
                segments += 1
                if cur in keep:
                    break
                prev, cur = cur, next(x for x in G.neighbors(cur) if x != prev) # cur has degree 2
            S.add_edge(n, cur, segments = segments, **sums)
    return S

def getReferenceNetworks(S, dim = 3):
    # Delaunay triangulation and Euclidean minimum spanning tree on the nodes of S (use simplifyGraph first),
    # dim = 2: triangulation of the x-y positions (dish seen from above), dim = 3: tetrahedra in x, y, z
    # edges get "length" and "length_xy" (straight lines), the MST minimises lengthWeight(dim)
    nodes = list(S.nodes)
    pos = np.array([S.nodes[n]["pos"] for n in nodes])
    tri = Delaunay(pos[:, :dim], qhull_options = "QJ") # QJ: robust against duplicate and coplanar points
    D = nx.Graph()
    for n in nodes:
        D.add_node(n, **S.nodes[n])
    for simplex in tri.simplices:
        for a in range(len(simplex)):
            for b in range(a + 1, len(simplex)):
                i, j = simplex[a], simplex[b]
                D.add_edge(nodes[i], nodes[j], length = float(np.linalg.norm(pos[i] - pos[j])), length_xy = float(np.linalg.norm(pos[i, :2] - pos[j, :2])))
    M = nx.minimum_spanning_tree(D, weight = lengthWeight(dim)) # the Euclidean MST is always part of the Delaunay graph
    return D, M

def weightedEfficiency(G, dim = 3):
    # global efficiency of Latora & Marchiori (2001): sum of 1/(shortest path length) over all node pairs,
    # normalised by sum of 1/(straight line distance), 1 = every pair connected by a straight line;
    # unconnected pairs contribute 0, pairs at the same position are skipped
    
    weight = lengthWeight(dim)
    nodes = list(G.nodes)
    pos = np.array([G.nodes[n]["pos"] for n in nodes])[:, :dim]
    eucl = pdist(pos)
    inv_eucl = 2 * np.sum(1. / eucl[eucl > 0])
    inv_path = 0.
    for s in nodes:
        dist = nx.single_source_dijkstra_path_length(G, s, weight = weight)
        inv_path += sum(1. / d for t, d in dist.items() if t != s and d > 0)
    return float(inv_path / inv_eucl)

def getSourcesTargets(G, barrier_thickness, tol = 1.e-5):
    # sources: nodes where hyphae cross into the hyphal compartment, i.e. the nodes crop(small_hyphae_dish) created within the
    # barrier (x <= barrier_thickness/2): on the cut plane x = -barrier_thickness/2, or with makedishes(hyphae_layer = ...)
    # also at the lower side of the layer inside the opening (degree 1); targets: all other nodes of degree 1 (tips and ends)
    ends = [n for n, d in G.degree() if d == 1]
    if "n_plant_nodes" in G.graph:
        sources = [n for n in ends if n >= G.graph["n_plant_nodes"] and G.nodes[n]["pos"][0] <= barrier_thickness / 2 + tol]
    else: # graph without this information, only the cut plane
        sources = [n for n in ends if abs(G.nodes[n]["pos"][0] + barrier_thickness / 2) < tol]
    s = set(sources)
    targets = [n for n in ends if n not in s]
    return sources, targets

def weightedEfficiencySubset(G, sources, targets = None, dim = 3):
    # efficiency of Latora & Marchiori (2001) only for the pairs (source, target), or for all pairs within sources if targets is None:
    # sum of 1/(shortest path length) normalised by sum of 1/(straight line distance), unconnected pairs contribute 0
    
    weight = lengthWeight(dim)
    sources = list(sources)
    if targets is None:
        targets = sources
    else:
        targets = [t for t in targets if t not in set(sources)]
    pos_s = np.array([G.nodes[n]["pos"] for n in sources])[:, :dim]
    pos_t = np.array([G.nodes[n]["pos"] for n in targets])[:, :dim]
    eucl = cdist(pos_s, pos_t) # pairs within one set are counted twice, as are the paths below
    inv_eucl = np.sum(1. / eucl[eucl > 0])
    inv_path = 0.
    for s in sources:
        dist = nx.single_source_dijkstra_path_length(G, s, weight = weight)
        inv_path += sum(1. / dist[t] for t in targets if t != s and t in dist and dist[t] > 0)
    return float(inv_path / inv_eucl)

def getTerminalRootNodes(mycp, S, barrier_thickness, root_dist = 0.5, t_max = None):
    # node sets as in Oyarte Galvez et al. (2025), both restricted to the nodes of S (e.g. simplifyGraph(getNetworkX(mycp, small_hyphae_dish))):
    # R (root nodes): all nodes within root_dist [cm] (5 mm) of the barrier face on the hyphal side (x = barrier_thickness/2),
    #                 i.e. sites of crossing from the root compartment (5 mm as defined by Oyarte Galvez et al. 2025)
    # T (terminal nodes): growing tips (last node of active non-BAS hyphae; also tips frozen because an ancestor died, see
    #                 sa_setup/model_hook.frozen_diagnostics) and attachment sites of BAS (first node of a BAS, subType 2,
    #                 whose parent is no BAS); nodes in R are not in T
    # t_max [d]: the sets at that time (S built with the same t_max): the tip of a runner is its last node created until t_max,
    #                 counted if the hypha still grew after t_max or is active at the end; BAS attachment sites created until t_max
    x_root = barrier_thickness / 2 + root_dist
    R = [n for n in S.nodes if S.nodes[n]["pos"][0] <= x_root]
    in_R = set(R)
    T = set()
    for h in mycp.getOrgans(pb.hyphae):
        if t_max is not None and h.getNodeCT(0) > t_max + 1.e-9:
            continue
        if int(h.getParameter("subType")) == 2:
            p = h.getParent()
            if not (int(p.getParameter("organType")) == 5 and int(p.getParameter("subType")) == 2):
                T.add(h.getNodeId(0))
        elif t_max is None:
            if h.isActive():
                T.add(h.getNodeId(h.getNumberOfNodes() - 1))
        else:
            n = h.getNumberOfNodes()
            cts = np.array([h.getNodeCT(i) for i in range(n)])
            k = int(np.searchsorted(cts, t_max + 1.e-9, side="right")) - 1
            if k >= 0 and (cts[-1] > t_max + 1.e-9 or h.isActive()):
                T.add(h.getNodeId(k))
    T = [n for n in T if n in S and n not in in_R]
    return T, R

def getEfficiencies(S, T, R, dim = 3):
    # (i) global efficiency: transport between arbitrary pairs of the transport relevant nodes T and R
    # (ii) root efficiency: transport between the terminal nodes T and the root nodes R
    return weightedEfficiencySubset(S, list(T) + list(R), None, dim), weightedEfficiencySubset(S, T, R, dim)

def relativeRootEfficiency(mycp, geometry, barrier_thickness, dim = 2, root_dist = 0.5, planar = True, t_max = None,
                           global_eff = False, max_sources = 300, seed = 0):
    # relative root efficiency as in Oyarte Galvez et al. (2025): E_r of the network between the terminal nodes T and the root nodes R
    # (getTerminalRootNodes), normalised between the minimum spanning tree and the Delaunay triangulation on the same nodes:
    # E_rel = (E_r - E_r,MST) / (E_r,DT - E_r,MST).
    # planar = True: graph as from the top view images of the experiment (getPlanarNetwork, fungal side x >= -barrier_thickness/2,
    # crossings are nodes), always 2D; planar = False: 3D graph cropped to geometry (small_hyphae_dish), dim = 2 or 3.
    # Dijkstra runs from the smaller of the two node sets (the efficiency is symmetric). Returns a dict with all four values and the set sizes
    # t_max [d] (planar only): the network as it was at that time
    # global_eff: additionally, on the same graphs (Oyarte Galvez et al. 2025, Ê_g and Ĉ; the normalisations are assumed,
    # the paper gives no details): the global efficiency between all pairs of T ∪ R, E_g_rel = (E_g - E_g,MST) / (E_g,DT - E_g,MST),
    # and the cost C_rel = (L - L_MST) / (L_DT - L_MST) with the total length L (x-y in 2D). For E_g at most max_sources nodes
    # of T ∪ R are used (random subset with seed, the same for the network and both references), otherwise it gets slow
    if planar:
        S, dim = simplifyGraph(getPlanarNetwork(mycp, -barrier_thickness / 2, t_max = t_max)), 2
    else:
        if t_max is not None:
            raise ValueError("relativeRootEfficiency: t_max only for planar = True")
        S = simplifyGraph(getNetworkX(mycp, geometry))
    T, R = getTerminalRootNodes(mycp, S, barrier_thickness, root_dist, t_max = t_max)
    out = {"n_T": len(T), "n_R": len(R), "E_r": np.nan, "E_r_DT": np.nan, "E_r_MST": np.nan, "E_rel": np.nan}
    if global_eff:
        out.update({k: np.nan for k in ("E_g", "E_g_DT", "E_g_MST", "E_g_rel", "L", "L_DT", "L_MST", "C_rel")})
    if len(T) == 0 or len(R) == 0 or S.number_of_nodes() < dim + 2:
        return out
    a, b = (T, R) if len(T) <= len(R) else (R, T)
    D, M = getReferenceNetworks(S, dim)
    out["E_r"] = weightedEfficiencySubset(S, a, b, dim)
    out["E_r_DT"] = weightedEfficiencySubset(D, a, b, dim)
    out["E_r_MST"] = weightedEfficiencySubset(M, a, b, dim)
    span = out["E_r_DT"] - out["E_r_MST"]
    out["E_rel"] = (out["E_r"] - out["E_r_MST"]) / span if span > 0 else np.nan
    if global_eff:
        w = lengthWeight(dim)
        nodes = list(T) + list(R)
        if max_sources is not None and len(nodes) > max_sources:
            nodes = list(np.random.default_rng(seed).choice(np.array(nodes), max_sources, replace = False))
        out["E_g"] = weightedEfficiencySubset(S, nodes, None, dim)
        out["E_g_DT"] = weightedEfficiencySubset(D, nodes, None, dim)
        out["E_g_MST"] = weightedEfficiencySubset(M, nodes, None, dim)
        span = out["E_g_DT"] - out["E_g_MST"]
        out["E_g_rel"] = (out["E_g"] - out["E_g_MST"]) / span if span > 0 else np.nan
        out["L"] = float(sum(d[w] for _, _, d in S.edges(data = True)))
        out["L_DT"] = float(sum(d[w] for _, _, d in D.edges(data = True)))
        out["L_MST"] = float(sum(d[w] for _, _, d in M.edges(data = True)))
        span = out["L_DT"] - out["L_MST"]
        out["C_rel"] = (out["L"] - out["L_MST"]) / span if span > 0 else np.nan
    return out

def betweennessSubset(G, sources, targets, dim = 3, edges = False):
    # fraction of all shortest (source, target) paths running through each node (or edge, if edges = True);
    # networkx counts every path in an undirected graph with 0.5 and its normalized=True divides by (n-1)(n-2), so rescale here
    weight = lengthWeight(dim)
    if edges:
        bc = nx.edge_betweenness_centrality_subset(G, sources, targets, normalized = False, weight = weight)
    else:
        bc = nx.betweenness_centrality_subset(G, sources, targets, normalized = False, weight = weight)
    scale = 2. / (len(sources) * len(targets)) if (sources and targets) else 0.
    return {k: v * scale for k, v in bc.items()}