import sys; sys.path.append("../.."); sys.path.append("../../src/")
import numpy as np
import plantbox as pb
import networkx as nx
from scipy.spatial import Delaunay
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

def getTerminalRootNodes(mycp, S, barrier_thickness, root_dist = 0.5):
    # node sets as in Oyarte Galvez et al. (2025), both restricted to the nodes of S (e.g. simplifyGraph(getNetworkX(mycp, small_hyphae_dish))):
    # R (root nodes): all nodes within root_dist [cm] (5 mm) of the barrier face on the hyphal side (x = barrier_thickness/2),
    #                 i.e. sites of crossing from the root compartment
    # T (terminal nodes): growing tips (last node of active non-BAS hyphae) and attachment sites of BAS (first node of a BAS,
    #                 subType 3, whose parent is no BAS); nodes in R are not in T
    x_root = barrier_thickness / 2 + root_dist
    R = [n for n in S.nodes if S.nodes[n]["pos"][0] <= x_root]
    in_R = set(R)
    T = set()
    for h in mycp.getOrgans(pb.hyphae):
        if int(h.getParameter("subType")) == 2:
            p = h.getParent()
            if not (int(p.getParameter("organType")) == 5 and int(p.getParameter("subType")) == 2):
                T.add(h.getNodeId(0))
        elif h.isActive():
            T.add(h.getNodeId(h.getNumberOfNodes() - 1))
    T = [n for n in T if n in S and n not in in_R]
    return T, R

def getEfficiencies(S, T, R, dim = 3):
    # (i) global efficiency: transport between arbitrary pairs of the transport relevant nodes T and R
    # (ii) root efficiency: transport between the terminal nodes T and the root nodes R
    return weightedEfficiencySubset(S, list(T) + list(R), None, dim), weightedEfficiencySubset(S, T, R, dim)

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