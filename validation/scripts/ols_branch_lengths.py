"""
Unconstrained ordinary-least-squares (OLS) branch-length fitting for a
fixed tree topology, used to compute the minimum-evolution (ME) criterion
the way Gascuel (1997) does: fit a candidate topology's branch lengths to
a (noisy) distance matrix by OLS, then sum the fitted lengths to get that
topology's ME "tree length" -- the same procedure applied to both the true
tree and the reconstructed (NJ/BioNJ) tree so their lengths are comparable.

This intentionally has no non-negativity constraint (matches the classic
OLS-ME framework of Rzhetsky & Nei 1993), unlike fastphylo.fit_branch_lengths
(L1, non-negative), which serves a different purpose (Tree.to_newick() with
guaranteed-plottable lengths).
"""
from __future__ import annotations

from collections import deque

import numpy as np


def _path_edges(adj, edge_idx, src, dst):
    parent_edge = {}
    visited = {src}
    queue = deque([src])
    while queue:
        node = queue.popleft()
        if node == dst:
            path = []
            cur = dst
            while cur != src:
                par, eidx = parent_edge[cur]
                path.append(eidx)
                cur = par
            return path
        for nb in adj.get(node, []):
            if nb not in visited:
                visited.add(nb)
                parent_edge[nb] = (node, edge_idx[(node, nb)])
                queue.append(nb)
    raise ValueError(f"no path from {src} to {dst}")


def ols_tree_length(edges, leaf_name_by_id, dm_names, dm) -> float:
    """OLS-fit branch lengths for `edges` (topology only; existing lengths
    on `edges` are ignored) against distance matrix `dm` (indexed by
    position in `dm_names`, matched to tree leaves by name). Returns the
    sum of fitted lengths -- the topology's ME criterion value.
    """
    name_to_vertex = {name: vid for vid, name in leaf_name_by_id.items()}
    b = len(edges)
    edge_idx: dict[tuple, int] = {}
    adj: dict[int, list[int]] = {}
    for k, (u, v, _w) in enumerate(edges):
        edge_idx[(u, v)] = k
        edge_idx[(v, u)] = k
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)

    n = len(dm_names)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    P = np.zeros((m, b))
    d_obs = np.zeros(m)
    for k, (i, j) in enumerate(pairs):
        vi = name_to_vertex[dm_names[i]]
        vj = name_to_vertex[dm_names[j]]
        for eidx in _path_edges(adj, edge_idx, vi, vj):
            P[k, eidx] = 1.0
        d_obs[k] = dm[i, j]

    x, *_ = np.linalg.lstsq(P, d_obs, rcond=None)
    return float(np.sum(x))
