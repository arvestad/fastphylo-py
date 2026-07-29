"""
Tree-comparison utilities for the reconstruction-accuracy benchmark:
Robinson-Foulds distance (topology) and leaf-pairwise path distances
(branch lengths). Both are keyed by leaf NAME rather than vertex ID, since
two trees built independently (e.g. the true tree and a reconstruction)
number their internal vertices differently.
"""
from __future__ import annotations

from collections import defaultdict


def _adjacency(edges):
    adj = defaultdict(list)
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))
    return adj


def bipartitions(edges, leaf_name_by_id, ref_name):
    """Non-trivial bipartitions induced by internal edges, as leaf-name
    frozensets. Each split is canonicalised to the side NOT containing
    ``ref_name``, so splits from two different trees over the same leaf
    set are directly comparable.
    """
    adj = _adjacency(edges)
    leaf_ids = set(leaf_name_by_id)
    all_names = frozenset(leaf_name_by_id.values())
    root = next(iter(leaf_ids))

    # Single DFS: `order` is parent-before-child, so reversing it gives a
    # valid bottom-up processing order for subtree-leafset aggregation.
    parent = {root: None}
    order = [root]
    stack = [root]
    while stack:
        cur = stack.pop()
        for nb, _w in adj[cur]:
            if nb not in parent:
                parent[nb] = cur
                order.append(nb)
                stack.append(nb)

    children = defaultdict(list)
    for node, p in parent.items():
        if p is not None:
            children[p].append(node)

    leafset = {}
    for node in reversed(order):
        if node in leaf_ids:
            leafset[node] = frozenset((leaf_name_by_id[node],))
        else:
            s = frozenset()
            for c in children[node]:
                s = s | leafset[c]
            leafset[node] = s

    splits = set()
    for node in order:
        if parent[node] is None:
            continue  # root has no parent edge
        side_names = leafset[node]
        if len(side_names) <= 1 or len(side_names) >= len(all_names) - 1:
            continue  # trivial (pendant-edge) split
        if ref_name in side_names:
            side_names = all_names - side_names
        splits.add(side_names)
    return splits


def rf_distance(edges_a, leaves_a, edges_b, leaves_b) -> float:
    """Normalised Robinson-Foulds distance in [0, 1] between two unrooted
    binary trees over the same leaf set (leaves matched by name)."""
    names = sorted(leaves_a.values())
    ref = names[0]
    n = len(names)
    sa = bipartitions(edges_a, leaves_a, ref)
    sb = bipartitions(edges_b, leaves_b, ref)
    max_splits = max(2 * (n - 3), 1)
    return len(sa ^ sb) / max_splits


def path_distances(edges, leaf_name_by_id) -> dict[tuple[str, str], float]:
    """All leaf-pairwise path (patristic) distances, keyed by a
    (name_i, name_j) tuple with name_i < name_j."""
    adj = _adjacency(edges)
    leaf_ids = sorted(leaf_name_by_id)
    out = {}
    for u in leaf_ids:
        dist = {u: 0.0}
        stack = [u]
        while stack:
            cur = stack.pop()
            for nb, w in adj[cur]:
                if nb not in dist:
                    dist[nb] = dist[cur] + w
                    stack.append(nb)
        for v in leaf_ids:
            if v <= u:
                continue
            ni, nj = leaf_name_by_id[u], leaf_name_by_id[v]
            if ni > nj:
                ni, nj = nj, ni
            out[(ni, nj)] = dist[v]
    return out
