"""
Generate random unrooted binary tree topologies with i.i.d. edge lengths,
for use in reconstruction-accuracy benchmarks.

Leaf vertex IDs are 0..n-1 (named "L0".."L{n-1}"); internal IDs start at n,
matching the convention used by fastphylo.Tree.
"""
from __future__ import annotations


def random_topology(n: int, rng) -> list[list[int]]:
    """Build a random unrooted binary tree topology on n leaves by
    sequential random edge subdivision (each new leaf is attached, via a
    fresh internal node, to a uniformly random existing edge).

    Returns a list of mutable [u, v] pairs (undirected edges, lengths not
    yet assigned).
    """
    if n < 3:
        raise ValueError("random_topology requires n >= 3")

    next_internal = n
    center = next_internal
    next_internal += 1
    edges = [[0, center], [1, center], [2, center]]

    for leaf in range(3, n):
        idx = int(rng.integers(len(edges)))
        u, v = edges[idx]
        w = next_internal
        next_internal += 1
        edges[idx] = [u, w]
        edges.append([w, v])
        edges.append([w, leaf])

    return edges


def random_tree(n: int, rng, lo: float = 0.02, hi: float = 0.1):
    """Random unrooted binary tree: topology + i.i.d. Uniform(lo, hi) edge
    lengths.

    Returns (edges, leaves) in the same format as fastphylo.Tree:
      edges  — list of (u, v, length)
      leaves — dict[int, str] mapping leaf vertex id -> name "L{i}"
    """
    topo = random_topology(n, rng)
    edges = [(u, v, float(rng.uniform(lo, hi))) for u, v in topo]
    leaves = {i: f"L{i}" for i in range(n)}
    return edges, leaves
