"""Tree data structure and NJ/FNJ/BioNJ reconstruction functions."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from . import _fastphylo


# ---------------------------------------------------------------------------
# Protocol for duck-typed distance matrices
# ---------------------------------------------------------------------------

@runtime_checkable
class DistanceProtocol(Protocol):
    def __len__(self) -> int: ...
    def __getitem__(self, key: tuple[int, int]) -> float: ...
    def names(self) -> list[str]: ...


# ---------------------------------------------------------------------------
# Tree
# ---------------------------------------------------------------------------

class Tree:
    """Unrooted phylogenetic tree stored as an edge set.

    Leaf vertex IDs 0…N-1 match the distance matrix row order.
    Internal node IDs start at N.
    Branch lengths are -1.0 for NJ/FNJ (not computed by FastPhylo);
    BioNJ sets real values (possibly negative for short branches).
    """

    def __init__(
        self,
        edges: list[tuple[int, int, float]],
        leaves: dict[int, str],
        num_vertices: int,
    ) -> None:
        self.edges = edges
        self.leaves = leaves          # vertex_id -> accession
        self.num_vertices = num_vertices

    def to_newick(self) -> str:
        """Serialize to Newick format.

        Branch lengths are omitted when all edges carry the sentinel -1.0
        (NJ/FNJ). BioNJ branch lengths are always included, even if negative.
        """
        if not self.edges:
            if self.leaves:
                return next(iter(self.leaves.values())) + ";"
            return "();"

        # Build undirected adjacency list: node → [(neighbour, branch_len)]
        adj: dict[int, list[tuple[int, float]]] = {}
        for u, v, w in self.edges:
            adj.setdefault(u, []).append((v, w))
            adj.setdefault(v, []).append((u, w))

        # Root = node never appearing as the child (first element) of an edge
        children_set = {u for u, v, _w in self.edges}
        all_nodes = set(adj)
        roots = all_nodes - children_set
        root = next(iter(roots))

        include_lengths = not all(w == -1.0 for _u, _v, w in self.edges)

        def _subtree(node: int, parent: int | None) -> str:
            kids = [(nb, w) for nb, w in adj.get(node, []) if nb != parent]
            if not kids:
                return self.leaves.get(node, str(node))
            parts = []
            for child, branch_len in kids:
                sub = _subtree(child, node)
                if include_lengths:
                    parts.append(f"{sub}:{branch_len:.6g}")
                else:
                    parts.append(sub)
            inner = f"({','.join(parts)})"
            return inner

        return _subtree(root, None) + ";"

    def merge(self, other: "Tree") -> "Tree":
        """Return a new Tree that is the union of both edge sets.

        All vertex IDs in *other* are offset by ``self.num_vertices`` so that
        no IDs collide. Leaf accessions identify taxa across trees.
        """
        offset = self.num_vertices

        new_edges = list(self.edges) + [
            (u + offset, v + offset, w) for u, v, w in other.edges
        ]
        new_leaves = dict(self.leaves)
        new_leaves.update({vid + offset: name for vid, name in other.leaves.items()})
        new_num_vertices = self.num_vertices + other.num_vertices

        return Tree(new_edges, new_leaves, new_num_vertices)

    def __repr__(self) -> str:
        return (
            f"Tree(leaves={len(self.leaves)}, "
            f"edges={len(self.edges)}, "
            f"vertices={self.num_vertices})"
        )


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _cpp_dm_from_protocol(dm: DistanceProtocol) -> _fastphylo.DistMatrix:
    """Build a C++ DistMatrix from any DistanceProtocol object."""
    n = len(dm)
    names = dm.names()
    cpp = _fastphylo.DistMatrix(n)
    for i, name in enumerate(names):
        cpp.set_name(i, name)
    for i in range(n):
        for j in range(n):
            cpp.set(i, j, dm[i, j])
    return cpp


def _raw_cpp_dm(dm: DistanceProtocol) -> _fastphylo.DistMatrix:
    """Return the underlying C++ DistMatrix, or build one from the protocol."""
    if hasattr(dm, "_cpp"):
        return dm._cpp  # type: ignore[attr-defined]
    return _cpp_dm_from_protocol(dm)


def _build_tree(
    cpp_fn,
    dm: DistanceProtocol,
) -> Tree:
    cpp_dm = _raw_cpp_dm(dm)
    edges_raw, leaf_names = cpp_fn(cpp_dm)

    n = len(leaf_names)
    leaves = {i: name for i, name in enumerate(leaf_names)}

    # num_vertices = max vertex ID + 1 (handles both small and large trees)
    all_ids = {u for u, v, _w in edges_raw} | {v for _u, v, _w in edges_raw}
    num_vertices = (max(all_ids) + 1) if all_ids else n

    return Tree(list(edges_raw), leaves, num_vertices)


# ---------------------------------------------------------------------------
# Public tree functions
# ---------------------------------------------------------------------------

def nj(dm: DistanceProtocol) -> Tree:
    """Reconstruct a Neighbour-Joining tree. Branch lengths are not computed (-1.0)."""
    return _build_tree(_fastphylo.nj_tree, dm)


def fnj(dm: DistanceProtocol) -> Tree:
    """Reconstruct a Fast Neighbour-Joining tree. Branch lengths are not computed (-1.0)."""
    return _build_tree(_fastphylo.fnj_tree, dm)


def bionj(dm: DistanceProtocol) -> Tree:
    """Reconstruct a BioNJ tree with real branch lengths."""
    return _build_tree(_fastphylo.bionj_tree, dm)


def fit_branch_lengths(tree: Tree, dm: DistanceProtocol) -> Tree:
    """Estimate branch lengths by L1-minimisation against a distance matrix.

    Returns a new Tree with the same topology but branch lengths chosen to
    minimise sum |path_dist(i,j) - dm[i,j]| over all leaf pairs, subject to
    non-negative branch lengths (requires scipy).
    """
    try:
        import numpy as np
        from scipy.optimize import linprog
    except ImportError as exc:
        raise ImportError(
            "fit_branch_lengths requires scipy: pip install scipy"
        ) from exc

    dm_names = dm.names()
    n = len(dm_names)
    name_to_vertex = {name: vid for vid, name in tree.leaves.items()}

    missing = [name for name in dm_names if name not in name_to_vertex]
    if missing:
        raise ValueError(f"Taxa in dm not found in tree: {missing}")

    # Directed edges → undirected index
    b = len(tree.edges)
    edge_idx: dict[tuple[int, int], int] = {}
    for k, (u, v, _) in enumerate(tree.edges):
        edge_idx[(u, v)] = k
        edge_idx[(v, u)] = k

    adj: dict[int, list[int]] = {}
    for u, v, _ in tree.edges:
        adj.setdefault(u, []).append(v)
        adj.setdefault(v, []).append(u)

    def _path_edges(src: int, dst: int) -> list[int]:
        from collections import deque
        parent_edge: dict[int, tuple[int, int]] = {}
        queue = deque([src])
        visited = {src}
        while queue:
            node = queue.popleft()
            if node == dst:
                path: list[int] = []
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
        return []

    # Path matrix P (m × b) and observed distances
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    P = np.zeros((m, b))
    d_obs = np.zeros(m)
    for k, (i, j) in enumerate(pairs):
        for eidx in _path_edges(name_to_vertex[dm_names[i]],
                                name_to_vertex[dm_names[j]]):
            P[k, eidx] = 1.0
        d_obs[k] = dm[i, j]

    # LP: minimise sum(s)  s.t.  P@x - s <= d,  -P@x - s <= -d,  x,s >= 0
    I_m = np.eye(m)
    result = linprog(
        c=np.concatenate([np.zeros(b), np.ones(m)]),
        A_ub=np.block([[P, -I_m], [-P, -I_m]]),
        b_ub=np.concatenate([d_obs, -d_obs]),
        bounds=[(0, None)] * (b + m),
        method="highs",
    )
    if not result.success:
        raise RuntimeError(f"Branch length LP did not converge: {result.message}")

    x_opt = result.x[:b]
    new_edges = [(u, v, float(x_opt[edge_idx[(u, v)]])) for u, v, _ in tree.edges]
    return Tree(new_edges, dict(tree.leaves), tree.num_vertices)
