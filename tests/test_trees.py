import pytest
from conftest import DATA

import fastphylo
from fastphylo import read_stockholm, distance_matrix, nj, fnj, bionj, fit_branch_lengths, Tree, DistanceMatrix


@pytest.fixture
def dm():
    aln = read_stockholm(DATA / "dna_aligned.sto")
    return distance_matrix(aln, model="k2p")


def test_nj_produces_tree(dm):
    t = nj(dm)
    assert isinstance(t, Tree)
    assert len(t.leaves) == 4
    assert len(t.edges) > 0


def test_fnj_produces_tree(dm):
    t = fnj(dm)
    assert isinstance(t, Tree)
    assert len(t.leaves) == 4


def test_bionj_produces_tree(dm):
    t = bionj(dm)
    assert isinstance(t, Tree)
    assert len(t.leaves) == 4


def test_tree_leaf_count(dm):
    t = fnj(dm)
    # 4 taxa → 4 leaves with IDs 0-3
    assert len(t.leaves) == 4
    assert set(t.leaves.keys()) == {0, 1, 2, 3}


def test_tree_edge_count(dm):
    # Unrooted binary tree on N leaves has 2N-3 edges
    t = fnj(dm)
    n = len(t.leaves)
    assert len(t.edges) == 2 * n - 3


def test_tree_leaf_names_match_input(dm):
    t = fnj(dm)
    input_names = set(dm.names())
    assert set(t.leaves.values()) == input_names


def test_tree_to_newick_fnj(dm):
    t = fnj(dm)
    nwk = t.to_newick()
    # Valid Newick: ends with semicolon
    assert nwk.endswith(";")
    # All leaf names present
    for name in dm.names():
        assert name in nwk
    # FNJ has no branch lengths → no colons
    assert ":" not in nwk


def test_tree_to_newick_bionj(dm):
    t = bionj(dm)
    nwk = t.to_newick()
    assert nwk.endswith(";")
    for name in dm.names():
        assert name in nwk
    # BioNJ includes branch lengths → colons present
    assert ":" in nwk


def test_tree_to_newick():
    # Hand-build a known tree and verify the Newick string
    # Two-leaf star: leaf 0 and leaf 1 both connected to internal node 2
    # edges as (child, parent, w): (0, 2, 0.1), (1, 2, 0.2)
    # root = 2 (never appears as child)
    t = Tree(
        edges=[(0, 2, 0.1), (1, 2, 0.2)],
        leaves={0: "A", 1: "B"},
        num_vertices=3,
    )
    nwk = t.to_newick()
    assert nwk.endswith(";")
    assert "A" in nwk and "B" in nwk
    assert ":" in nwk  # real branch lengths → included


def test_tree_merge(dm):
    t1 = fnj(dm)
    t2 = fnj(dm)
    merged = t1.merge(t2)

    # Merged tree has combined leaf and edge counts
    assert len(merged.leaves) == len(t1.leaves) + len(t2.leaves)
    assert len(merged.edges) == len(t1.edges) + len(t2.edges)
    assert merged.num_vertices == t1.num_vertices + t2.num_vertices

    # t1 vertex IDs are unchanged
    for vid, name in t1.leaves.items():
        assert merged.leaves[vid] == name

    # t2 vertex IDs are offset by t1.num_vertices
    offset = t1.num_vertices
    for vid, name in t2.leaves.items():
        assert merged.leaves[vid + offset] == name

    # No vertex ID collisions
    t1_ids = {u for u, v, _ in t1.edges} | {v for u, v, _ in t1.edges}
    t2_ids_shifted = {u + offset for u, v, _ in t2.edges} | {v + offset for u, v, _ in t2.edges}
    assert t1_ids.isdisjoint(t2_ids_shifted)


def test_fit_branch_lengths_basic(dm):
    scipy = pytest.importorskip("scipy")
    t = fnj(dm)
    # All branch lengths are -1 (not computed by FNJ)
    assert all(w == -1.0 for _, _, w in t.edges)

    t_fit = fit_branch_lengths(t, dm)

    # Same topology
    assert len(t_fit.edges) == len(t.edges)
    assert t_fit.leaves == t.leaves
    # All branch lengths are non-negative
    assert all(w >= 0.0 for _, _, w in t_fit.edges)
    # Newick now includes branch lengths
    assert ":" in t_fit.to_newick()


def test_fit_branch_lengths_path_distances_close(dm):
    """Path distances in the fitted tree should be close to dm values."""
    scipy = pytest.importorskip("scipy")
    import numpy as np
    from collections import deque

    t = fit_branch_lengths(fnj(dm), dm)
    n = len(dm)

    # Build adjacency with weights
    adj: dict[int, list[tuple[int, float]]] = {}
    for u, v, w in t.edges:
        adj.setdefault(u, []).append((v, w))
        adj.setdefault(v, []).append((u, w))

    name_to_vid = {name: vid for vid, name in t.leaves.items()}

    def path_dist(src, dst):
        queue = deque([(src, 0.0)])
        visited = {src}
        while queue:
            node, dist = queue.popleft()
            if node == dst:
                return dist
            for nb, w in adj.get(node, []):
                if nb not in visited:
                    visited.add(nb)
                    queue.append((nb, dist + w))
        return float("inf")

    names = dm.names()
    max_err = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            pd = path_dist(name_to_vid[names[i]], name_to_vid[names[j]])
            max_err = max(max_err, abs(pd - dm[i, j]))

    # The LP minimises total deviation; individual errors should be small
    assert max_err < 0.01


def test_fit_branch_lengths_exact_recovery():
    """On a perfectly tree-additive distance matrix, LP recovers exact branch lengths."""
    scipy = pytest.importorskip("scipy")

    # Star topology: leaves 0,1,2,3 connected to internal node 4.
    # True branch lengths: a=0.1, b=0.2, c=0.15, d=0.3
    # Path distances: d(0,1)=a+b, d(0,2)=a+c, d(0,3)=a+d, d(1,2)=b+c, d(1,3)=b+d, d(2,3)=c+d
    a, b, c, d = 0.1, 0.2, 0.15, 0.3
    names = ["A", "B", "C", "D"]
    dm_star = DistanceMatrix.zeros(names)
    dm_star["A", "B"] = dm_star["B", "A"] = a + b
    dm_star["A", "C"] = dm_star["C", "A"] = a + c
    dm_star["A", "D"] = dm_star["D", "A"] = a + d
    dm_star["B", "C"] = dm_star["C", "B"] = b + c
    dm_star["B", "D"] = dm_star["D", "B"] = b + d
    dm_star["C", "D"] = dm_star["D", "C"] = c + d

    # Build the star tree manually (all leaves point to internal node 4)
    star = Tree(
        edges=[(0, 4, -1.0), (1, 4, -1.0), (2, 4, -1.0), (3, 4, -1.0)],
        leaves={0: "A", 1: "B", 2: "C", 3: "D"},
        num_vertices=5,
    )
    t_fit = fit_branch_lengths(star, dm_star)

    weight = {vid: w for vid, _, w in t_fit.edges}
    assert weight[0] == pytest.approx(a, abs=1e-6)
    assert weight[1] == pytest.approx(b, abs=1e-6)
    assert weight[2] == pytest.approx(c, abs=1e-6)
    assert weight[3] == pytest.approx(d, abs=1e-6)


def test_fit_branch_lengths_unknown_taxon(dm):
    scipy = pytest.importorskip("scipy")
    t = fnj(dm)
    dm2 = DistanceMatrix.zeros(["seq1", "seq2", "no_such_seq"])
    with pytest.raises(ValueError, match="not found in tree"):
        fit_branch_lengths(t, dm2)


def test_distance_protocol_numpy(dm):
    numpy = pytest.importorskip("numpy")

    # Wrap DistanceMatrix in a plain object that satisfies DistanceProtocol
    # (no ._cpp attribute) — exercises the fallback build path
    class NumpyDM:
        def __init__(self, src: DistanceMatrix):
            self._arr = src.to_numpy()
            self._names = src.names()

        def __len__(self):
            return len(self._names)

        def __getitem__(self, key):
            i, j = key
            return float(self._arr[i, j])

        def names(self):
            return list(self._names)

    proto_dm = NumpyDM(dm)

    # Confirm no ._cpp shortcut is taken
    assert not hasattr(proto_dm, "_cpp")

    t = fnj(proto_dm)
    assert isinstance(t, Tree)
    assert set(t.leaves.values()) == set(dm.names())
