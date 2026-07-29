"""
The six 12-taxon model trees (A-F) from Gascuel (1997), "BIONJ: An Improved
Version of the NJ Algorithm Based on a Simple Model of Sequence Data",
Mol. Biol. Evol. 14(7):685-695, Figure 3.

Each tree is two copies of the same six-taxon tree (Kumar 1996 / Saitou &
Imanishi 1989), joined by one internal edge. Interior branches are one unit
long (`a` for the constant-rate trees A/B, `b` for the variable-rate trees
C-F); pendant (external) branch lengths are given as multiples of that unit.

Trees A and B are ultrametric (molecular clock): all 12 leaves are
equidistant from the tree's center. The labeled cherry pendant length (4a
for A, 6a for B) combined with unit-length interior branches determines
every other pendant length by the ultrametric constraint (see
_ultrametric_caterpillar / _ultrametric_balanced below).

Trees C and D have explicit, non-ultrametric branch lengths (varying
substitution rates among lineages). Trees E and F are "identical to C and
D... except that the short external branches with length b have been
replaced by longer branches having length 3b" (Gascuel 1997, p. 689) --
read here as the cherry's default-length member (t2) becoming 3b. This E/F
substitution is an extrapolation, not confirmed against the figure directly
-- flagged here for correction if wrong.

All six trees are built as two EXACT copies of their six-taxon half,
joined by one interior (unit-length) edge, per "each consists of two
copies of the same six-taxon tree."
"""
from __future__ import annotations


# ---------------------------------------------------------------------------
# Six-taxon half-tree builders.
# Each returns (edges, root) where edges is a list of (u, v, length) with
# u, v either leaf ids (from `ids`) or string labels for internal nodes,
# and `root` is the half's attachment point (joined to the other half).
# ---------------------------------------------------------------------------

def _ultrametric_caterpillar(unit: float, cherry_len: float, ids):
    """Fully pectinate 6-taxon tree, ultrametric, unit-length interior edges.

    Topology: (((((t1,t2),t3),t4),t5),t6)
    t1 and t2 (the shallowest cherry) have pendant length `cherry_len`;
    each subsequent sibling's pendant length increases by `unit` (the only
    way to keep the tree ultrametric with unit-length interior branches).
    """
    t1, t2, t3, t4, t5, t6 = ids
    n12    = f"n12_{t1}_{t2}"
    n123   = f"n123_{t1}_{t2}"
    n1234  = f"n1234_{t1}_{t2}"
    n12345 = f"n12345_{t1}_{t2}"
    root   = f"root_{t1}_{t2}"

    edges = [
        (t1, n12, cherry_len),
        (t2, n12, cherry_len),
        (n12, n123, unit),
        (t3, n123, cherry_len + unit),
        (n123, n1234, unit),
        (t4, n1234, cherry_len + 2 * unit),
        (n1234, n12345, unit),
        (t5, n12345, cherry_len + 3 * unit),
        (n12345, root, unit),
        (t6, root, cherry_len + 4 * unit),
    ]
    return edges, root


def _ultrametric_balanced(unit: float, cherry_len: float, ids):
    """Balanced 6-taxon tree (two mirrored 3-leaf pectinate triples),
    ultrametric, unit-length interior edges.

    Topology: (((t1,t2),t3), ((t4,t5),t6))
    Both cherries have pendant length `cherry_len`; both "third" leaves
    (t3, t6) have pendant length `cherry_len + unit` (ultrametric
    constraint, same reasoning as the caterpillar case).
    """
    t1, t2, t3, t4, t5, t6 = ids
    n12  = f"n12_{t1}_{t2}"
    tr1  = f"tr1_{t1}_{t2}"
    n45  = f"n45_{t1}_{t2}"
    tr2  = f"tr2_{t1}_{t2}"
    root = f"root_{t1}_{t2}"

    edges = [
        (t1, n12, cherry_len),
        (t2, n12, cherry_len),
        (n12, tr1, unit),
        (t3, tr1, cherry_len + unit),
        (t4, n45, cherry_len),
        (t5, n45, cherry_len),
        (n45, tr2, unit),
        (t6, tr2, cherry_len + unit),
        (tr1, root, unit),
        (tr2, root, unit),
    ]
    return edges, root


def _explicit_half(unit: float, lengths, shape: str, ids):
    """Non-ultrametric 6-taxon half with explicit pendant lengths (used for
    C, D, E, F). `lengths` gives (t1..t6) pendant lengths. All interior
    edges = unit.

    shape "C": (((t1,t2),t3),(t4,t5)),t6)   -- C/E topology
    shape "D": (t1,(t2,((t3,t4),(t5,t6))))  -- D/F topology
    """
    t1, t2, t3, t4, t5, t6 = ids
    L1, L2, L3, L4, L5, L6 = lengths

    if shape == "C":
        n12    = f"n12_{t1}"
        clade1 = f"clade1_{t1}"
        n45    = f"n45_{t1}"
        clade2 = f"clade2_{t1}"
        root   = f"root_{t1}"
        edges = [
            (t1, n12, L1),
            (t2, n12, L2),
            (n12, clade1, unit),
            (t3, clade1, L3),
            (t4, n45, L4),
            (t5, n45, L5),
            (n45, clade2, unit),
            (clade1, clade2, unit),
            (clade2, root, unit),
            (t6, root, L6),
        ]
        return edges, root

    if shape == "D":
        n34  = f"n34_{t1}"
        n56  = f"n56_{t1}"
        y    = f"y_{t1}"
        x    = f"x_{t1}"
        root = f"root_{t1}"
        edges = [
            (t3, n34, L3),
            (t4, n34, L4),
            (t5, n56, L5),
            (t6, n56, L6),
            (n34, y, unit),
            (n56, y, unit),
            (t2, x, L2),
            (y, x, unit),
            (t1, root, L1),
            (x, root, unit),
        ]
        return edges, root

    raise ValueError(f"unknown shape {shape!r}")


# ---------------------------------------------------------------------------
# Full 12-taxon trees: two exact copies of a half, joined by one unit edge
# ---------------------------------------------------------------------------

def _mirror(half_builder, unit: float, *args):
    """Build a 12-taxon tree from two exact copies of a 6-taxon half.

    Leaf ids 0-5 are the first copy, 6-11 the second (both use the SAME
    `half_builder(unit, *args, ids)` lengths, since Gascuel's trees are two
    copies of the same six-taxon tree). Internal-node string labels are
    generated from each half's own leaf ids, so the two copies' labels
    never collide without needing an extra suffixing step.
    """
    leaves = {i: f"t{i + 1}" for i in range(12)}

    edges1, root1 = half_builder(unit, *args, ids=list(range(0, 6)))
    edges2, root2 = half_builder(unit, *args, ids=list(range(6, 12)))
    edges = edges1 + edges2 + [(root1, root2, unit)]

    # Map string internal-node labels to contiguous integer ids >= 12.
    node_ids: dict = {}
    next_id = 12
    final_edges = []
    for u, v, w in edges:
        if not isinstance(u, int):
            node_ids.setdefault(u, next_id)
            if node_ids[u] == next_id:
                next_id += 1
            u = node_ids[u]
        if not isinstance(v, int):
            node_ids.setdefault(v, next_id)
            if node_ids[v] == next_id:
                next_id += 1
            v = node_ids[v]
        final_edges.append((u, v, float(w)))

    return final_edges, leaves


def tree_A(a: float):
    """Constant-rate, pectinate, ultrametric. Cherry pendant = 4a."""
    return _mirror(_ultrametric_caterpillar, a, 4 * a)


def tree_B(a: float):
    """Constant-rate, balanced, ultrametric. Cherry pendant = 6a."""
    return _mirror(_ultrametric_balanced, a, 6 * a)


def tree_C(b: float):
    """Highly-varying-rate. t1=9b,t2=1b,t3=1b,t4=1b,t5=9b,t6=8b."""
    lengths = (9 * b, 1 * b, 1 * b, 1 * b, 9 * b, 8 * b)
    return _mirror(_explicit_half, b, lengths, "C")


def tree_D(b: float):
    """Highly-varying-rate. t1=8b,t2=1b,t3=1b,t4=9b,t5=1b,t6=9b."""
    lengths = (8 * b, 1 * b, 1 * b, 9 * b, 1 * b, 9 * b)
    return _mirror(_explicit_half, b, lengths, "D")


def tree_E(b: float):
    """Intermediate: tree C with the cherry's default-length member (t2)
    elevated from 1b to 3b. EXTRAPOLATED from Gascuel's description, not
    read directly off Fig. 3 -- verify against the paper if precision
    matters."""
    lengths = (9 * b, 3 * b, 1 * b, 1 * b, 9 * b, 8 * b)
    return _mirror(_explicit_half, b, lengths, "C")


def tree_F(b: float):
    """Intermediate: tree D with t2 elevated from 1b to 3b. EXTRAPOLATED,
    see tree_E docstring."""
    lengths = (8 * b, 3 * b, 1 * b, 9 * b, 1 * b, 9 * b)
    return _mirror(_explicit_half, b, lengths, "D")


MODEL_TREES = {
    "A": tree_A, "B": tree_B, "C": tree_C,
    "D": tree_D, "E": tree_E, "F": tree_F,
}
