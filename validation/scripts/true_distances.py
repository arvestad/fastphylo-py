"""
Parse a Newick treefile and compute all pairwise path-length distances.

Usage: true_distances.py <treefile> <output.tsv>

Output TSV columns: taxon_i, taxon_j, true_distance
(upper triangle only, i < j by name sort)
"""
import sys
import re
from collections import defaultdict, deque


# ---------------------------------------------------------------------------
# Lightweight Newick parser — no third-party deps
# ---------------------------------------------------------------------------

class _Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name="", length=0.0):
        self.name     = name
        self.length   = length
        self.children = []
        self.parent   = None


def _parse_newick(text: str) -> _Node:
    """Parse a Newick string into a tree of _Node objects."""
    text = text.strip().rstrip(";")
    root = _Node()
    stack = [root]
    current = root
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c == "(":
            child = _Node()
            child.parent = current
            current.children.append(child)
            stack.append(child)
            current = child
            i += 1
        elif c == ",":
            current = stack[-1]   # back to parent
            child = _Node()
            child.parent = current
            current.children.append(child)
            stack.append(child)
            current = child
            i += 1
        elif c == ")":
            stack.pop()
            current = stack[-1]
            i += 1
        elif c == ":":
            # read branch length
            j = i + 1
            while j < n and text[j] not in "(),;":
                j += 1
            current.length = float(text[i+1:j])
            i = j
        else:
            # read label
            j = i
            while j < n and text[j] not in "():,;":
                j += 1
            current.name = text[i:j]
            i = j
    return root


def _collect_leaves(node: _Node):
    if not node.children:
        yield node
    for child in node.children:
        yield from _collect_leaves(child)


def _path_distance(a: _Node, b: _Node) -> float:
    """Sum of branch lengths on the unique path between two nodes."""
    # BFS upward from each node, collecting distance to each ancestor
    def ancestors(node):
        dist = {}
        d = 0.0
        cur = node
        while cur is not None:
            dist[id(cur)] = d
            d += cur.length
            cur = cur.parent
        return dist

    anc_a = ancestors(a)
    anc_b = ancestors(b)
    # LCA is the first ancestor of b found in anc_a
    cur = b
    while cur is not None:
        if id(cur) in anc_a:
            return anc_a[id(cur)] + anc_b[id(cur)]
        cur = cur.parent
    raise ValueError("Nodes do not share a common ancestor")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <treefile> <output.tsv>")
        sys.exit(1)

    treefile, out_tsv = sys.argv[1], sys.argv[2]

    with open(treefile) as fh:
        newick = fh.read().strip()

    root   = _parse_newick(newick)
    leaves = sorted(_collect_leaves(root), key=lambda n: n.name)

    with open(out_tsv, "w") as out:
        out.write("taxon_i\ttaxon_j\ttrue_distance\n")
        for i, leaf_i in enumerate(leaves):
            for leaf_j in leaves[i+1:]:
                d = _path_distance(leaf_i, leaf_j)
                out.write(f"{leaf_i.name}\t{leaf_j.name}\t{d:.10f}\n")


if __name__ == "__main__":
    main()
