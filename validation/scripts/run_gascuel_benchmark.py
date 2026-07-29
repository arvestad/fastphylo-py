"""
Recreate Gascuel (1997)'s BIONJ simulation study: for each of the six
12-taxon model trees (A-F), four substitution-rate conditions, and two
sequence lengths, simulate replicate alignments under the Kimura (1980)
two-parameter model (ts/tv = 2, via kappa = 4) using IQ-TREE2's AliSim,
estimate K2P pairwise distances, reconstruct NJ and BioNJ trees, and
compare:
  - RF distance to the true tree topology (both algorithms)
  - ME criterion (OLS-refit tree length) for the true/NJ/BioNJ topologies,
    fit to the (noisy) estimated distance matrix (see ols_branch_lengths.py)

Rate conditions (interior-branch unit values, Gascuel 1997 p.690):
  low:      a=0.0055  b=0.004   (max pairwise divergence ~0.1)
  middle:   a=0.0275  b=0.02    (~0.5)
  high:     a=0.055   b=0.04    (~1.0)
  high_low: half the sites at "low" rate, half at "high" rate (~0.55)

Trees A, B use the constant-rate unit `a`; trees C-F use the
variable-rate unit `b` (see gascuel_trees.py). The ME criterion only
needs the true TOPOLOGY (branch lengths are always OLS-refit against the
estimated distance matrix, per Gascuel's own methodology), so a single
condition-independent reference tree per model is used for both RF and ME
comparisons -- only the *simulation* step needs real branch lengths.

Usage:
  run_gascuel_benchmark.py <output.tsv> [--replicates 50] [--seed 0]
      [--trees A,B,C,D,E,F] [--conditions low,middle,high,high_low]
      [--lengths 300,600] [--workdir DIR] [--iqtree iqtree2]
"""
import argparse
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))

from gascuel_trees import MODEL_TREES
from tree_metrics import rf_distance
from ols_branch_lengths import ols_tree_length

RATE_UNITS = {
    "low":    {"a": 0.0055, "b": 0.004},
    "middle": {"a": 0.0275, "b": 0.02},
    "high":   {"a": 0.055,  "b": 0.04},
}
TREE_UNIT_KIND = {"A": "a", "B": "a", "C": "b", "D": "b", "E": "b", "F": "b"}
ALL_CONDITIONS = ["low", "middle", "high", "high_low"]
KAPPA = 4  # K80 kappa calibrated to ts/tv = 2 (kappa = 2*ts/tv for equal base freqs)


# ---------------------------------------------------------------------------
# Tree -> Newick, sequence simulation
# ---------------------------------------------------------------------------

def _to_newick(edges, leaves) -> str:
    adj = defaultdict(list)
    for u, v, w in edges:
        adj[u].append((v, w))
        adj[v].append((u, w))
    root = next(n for n in adj if len(adj[n]) > 1)  # avoid leaf-root double-wrap

    def rec(node, parent):
        kids = [(nb, w) for nb, w in adj[node] if nb != parent]
        if not kids:
            return leaves.get(node, str(node))
        parts = [f"{rec(c, node)}:{w:.8f}" for c, w in kids]
        return "(" + ",".join(parts) + ")"

    return rec(root, None) + ";"


def _run_alisim(newick: str, length: int, workdir: str, prefix: str,
                 iqtree_bin: str, seed: int) -> str:
    tree_path = os.path.join(workdir, f"{prefix}.nwk")
    with open(tree_path, "w") as fh:
        fh.write(newick)
    out_prefix = os.path.join(workdir, prefix)
    cmd = [
        iqtree_bin, "--alisim", out_prefix,
        "-t", tree_path,
        "-m", f"K80{{{KAPPA}}}",
        "--seqtype", "DNA",
        "--length", str(length),
        "--out-format", "fasta",
        "--redo", "-quiet",
        "--seed", str(seed),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    fa_path = out_prefix + ".fa"
    if result.returncode != 0 or not os.path.exists(fa_path):
        raise RuntimeError(f"AliSim failed ({' '.join(cmd)}):\n{result.stderr}")
    return fa_path


def _read_fasta_dict(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            else:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def simulate_alignment(tree_name: str, condition: str, length: int, workdir: str,
                        rep_seed: int, iqtree_bin: str) -> dict[str, str]:
    """Return dict name -> sequence for one replicate."""
    kind = TREE_UNIT_KIND[tree_name]
    builder = MODEL_TREES[tree_name]

    if condition != "high_low":
        unit = RATE_UNITS[condition][kind]
        edges, leaves = builder(unit)
        newick = _to_newick(edges, leaves)
        fa = _run_alisim(newick, length, workdir, f"{tree_name}_{condition}_{rep_seed}",
                          iqtree_bin, rep_seed)
        return _read_fasta_dict(fa)

    # high_low: half the sites at "low" rate, half at "high" rate, same topology.
    half = length // 2
    parts = []
    for tag, half_len, rate, seed_off in [("lo", half, "low", 0),
                                           ("hi", length - half, "high", 1)]:
        unit = RATE_UNITS[rate][kind]
        edges, leaves = builder(unit)
        newick = _to_newick(edges, leaves)
        fa = _run_alisim(newick, half_len, workdir, f"{tree_name}_hl{tag}_{rep_seed}",
                          iqtree_bin, rep_seed * 2 + seed_off)
        parts.append(_read_fasta_dict(fa))

    names = parts[0].keys()
    return {n: parts[0][n] + parts[1][n] for n in names}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("output")
    ap.add_argument("--replicates", type=int, default=50)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--trees", default="A,B,C,D,E,F")
    ap.add_argument("--conditions", default=",".join(ALL_CONDITIONS))
    ap.add_argument("--lengths", default="300,600")
    ap.add_argument("--workdir", default=None)
    ap.add_argument("--iqtree", default="iqtree2")
    args = ap.parse_args()

    tree_names = args.trees.split(",")
    conditions = args.conditions.split(",")
    lengths = [int(x) for x in args.lengths.split(",")]

    import fastphylo

    workdir = args.workdir or tempfile.mkdtemp(prefix="gascuel_")
    os.makedirs(workdir, exist_ok=True)

    # Condition-independent reference topology: RF and ME both only need
    # the topology (ME always refits lengths via OLS), so one unit-length
    # reference per tree suffices for every condition.
    ref_topology = {t: MODEL_TREES[t](1.0) for t in tree_names}

    rows = []
    for tree_name in tree_names:
        true_edges, true_leaves = ref_topology[tree_name]
        for condition in conditions:
            for length in lengths:
                for rep in range(args.replicates):
                    rep_seed = (args.seed * 104729
                                + hash((tree_name, condition, length, rep)) % 1_000_000)

                    seqs = simulate_alignment(tree_name, condition, length, workdir,
                                               rep_seed, args.iqtree)

                    fa_path = os.path.join(workdir, f"{tree_name}_{condition}_{length}_{rep}.fa")
                    with open(fa_path, "w") as fh:
                        for name, seq in seqs.items():
                            fh.write(f">{name}\n{seq}\n")

                    aln = fastphylo.read_fasta(fa_path)
                    dm = fastphylo.distance_matrix(aln, model="k2p")
                    dm_names = dm.names()

                    me_true = ols_tree_length(true_edges, true_leaves, dm_names, dm)

                    row = dict(tree=tree_name, condition=condition, length=length,
                               replicate=rep, me_true=me_true)

                    for method in ("nj", "bionj"):
                        t = getattr(fastphylo, method)(dm)
                        rf = rf_distance(true_edges, true_leaves, t.edges, t.leaves)
                        me = ols_tree_length(t.edges, t.leaves, dm_names, dm)
                        row[f"rf_{method}"] = rf
                        row[f"me_{method}"] = me

                    rows.append(row)
                    os.remove(fa_path)

                print(f"{tree_name} {condition} L={length}: {args.replicates} reps done",
                      file=sys.stderr)

    import pandas as pd
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()
