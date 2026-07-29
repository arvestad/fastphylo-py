"""
Benchmark: how well NJ, FNJ, and BioNJ reconstruct a known random tree
from a noisy pairwise-distance matrix.

For each number of leaves n and noise level sigma, random unrooted binary
trees are generated with i.i.d. Uniform(0.02, 0.1) edge lengths. The true
pairwise (patristic) distances are computed, then independently perturbed
with N(0, sigma) noise (clipped at 0 — real distance estimators cannot
return negative distances) before being handed to each reconstruction
method.

Metrics (see tree_metrics.py):
  rf         — normalised Robinson-Foulds distance to the true topology,
               in [0, 1]. Reported for all three methods.
  path_rmse  — RMSE between the reconstructed tree's leaf-pairwise path
               distances and the TRUE (noise-free) distances. Only
               meaningful for NJ and BioNJ, which produce real branch
               lengths; NaN for FNJ.

Usage:
  run_benchmark.py <output.tsv> [--replicates 30] [--seed 0]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))

from random_tree import random_tree
from tree_metrics import rf_distance, path_distances

N_VALUES     = [10, 20, 40, 80, 160, 320]
SIGMA_VALUES = [0.01, 0.025, 0.05, 0.075]
METHODS      = ["nj", "fnj", "bionj"]


def _noisy_matrix(true_dist, names, sigma, rng):
    import fastphylo
    n = len(names)
    dm = fastphylo.DistanceMatrix.zeros(names)
    for i in range(n):
        for j in range(i + 1, n):
            ni, nj = names[i], names[j]
            key = (ni, nj) if ni < nj else (nj, ni)
            d = true_dist[key]
            noisy = max(d + rng.normal(0.0, sigma), 0.0)
            dm[i, j] = noisy
            dm[j, i] = noisy
    return dm


def _path_rmse(recon_edges, recon_leaves, true_dist):
    import numpy as np
    recon_dist = path_distances(recon_edges, recon_leaves)
    errs = [recon_dist[k] - true_dist[k] for k in true_dist]
    return float(np.sqrt(np.mean(np.square(errs))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("output")
    ap.add_argument("--replicates", type=int, default=30)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import numpy as np
    import fastphylo

    rng = np.random.default_rng(args.seed)

    rows = []
    for n in N_VALUES:
        for rep in range(args.replicates):
            true_edges, true_leaves = random_tree(n, rng)
            names = [true_leaves[i] for i in range(n)]
            true_dist = path_distances(true_edges, true_leaves)

            for sigma in SIGMA_VALUES:
                dm = _noisy_matrix(true_dist, names, sigma, rng)

                for method in METHODS:
                    fn = getattr(fastphylo, method)
                    t = fn(dm)

                    rf = rf_distance(true_edges, true_leaves, t.edges, t.leaves)
                    rmse = (float("nan") if method == "fnj"
                            else _path_rmse(t.edges, t.leaves, true_dist))

                    rows.append(dict(
                        n=n, replicate=rep, sigma=sigma, method=method,
                        rf=rf, path_rmse=rmse,
                    ))

        print(f"n={n} done", file=sys.stderr)

    import pandas as pd
    df = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output}")


if __name__ == "__main__":
    main()
