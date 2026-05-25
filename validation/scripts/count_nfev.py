"""
Count likelihood function evaluations for the ML distance estimator.

Reproduces the bounded Brent minimizer (xtol=0.02, maxiter=20) used by the
C++ backend and records scipy's nfev counter for every sequence pair.

Usage: count_nfev.py <alignment.fa> <model> <output.tsv>
"""
import sys
import argparse


def _neg_log_lik(d, N, model):
    import numpy as np
    P = model.get_replacement_probs(d)
    P_safe = np.where(P > 0, P, 1e-300)
    return -float((N * _np_log(P_safe)).sum())


def _np_log(x):
    import numpy as np
    return np.log(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fa_path")
    ap.add_argument("model")
    ap.add_argument("out_tsv")
    args = ap.parse_args()

    import math
    import numpy as np
    from scipy.optimize import minimize_scalar
    import fastphylo
    from fastphylo.protein import (
        RateMatrix, _aa_count_matrix, kimura_protein_distance,
        _DELTA, _MAX_DIST,
    )

    aln   = fastphylo.read_fasta(args.fa_path)
    model = RateMatrix.instantiate(args.model)
    names = [s.accession for s in aln]
    seqs  = [s.data      for s in aln]
    n     = len(names)

    rows = []
    for i in range(n):
        for j in range(i + 1, n):
            N = _aa_count_matrix(seqs[i], seqs[j])
            d_k = kimura_protein_distance(seqs[i], seqs[j])
            if not math.isfinite(d_k) or d_k <= 0:
                d_k = 0.5

            def obj(d):
                P = model.get_replacement_probs(d)
                P_safe = np.where(P > 0, P, 1e-300)
                return -float(np.sum(N * np.log(P_safe)))

            result = minimize_scalar(
                obj,
                bounds=(_DELTA, _MAX_DIST),
                method="bounded",
                options={"xatol": 0.02, "maxiter": 20},
            )
            rows.append((names[i], names[j], result.nfev, result.x, d_k))

    with open(args.out_tsv, "w") as f:
        f.write("taxon_i\ttaxon_j\tnfev\td_ml\td_kimura\n")
        for ti, tj, nfev, d_ml, d_k in rows:
            f.write(f"{ti}\t{tj}\t{nfev}\t{d_ml:.6f}\t{d_k:.6f}\n")


if __name__ == "__main__":
    main()
