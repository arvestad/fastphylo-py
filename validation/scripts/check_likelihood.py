"""
For every sequence pair in one replicate, compute the log-likelihood at the
estimated distance and at the true (tree-derived) distance.

A "violation" is when LL(t_true) > LL(t_estimated) + TOL — which would mean
the optimizer failed to find the maximum.

Usage:
  check_likelihood.py <alignment.fa> <estimated.tsv> <true.tsv> <model> <output.tsv>

Output TSV columns:
  taxon_i  taxon_j  t_true  t_est  ll_true  ll_est  violation
"""
import sys
import math
import numpy as np

TOL = 1e-6   # numerical tolerance for declaring a violation

# ---------------------------------------------------------------------------
# DNA helpers
# ---------------------------------------------------------------------------
_DNA_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
             'a': 0, 'c': 1, 'g': 2, 't': 3}
_TRANSITIONS = {(0, 2), (2, 0), (1, 3), (3, 1)}   # A↔G, C↔T


def _dna_count_matrix(s1, s2):
    N = np.zeros((4, 4))
    for a, b in zip(s1, s2):
        i = _DNA_IDX.get(a, -1)
        j = _DNA_IDX.get(b, -1)
        if i >= 0 and j >= 0:
            N[i, j] += 1.0
    return N


def _ll(N, P):
    """log-likelihood given count matrix N and probability matrix P."""
    ll = 0.0
    mask = N > 0
    P_safe = np.where(P > 0, P, 1e-300)
    ll = float(np.sum(N[mask] * np.log(P_safe[mask])))
    return ll


def _jc_matrix(t):
    if t <= 0:
        return np.eye(4)
    e = math.exp(-4 * t / 3)
    P = np.full((4, 4), 0.25 - 0.25 * e)
    np.fill_diagonal(P, 0.25 + 0.75 * e)
    return P


def _k2p_matrix(t, kappa):
    if t <= 0:
        return np.eye(4)
    beta  = 1.0 / (kappa + 2)
    alpha = kappa * beta
    e1 = math.exp(-4 * beta * t)
    e2 = math.exp(-2 * (alpha + beta) * t)
    p_same  = 0.25 + 0.25 * e1 + 0.5 * e2
    p_trans = 0.25 + 0.25 * e1 - 0.5 * e2
    p_transv = 0.25 - 0.25 * e1
    P = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            if i == j:
                P[i, j] = p_same
            elif (i, j) in _TRANSITIONS:
                P[i, j] = p_trans
            else:
                P[i, j] = p_transv
    return P


def _tn93_matrix(t, kappa_r, kappa_y, freqs):
    """TN93 via matrix exponentiation (scipy.linalg.expm)."""
    from scipy.linalg import expm
    piA, piC, piG, piT = freqs
    beta = 1.0  # free scale; we fix overall rate via normalization below
    # Rate matrix Q (A=0, C=1, G=2, T=3)
    Q = np.array([
        [0,       piC,          piG * kappa_r, piT       ],
        [piA,     0,            piG,           piT * kappa_y],
        [piA * kappa_r, piC,   0,             piT       ],
        [piA,     piC * kappa_y, piG,          0         ],
    ])
    np.fill_diagonal(Q, 0)
    Q -= np.diag(Q.sum(axis=1))
    # Normalize: mean rate = -sum_i pi_i Q_ii = 1
    pi = np.array(freqs)
    mean_rate = -float(pi @ np.diag(Q))
    if mean_rate > 0:
        Q /= mean_rate
    return expm(Q * t)


def _profile_ll_k2p(N, t):
    """Profile log-likelihood for K2P at fixed t, optimised over kappa."""
    from scipy.optimize import minimize_scalar
    def neg_ll(kappa):
        return -_ll(N, _k2p_matrix(t, kappa))
    res = minimize_scalar(neg_ll, bounds=(0.01, 200.0), method='bounded')
    return -res.fun


def _profile_ll_tn93(N, t, freqs):
    """Profile log-likelihood for TN93 at fixed t, optimised over kappa_R, kappa_Y."""
    from scipy.optimize import minimize
    def neg_ll(x):
        kr, ky = x
        return -_ll(N, _tn93_matrix(t, kr, ky, freqs))
    res = minimize(neg_ll, x0=[2.0, 2.0], bounds=[(0.01, 100), (0.01, 100)],
                   method='L-BFGS-B')
    return -res.fun


def _dna_profile_ll(model, N, t, freqs):
    if model == "jc":
        return _ll(N, _jc_matrix(t))
    elif model == "k2p":
        return _profile_ll_k2p(N, t)
    elif model == "tn93":
        return _profile_ll_tn93(N, t, freqs)
    else:
        raise ValueError(f"Unknown DNA model: {model}")


# ---------------------------------------------------------------------------
# Protein helpers
# ---------------------------------------------------------------------------
def _protein_count_matrix(s1, s2, aa_idx):
    N = np.zeros((20, 20))
    for a, b in zip(s1, s2):
        i = aa_idx.get(a.upper(), -1)
        j = aa_idx.get(b.upper(), -1)
        if i >= 0 and j >= 0:
            N[i, j] += 1.0
    return N


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <alignment.fa> <estimated.tsv> <true.tsv> <model> <output.tsv>")
        sys.exit(1)

    fa_path, est_path, true_path, model, out_path = sys.argv[1:]

    import fastphylo
    import pandas as pd

    aln      = fastphylo.read_fasta(fa_path)
    seq_dict = {s.accession: s.data for s in aln}

    est_df  = pd.read_csv(est_path,  sep="\t")
    true_df = pd.read_csv(true_path, sep="\t")
    merged  = true_df.merge(est_df, on=["taxon_i", "taxon_j"])

    # Detect dtype from model name
    dna_models = {"jc", "k2p", "tn93"}
    is_dna = model.lower() in dna_models

    if is_dna:
        # Empirical base frequencies (for TN93)
        counts = np.zeros(4)
        for s in seq_dict.values():
            for c in s:
                i = _DNA_IDX.get(c, -1)
                if i >= 0:
                    counts[i] += 1
        total = counts.sum()
        freqs = (counts / total).tolist() if total > 0 else [0.25] * 4
    else:
        from fastphylo.protein import RateMatrix, AA_ORDER
        rate_model = RateMatrix.instantiate(model)
        aa_idx = {aa: i for i, aa in enumerate(AA_ORDER)}
        aa_idx.update({aa.lower(): i for i, aa in enumerate(AA_ORDER)})

    rows = []
    for _, row in merged.iterrows():
        ti, tj   = row["taxon_i"], row["taxon_j"]
        t_true   = float(row["true_distance"])
        t_est    = float(row["estimated_distance"])

        if not math.isfinite(t_est) or not math.isfinite(t_true):
            continue

        s1 = seq_dict.get(ti, "")
        s2 = seq_dict.get(tj, "")
        if not s1 or not s2:
            continue

        if is_dna:
            N  = _dna_count_matrix(s1, s2)
            ll_est  = _dna_profile_ll(model, N, t_est,  freqs)
            ll_true = _dna_profile_ll(model, N, t_true, freqs)
        else:
            N  = _protein_count_matrix(s1, s2, aa_idx)
            P_est  = rate_model.get_replacement_probs(t_est)
            P_true = rate_model.get_replacement_probs(t_true)
            ll_est  = _ll(N, np.array(P_est))
            ll_true = _ll(N, np.array(P_true))

        violation = int(ll_true > ll_est + TOL)
        rows.append((ti, tj, t_true, t_est, ll_true, ll_est, violation))

    df_out = pd.DataFrame(rows, columns=[
        "taxon_i", "taxon_j", "t_true", "t_est", "ll_true", "ll_est", "violation"
    ])
    df_out.to_csv(out_path, sep="\t", index=False, float_format="%.8f")


if __name__ == "__main__":
    main()
