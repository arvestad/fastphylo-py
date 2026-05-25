"""
Estimate pairwise distances from a FASTA alignment using fastphylo.

Usage: estimate_distances.py <alignment.fa> <output.tsv> <model>
                             [--method ml|expected] [--n-points N]

Output TSV columns: taxon_i, taxon_j, estimated_distance
(upper triangle only, same order as true_distances.tsv)
"""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fa_path")
    ap.add_argument("out_tsv")
    ap.add_argument("model")
    ap.add_argument("--method",   default="ml", choices=["ml", "expected", "expected_noprior"])
    ap.add_argument("--n-points", type=int, default=15)
    args = ap.parse_args()

    import fastphylo
    aln = fastphylo.read_fasta(args.fa_path)
    dm  = fastphylo.distance_matrix(aln, model=args.model,
                                    method=args.method,
                                    n_points=args.n_points)

    names = dm.names()
    n     = len(names)

    with open(args.out_tsv, "w") as out:
        out.write("taxon_i\ttaxon_j\testimated_distance\n")
        for i in range(n):
            for j in range(i + 1, n):
                out.write(f"{names[i]}\t{names[j]}\t{dm[i, j]}\n")


if __name__ == "__main__":
    main()
