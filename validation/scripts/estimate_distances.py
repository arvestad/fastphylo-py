"""
Estimate pairwise distances from a FASTA alignment using fastphylo.

Usage: estimate_distances.py <alignment.fa> <output.tsv> <model>

Output TSV columns: taxon_i, taxon_j, estimated_distance
(upper triangle only, same order as true_distances.tsv)
"""
import math
import sys


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <alignment.fa> <output.tsv> <model>")
        sys.exit(1)

    fa_path, out_tsv, model = sys.argv[1], sys.argv[2], sys.argv[3]

    import fastphylo
    aln = fastphylo.read_fasta(fa_path)
    dm  = fastphylo.distance_matrix(aln, model=model)

    names = dm.names()
    n     = len(names)

    with open(out_tsv, "w") as out:
        out.write("taxon_i\ttaxon_j\testimated_distance\n")
        for i in range(n):
            for j in range(i + 1, n):
                d = dm[i, j]
                out.write(f"{names[i]}\t{names[j]}\t{d}\n")


if __name__ == "__main__":
    main()
