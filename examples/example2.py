"""
example2.py — read unaligned sequences, align, compute distances, output a Newick tree.

Requires the optional alignment dependency:
    pip install fastphylo[align]

Usage:
    python example2.py <sequences-file> [distance-method]

The input file can be FASTA, Stockholm, or Phylip format; the format is
detected automatically. Sequences do not need to be the same length.

Distance methods:
    DNA/RNA : hamming | jc | k2p (default) | tn93
    Protein : WAG (default) | LG | JTT | Dayhoff | ...

Examples:
    python example2.py sequences.fasta
    python example2.py sequences.fasta tn93
    python example2.py proteins.fasta LG
"""

import sys
import fastphylo

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    method = sys.argv[2] if len(sys.argv) > 2 else None

    # Read unaligned sequences.
    # Returns a SequenceCollection (sequences may have different lengths).
    collection = fastphylo.read(path)

    # Align using FAMSA. The scoring matrix is pinned explicitly so results
    # are reproducible regardless of the installed pyfamsa version.
    alignment = fastphylo.align(collection, scoring_matrix="PFASUM43")

    # Compute pairwise distance matrix and reconstruct tree.
    dm = fastphylo.distance_matrix(alignment, method=method)
    tree = fastphylo.fnj(dm)

    print(tree.to_newick())


if __name__ == "__main__":
    main()
