"""
example1.py — read an alignment, compute distances, output a Newick tree.

Usage:
    python example1.py <alignment-file> [distance-method]

The alignment file can be FASTA, Stockholm, or Phylip format; the format is
detected automatically from the file extension and content.

Distance methods:
    DNA/RNA : hamming | jc | k2p (default) | tn93
    Protein : WAG (default) | LG | JTT | Dayhoff | ...

Examples:
    python example1.py sequences.sto
    python example1.py sequences.fasta k2p
    python example1.py proteins.phy WAG
"""

import sys
import fastphylo

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    path = sys.argv[1]
    method = sys.argv[2] if len(sys.argv) > 2 else None

    # Read the alignment — format auto-detected from extension and content.
    # Returns an Alignment (all sequences same length, gaps allowed).
    alignment = fastphylo.read(path)

    # Compute pairwise distance matrix.
    # If method is not given, fastphylo picks a default based on sequence type:
    #   DNA/RNA → k2p,  Protein → WAG
    dm = fastphylo.distance_matrix(alignment, method=method)

    # Reconstruct the tree with Fast Neighbor Joining.
    tree = fastphylo.fnj(dm)

    print(tree.to_newick())


if __name__ == "__main__":
    main()
