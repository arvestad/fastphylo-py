"""fastphylo — molecular sequence analysis."""

from .sequences import Alignment, Sequence, SequenceCollection, SeqType
from .io import read, read_fasta, read_stockholm, read_phylip, write_fasta, write_phylip
from .distances import DistanceMatrix, distance_matrix
from .trees import Tree, DistanceProtocol, nj, fnj, bionj, fit_branch_lengths
from .alignment import align

__version__ = "0.1.0"

__all__ = [
    "Alignment",
    "Sequence",
    "SequenceCollection",
    "SeqType",
    "read",
    "read_fasta",
    "read_stockholm",
    "read_phylip",
    "write_fasta",
    "write_phylip",
    "DistanceMatrix",
    "distance_matrix",
    "Tree",
    "DistanceProtocol",
    "nj",
    "fnj",
    "bionj",
    "fit_branch_lengths",
    "align",
]
