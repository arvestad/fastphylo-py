"""fastphylo — molecular sequence analysis."""

from importlib.metadata import version as _version

from .sequences import Alignment, Sequence, SequenceCollection, SeqType
from .io import read, read_fasta, read_stockholm, read_phylip, write_fasta, write_phylip
from .distances import DistanceMatrix, distance_matrix
from .trees import Tree, DistanceProtocol, nj, fnj, bionj, fit_branch_lengths
from .alignment import align

# Read from the installed package's own metadata (pyproject.toml's
# version, at build time) instead of a second, hand-maintained literal -
# the two had already drifted (this used to say "1.0.0" while
# pyproject.toml had moved on to "1.1.0") before anyone noticed.
__version__ = _version("fastphylo")

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
