"""Distance computation: DNA distance wrappers and the DistanceMatrix class."""

from __future__ import annotations

import numpy as np

from . import _fastphylo
from .sequences import Alignment, SeqType

_DNA_MODELS = frozenset(("hamming", "jc", "k2p", "tn93"))
_DEFAULT_DNA_MODEL = "k2p"
_DEFAULT_PROTEIN_MODEL = "WAG"


# ---------------------------------------------------------------------------
# DistanceMatrix  (Subplan G)
# ---------------------------------------------------------------------------

class DistanceMatrix:
    """Python wrapper around _fastphylo.DistMatrix."""

    def __init__(self, cpp: _fastphylo.DistMatrix) -> None:
        self._cpp = cpp

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def zeros(cls, names: list[str]) -> "DistanceMatrix":
        """Create an all-zeros distance matrix with the given taxon names."""
        n = len(names)
        cpp = _fastphylo.DistMatrix(n)
        for i, name in enumerate(names):
            cpp.set_name(i, name)
        return cls(cpp)

    # ------------------------------------------------------------------
    # Core accessors
    # ------------------------------------------------------------------

    def __len__(self) -> int:
        return self._cpp.size()

    def _resolve(self, k: int | str) -> int:
        """Resolve a taxon name or integer index to a row/column index."""
        if isinstance(k, str):
            names = self._cpp.names()
            try:
                return names.index(k)
            except ValueError:
                raise KeyError(k) from None
        return k

    def __getitem__(self, key: tuple) -> float:
        i, j = key
        return self._cpp.get(self._resolve(i), self._resolve(j))

    def __setitem__(self, key: tuple, val: float) -> None:
        i, j = key
        self._cpp.set(self._resolve(i), self._resolve(j), val)

    def names(self) -> list[str]:
        return self._cpp.names()

    def copy(self) -> "DistanceMatrix":
        """Return a deep copy of this distance matrix."""
        n = len(self)
        cpp = _fastphylo.DistMatrix(n)
        for i in range(n):
            cpp.set_name(i, self._cpp.name(i))
            for j in range(n):
                cpp.set(i, j, self._cpp.get(i, j))
        return DistanceMatrix(cpp)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def to_numpy(self) -> np.ndarray:
        n = len(self)
        arr = np.empty((n, n), dtype=float)
        for i in range(n):
            for j in range(n):
                arr[i, j] = self._cpp.get(i, j)
        return arr

    def to_phylip(self) -> str:
        """Return the distance matrix as a Phylip-format string."""
        n = len(self)
        lines = [f"{n:6}"]
        for i in range(n):
            name = self._cpp.name(i)[:10].ljust(10)
            dists = "  ".join(f"{self._cpp.get(i, j):.6f}" for j in range(n))
            lines.append(f"{name}  {dists}")
        return "\n".join(lines) + "\n"

    def __repr__(self) -> str:
        return f"DistanceMatrix({len(self)} taxa)"


# ---------------------------------------------------------------------------
# distance_matrix  (Subplan F)
# ---------------------------------------------------------------------------

def distance_matrix(
    alignment: Alignment,
    model: str | None = None,
    method: str = "ml",
    n_points: int = 15,
) -> DistanceMatrix:
    """Compute pairwise distances from an Alignment.

    For DNA/RNA: model ∈ {'hamming', 'jc', 'k2p', 'tn93'} (default: 'k2p').
    For Protein: model ∈ {'WAG', 'LG', 'JTT', ...} (default: 'WAG').
    Pass model=None to auto-select based on sequence type.

    method='ml' (default): maximum-likelihood estimate via Brent minimisation.
    method='expected': posterior mean E[d|alignment] with FastPhylo prior.
    method='expected_noprior': same but with a flat prior (for comparison).
    Both expected methods use a log-spaced grid of *n_points* distances (protein only).
    """
    if not alignment:
        raise ValueError("alignment is empty")

    # Infer sequence type from first sequence
    seq_type = alignment[0].seq_type

    if model is None:
        if seq_type in (SeqType.DNA, SeqType.RNA, SeqType.Unknown):
            model = _DEFAULT_DNA_MODEL
        else:
            model = _DEFAULT_PROTEIN_MODEL

    names = [s.accession for s in alignment]
    seqs = [s.data for s in alignment]

    if model in _DNA_MODELS:
        cpp_dm = _fastphylo.compute_dna_distances(names, seqs, model=model)
        return DistanceMatrix(cpp_dm)

    # Protein distance
    if method == "expected":
        return _protein_expected_distance_matrix(names, seqs, model, n_points, use_prior=True)
    if method == "expected_noprior":
        return _protein_expected_distance_matrix(names, seqs, model, n_points, use_prior=False)
    return _protein_distance_matrix(names, seqs, model)


def _protein_distance_matrix(
    names: list[str],
    seqs: list[str],
    model: str,
) -> DistanceMatrix:
    """Compute protein pairwise distances using a substitution model."""
    from .protein import compute_protein_distances
    cpp_dm = compute_protein_distances(names, seqs, model)
    return DistanceMatrix(cpp_dm)


def _protein_expected_distance_matrix(
    names: list[str],
    seqs: list[str],
    model: str,
    n_points: int,
    use_prior: bool = True,
) -> DistanceMatrix:
    """Compute expected protein distances E[d|alignment]."""
    from .protein import compute_protein_expected_distances
    cpp_dm = compute_protein_expected_distances(names, seqs, model, n_points, use_prior=use_prior)
    return DistanceMatrix(cpp_dm)
