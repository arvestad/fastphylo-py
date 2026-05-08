"""Multiple sequence alignment via pyfamsa (optional dependency)."""

from __future__ import annotations

from .sequences import Alignment, Sequence, SequenceCollection

_SCORING_MATRICES: dict[str, object] = {}  # populated lazily on first import


def _resolve_matrix(name: str) -> object:
    """Return the pyfamsa scoring matrix object for the given name."""
    if not _SCORING_MATRICES:
        try:
            import pyfamsa as _pf
            _SCORING_MATRICES.update({
                "PFASUM43": _pf.PFASUM43,
                "PFASUM60": _pf.PFASUM60,
                "PFASUM31": _pf.PFASUM31,
                "MIQS":     _pf.MIQS,
            })
        except ImportError:
            pass
    matrix = _SCORING_MATRICES.get(name)
    if matrix is None:
        available = list(_SCORING_MATRICES) or ["PFASUM43", "PFASUM60", "PFASUM31", "MIQS"]
        raise ValueError(
            f"Unknown scoring matrix {name!r}. "
            f"Available: {', '.join(available)}"
        )
    return matrix


def align(
    collection: SequenceCollection,
    *,
    scoring_matrix: str = "PFASUM43",
) -> Alignment:
    """Align an unaligned SequenceCollection using FAMSA.

    Requires the optional alignment dependency::

        pip install fastphylo[align]

    The scoring matrix is always specified explicitly so results are
    reproducible regardless of the pyfamsa version's default.

    Parameters
    ----------
    collection:
        Unaligned (or already aligned) sequences.
    scoring_matrix:
        Name of the FAMSA scoring matrix: ``"PFASUM43"`` (default),
        ``"PFASUM60"``, ``"PFASUM31"``, or ``"MIQS"``.

    Returns
    -------
    Alignment
        All sequences at equal length with gap characters inserted.
    """
    try:
        import pyfamsa
    except ModuleNotFoundError as exc:
        raise ImportError(
            "pyfamsa is required for alignment; "
            "install it with: pip install 'fastphylo[align]'"
        ) from exc

    matrix_obj = _resolve_matrix(scoring_matrix)

    # Build pyfamsa input sequences (bytes ids and sequences)
    famsa_seqs = [
        pyfamsa.Sequence(s.accession.encode(), s.data.encode())
        for s in collection
    ]

    aligner = pyfamsa.Aligner(scoring_matrix=matrix_obj)
    gapped = aligner.align(famsa_seqs)

    # Build a name → original Sequence map for metadata preservation
    orig: dict[str, Sequence] = {s.accession: s for s in collection}

    aligned_seqs: list[Sequence] = []
    for gs in gapped:
        accession = gs.id.decode()
        data = gs.sequence.decode()
        src = orig.get(accession)
        aligned_seqs.append(Sequence(
            accession,
            data,
            description=src.description if src else "",
            seq_type=src.seq_type if src else None,
            organism=src.organism if src else None,
            accession2=src.accession2 if src else None,
        ))

    description = collection.description if isinstance(collection, Alignment) else None
    return Alignment(aligned_seqs, description=description)
