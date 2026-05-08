"""Sequence types: Sequence, SequenceCollection, Alignment."""

from __future__ import annotations

import re
from enum import Enum, auto
from typing import Iterator, overload

# Characters that can appear in DNA/RNA sequences (IUPAC ambiguity codes included)
_DNA_RNA_CHARS = frozenset("ACGTURYMKSWHBVDNacgturymkswhbvdn.-")


class SeqType(Enum):
    DNA = auto()
    RNA = auto()
    Protein = auto()
    Unknown = auto()


def _detect_seq_type(data: str) -> SeqType:
    clean = data.replace("-", "").replace(".", "").replace(" ", "")
    if not clean:
        return SeqType.Unknown
    chars = set(clean)
    if chars <= _DNA_RNA_CHARS:
        if chars & set("Uu"):
            return SeqType.RNA
        return SeqType.DNA
    return SeqType.Protein


# Binomial name: first two whitespace-separated words, each capitalised
_BINOMIAL_RE = re.compile(r"([A-Z][a-z]+)\s+([a-z]+(?:\s+[a-z]+)?)")


class Sequence:
    """A single biological sequence with metadata."""

    __slots__ = ("accession", "description", "seq_type", "data", "organism", "accession2")

    def __init__(
        self,
        accession: str,
        data: str,
        *,
        description: str = "",
        seq_type: SeqType | None = None,
        organism: str | None = None,
        accession2: str | None = None,
    ) -> None:
        self.accession = accession
        self.description = description
        self.data = data
        self.seq_type = seq_type if seq_type is not None else _detect_seq_type(data)
        self.organism = organism
        self.accession2 = accession2

    def species(self) -> str | None:
        """Return binomial name extracted from organism field, or None."""
        if self.organism is None:
            return None
        m = _BINOMIAL_RE.search(self.organism)
        if m:
            return f"{m.group(1)} {m.group(2)}"
        return None

    def __len__(self) -> int:
        return len(self.data)

    def __repr__(self) -> str:
        return f"Sequence({self.accession!r}, len={len(self.data)}, type={self.seq_type.name})"


class SequenceCollection:
    """An ordered collection of sequences, indexable by position or accession."""

    def __init__(self, sequences: list[Sequence] | None = None) -> None:
        self._seqs: list[Sequence] = sequences or []
        self._index: dict[str, int] = {s.accession: i for i, s in enumerate(self._seqs)}

    def append(self, seq: Sequence) -> None:
        if seq.accession in self._index:
            raise ValueError(f"Duplicate accession: {seq.accession!r}")
        self._index[seq.accession] = len(self._seqs)
        self._seqs.append(seq)

    @overload
    def __getitem__(self, key: int) -> Sequence: ...
    @overload
    def __getitem__(self, key: str) -> Sequence: ...

    def __getitem__(self, key: int | str) -> Sequence:
        if isinstance(key, int):
            return self._seqs[key]
        return self._seqs[self._index[key]]

    def __len__(self) -> int:
        return len(self._seqs)

    def __iter__(self) -> Iterator[Sequence]:
        return iter(self._seqs)

    def __contains__(self, accession: str) -> bool:
        return accession in self._index

    def __repr__(self) -> str:
        return f"{type(self).__name__}({len(self._seqs)} sequences)"


class Alignment(SequenceCollection):
    """A sequence collection where all sequences have equal length (gaps allowed)."""

    def __init__(
        self,
        sequences: list[Sequence] | None = None,
        *,
        description: str | None = None,
    ) -> None:
        super().__init__(sequences)
        self.description = description
        if sequences:
            lengths = {len(s.data) for s in sequences}
            if len(lengths) > 1:
                raise ValueError(
                    f"Alignment requires equal-length sequences; got lengths {lengths}"
                )

    def append(self, seq: Sequence) -> None:
        if self._seqs:
            expected = len(self._seqs[0].data)
            if len(seq.data) != expected:
                raise ValueError(
                    f"Sequence {seq.accession!r} has length {len(seq.data)}; "
                    f"alignment expects {expected}"
                )
        super().append(seq)

    @property
    def alignment_length(self) -> int:
        return len(self._seqs[0].data) if self._seqs else 0
