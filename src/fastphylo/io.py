"""I/O: read/write FASTA, Stockholm, and Phylip sequence files."""

from __future__ import annotations

import io
import os
from pathlib import Path
from typing import IO, Union

from .sequences import Alignment, Sequence, SequenceCollection, SeqType

PathOrStream = Union[str, os.PathLike, IO[str]]


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------

def _detect_format(path: PathOrStream) -> str:
    if hasattr(path, "read"):
        # stream — sniff content
        return _sniff_stream(path)
    ext = Path(str(path)).suffix.lower()
    if ext in (".fa", ".fasta", ".fna"):
        return "fasta"
    if ext in (".sto", ".stk", ".stockholm"):
        return "stockholm"
    if ext in (".phy", ".phylip"):
        return "phylip"
    # Unknown extension — open and sniff
    with open(path) as fh:
        return _sniff_stream(fh)


def _sniff_stream(fh: IO[str]) -> str:
    pos = fh.tell() if fh.seekable() else None
    head = fh.read(4096)
    if pos is not None:
        fh.seek(pos)
    else:
        # Can't rewind; caller must handle
        pass
    for line in head.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            return "fasta"
        if line.startswith("# STOCKHOLM"):
            return "stockholm"
        # Phylip: first non-blank line is "N L"
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            return "phylip"
        break
    raise ValueError("Cannot detect sequence file format")


def _open(path_or_stream: PathOrStream, mode: str = "r") -> tuple[IO[str], bool]:
    """Return (file_handle, should_close)."""
    if hasattr(path_or_stream, "read") or hasattr(path_or_stream, "write"):
        return path_or_stream, False
    return open(path_or_stream, mode), True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def read(
    path: PathOrStream,
    *,
    format: str | None = None,
    seq_type: SeqType | None = None,
) -> SequenceCollection | Alignment:
    """Read sequences, auto-detecting format from extension or content."""
    fmt = format or _detect_format(path)
    if fmt == "fasta":
        return read_fasta(path, seq_type=seq_type)
    if fmt == "stockholm":
        return read_stockholm(path, seq_type=seq_type)
    if fmt == "phylip":
        return read_phylip(path, seq_type=seq_type)
    raise ValueError(f"Unknown format: {fmt!r}")


def read_fasta(
    path_or_stream: PathOrStream,
    *,
    seq_type: SeqType | None = None,
) -> SequenceCollection:
    fh, close = _open(path_or_stream)
    try:
        seqs: list[Sequence] = []
        accession = description = None
        chunks: list[str] = []

        for raw in fh:
            line = raw.rstrip("\n\r")
            if line.startswith(">"):
                if accession is not None:
                    seqs.append(Sequence(accession, "".join(chunks),
                                         description=description or "",
                                         seq_type=seq_type))
                header = line[1:]
                space = header.find(" ")
                if space == -1:
                    accession, description = header, ""
                else:
                    accession, description = header[:space], header[space + 1:]
                chunks = []
            else:
                if accession is not None:
                    chunks.append(line.strip())

        if accession is not None:
            seqs.append(Sequence(accession, "".join(chunks),
                                  description=description or "",
                                  seq_type=seq_type))
        return SequenceCollection(seqs)
    finally:
        if close:
            fh.close()


def read_stockholm(
    path_or_stream: PathOrStream,
    *,
    seq_type: SeqType | None = None,
) -> Alignment:
    fh, close = _open(path_or_stream)
    try:
        aln_description: str | None = None
        seq_data: dict[str, list[str]] = {}   # accession -> sequence chunks (ordered)
        seq_order: list[str] = []
        gs_de: dict[str, str] = {}
        gs_ac: dict[str, str] = {}
        gs_os: dict[str, str] = {}

        for raw in fh:
            line = raw.rstrip("\n\r")
            if line.startswith("# STOCKHOLM") or line.startswith("//") or not line:
                continue
            if line.startswith("#=GF DE"):
                aln_description = line[len("#=GF DE"):].strip()
            elif line.startswith("#=GS"):
                # #=GS <id> <tag> <value>
                parts = line.split(None, 3)
                if len(parts) == 4:
                    _, sid, tag, value = parts
                    if tag == "DE":
                        gs_de[sid] = value
                    elif tag == "AC":
                        gs_ac[sid] = value
                    elif tag == "OS":
                        gs_os[sid] = value
            elif line.startswith("#"):
                # All other annotation lines ignored
                continue
            else:
                # Sequence line: <name>  <residues>
                parts = line.split(None, 1)
                if len(parts) == 2:
                    name, residues = parts
                    if name not in seq_data:
                        seq_data[name] = []
                        seq_order.append(name)
                    seq_data[name].append(residues)

        seqs = []
        for name in seq_order:
            data = "".join(seq_data[name])
            seqs.append(Sequence(
                name, data,
                description=gs_de.get(name, ""),
                seq_type=seq_type,
                organism=gs_os.get(name),
                accession2=gs_ac.get(name),
            ))
        return Alignment(seqs, description=aln_description)
    finally:
        if close:
            fh.close()


def read_phylip(
    path_or_stream: PathOrStream,
    *,
    seq_type: SeqType | None = None,
) -> Alignment:
    """Read interleaved or sequential Phylip format."""
    fh, close = _open(path_or_stream)
    try:
        lines = [l.rstrip("\n\r") for l in fh]
    finally:
        if close:
            fh.close()

    # Skip blank leading lines
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        raise ValueError("Empty Phylip file")

    header = lines[i].split()
    if len(header) < 2 or not header[0].isdigit() or not header[1].isdigit():
        raise ValueError(f"Phylip header expected '<N> <L>', got: {lines[i]!r}")
    n_seqs = int(header[0])
    seq_len = int(header[1])
    i += 1

    names: list[str] = []
    chunks: list[list[str]] = [[] for _ in range(n_seqs)]
    seq_idx = 0  # which sequence we're reading (for interleaved)
    first_block = True

    while i < len(lines):
        line = lines[i]
        i += 1
        if not line.strip():
            # blank line separates interleaved blocks
            seq_idx = 0
            first_block = False
            continue
        if first_block and seq_idx < n_seqs and len(names) < n_seqs:
            # First block: name is first 10 chars (padded), rest is sequence
            name = line[:10].strip()
            residues = line[10:].replace(" ", "")
            names.append(name)
            chunks[seq_idx].append(residues)
        else:
            residues = line.replace(" ", "")
            chunks[seq_idx % n_seqs].append(residues)
        seq_idx += 1

    seqs = []
    for j, name in enumerate(names):
        data = "".join(chunks[j])
        seqs.append(Sequence(name, data, seq_type=seq_type))

    return Alignment(seqs)


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_fasta(seqs: SequenceCollection, path_or_stream: PathOrStream) -> None:
    fh, close = _open(path_or_stream, "w")
    try:
        for seq in seqs:
            header = seq.accession
            if seq.description:
                header += " " + seq.description
            fh.write(f">{header}\n")
            data = seq.data
            for start in range(0, len(data), 60):
                fh.write(data[start:start + 60] + "\n")
    finally:
        if close:
            fh.close()


def write_phylip(seqs: SequenceCollection, path_or_stream: PathOrStream) -> None:
    seq_list = list(seqs)
    if not seq_list:
        return
    n = len(seq_list)
    length = len(seq_list[0].data)
    fh, close = _open(path_or_stream, "w")
    try:
        fh.write(f" {n} {length}\n")
        for seq in seq_list:
            # Name field: exactly 10 chars, left-padded with spaces
            name = seq.accession[:10].ljust(10)
            fh.write(f"{name}{seq.data}\n")
    finally:
        if close:
            fh.close()
