import io
from pathlib import Path
from conftest import DATA
import pytest

import fastphylo
from fastphylo import (
    Alignment, Sequence, SequenceCollection, SeqType,
    read, read_fasta, read_stockholm, read_phylip, write_fasta, write_phylip,
)


def test_read_fasta_dna():
    col = read_fasta(DATA / "dna_aligned.fasta")
    assert isinstance(col, SequenceCollection)
    assert len(col) == 4
    assert col[0].accession == "seq1"
    assert col["seq1"].seq_type == SeqType.DNA
    assert col[0].description == "first sequence"
    assert len(col[0].data) == 30


def test_read_fasta_protein():
    buf = io.StringIO(
        ">prot1 a protein\nACDEFGHIKLMNPQRSTVWY\n"
        ">prot2 another\nMKLLVVIFATCLLLAAQTWA\n"
    )
    col = read_fasta(buf)
    assert len(col) == 2
    assert col[0].seq_type == SeqType.Protein
    assert col["prot1"].description == "a protein"


def test_read_stockholm():
    aln = read_stockholm(DATA / "dna_aligned.sto")
    assert isinstance(aln, Alignment)
    assert len(aln) == 4
    assert aln.description == "Small DNA test alignment"
    seq1 = aln["seq1"]
    assert seq1.organism == "Homo sapiens"
    assert seq1.description == "first sequence"
    assert seq1.seq_type == SeqType.DNA
    assert aln["seq2"].organism == "Mus musculus"
    assert aln.alignment_length == 30


def test_read_phylip():
    aln = read_phylip(DATA / "dna_aligned.phy")
    assert isinstance(aln, Alignment)
    assert len(aln) == 4
    assert aln["seq1"].seq_type == SeqType.DNA
    assert aln.alignment_length == 30
    assert aln[0].data == "ACGTACGTACGTACGTACGTACGTACGTAC"


def test_roundtrip_fasta():
    original = read_fasta(DATA / "dna_aligned.fasta")
    buf = io.StringIO()
    write_fasta(original, buf)
    buf.seek(0)
    recovered = read_fasta(buf)
    assert len(recovered) == len(original)
    for a, b in zip(original, recovered):
        assert a.accession == b.accession
        assert a.data == b.data
        assert a.description == b.description


def test_sequence_type_detection():
    dna = Sequence("d", "ACGTACGT")
    assert dna.seq_type == SeqType.DNA

    rna = Sequence("r", "ACGUACGU")
    assert rna.seq_type == SeqType.RNA

    prot = Sequence("p", "MKLLFVIAC")
    assert prot.seq_type == SeqType.Protein

    # Overridable
    forced = Sequence("f", "ACGTACGT", seq_type=SeqType.Protein)
    assert forced.seq_type == SeqType.Protein


def test_rna_mapped_to_dna():
    # SeqType.RNA should be detected from U/u presence
    rna = Sequence("r", "ACGUACGU")
    assert rna.seq_type == SeqType.RNA
    # The data is stored as-is; C++ boundary does the U→T mapping
    assert "U" in rna.data


def test_species_extraction():
    seq = Sequence("x", "ACGT", organism="Homo sapiens (human)")
    assert seq.species() == "Homo sapiens"

    seq2 = Sequence("y", "ACGT", organism="Mus musculus")
    assert seq2.species() == "Mus musculus"

    seq3 = Sequence("z", "ACGT")
    assert seq3.species() is None


def test_auto_format_detection_fasta():
    result = read(DATA / "dna_aligned.fasta")
    assert isinstance(result, SequenceCollection)
    assert len(result) == 4


def test_auto_format_detection_stockholm():
    result = read(DATA / "dna_aligned.sto")
    assert isinstance(result, Alignment)


def test_auto_format_detection_phylip():
    result = read(DATA / "dna_aligned.phy")
    assert isinstance(result, Alignment)


def test_alignment_enforces_equal_lengths():
    seqs = [
        Sequence("a", "ACGT"),
        Sequence("b", "ACGTTT"),
    ]
    with pytest.raises(ValueError, match="equal-length"):
        Alignment(seqs)


def test_collection_indexing():
    col = read_fasta(DATA / "dna_aligned.fasta")
    # Integer indexing
    assert col[0].accession == "seq1"
    assert col[3].accession == "seq4"
    # String indexing
    assert col["seq2"].accession == "seq2"
    # Out-of-range
    with pytest.raises(IndexError):
        _ = col[99]
    # Unknown accession
    with pytest.raises(KeyError):
        _ = col["nosuchseq"]


def test_stockholm_ac_annotation():
    sto = io.StringIO(
        "# STOCKHOLM 1.0\n"
        "#=GS s1 AC P12345\n"
        "#=GS s1 OS Homo sapiens\n"
        "s1  ACGTACGT\n"
        "//\n"
    )
    aln = read_stockholm(sto)
    assert aln["s1"].accession2 == "P12345"
    assert aln["s1"].organism == "Homo sapiens"
