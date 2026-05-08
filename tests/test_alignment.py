import pytest
from conftest import DATA

import fastphylo
from fastphylo import (
    read_fasta, read_stockholm, Alignment, SequenceCollection, SeqType, align,
)

pyfamsa = pytest.importorskip("pyfamsa", reason="pyfamsa not installed")


@pytest.fixture
def unaligned():
    return read_fasta(DATA / "dna_unaligned.fasta")


@pytest.fixture
def stockholm_aln():
    return read_stockholm(DATA / "dna_aligned.sto")


def test_align_returns_alignment(unaligned):
    aln = align(unaligned)
    assert isinstance(aln, Alignment)


def test_align_equal_lengths(unaligned):
    aln = align(unaligned)
    lengths = {len(s.data) for s in aln}
    assert len(lengths) == 1, f"Expected one unique length, got {lengths}"


def test_align_preserves_accessions(unaligned):
    aln = align(unaligned)
    assert {s.accession for s in aln} == {s.accession for s in unaligned}


def test_align_preserves_seq_type(unaligned):
    aln = align(unaligned)
    for seq in aln:
        assert seq.seq_type == SeqType.DNA


def test_align_preserves_metadata(stockholm_aln):
    # Align an already-aligned Stockholm collection — metadata (OS, DE) must survive
    aln = align(stockholm_aln)
    assert isinstance(aln, Alignment)
    s1 = aln["seq1"]
    assert s1.organism == "Homo sapiens"
    assert s1.description == "first sequence"


def test_align_description_from_alignment(stockholm_aln):
    # When input is an Alignment with a description, it is propagated
    assert stockholm_aln.description == "Small DNA test alignment"
    aln = align(stockholm_aln)
    assert aln.description == "Small DNA test alignment"


def test_align_description_none_for_collection(unaligned):
    aln = align(unaligned)
    assert aln.description is None


def test_align_scoring_matrix_pfasum43(unaligned):
    aln = align(unaligned, scoring_matrix="PFASUM43")
    assert isinstance(aln, Alignment)


def test_align_invalid_scoring_matrix(unaligned):
    with pytest.raises(ValueError, match="Unknown scoring matrix"):
        align(unaligned, scoring_matrix="NOSUCHMATRIX")


def test_align_gaps_inserted(unaligned):
    # Unaligned sequences have different lengths; aligned versions must contain gaps
    aln = align(unaligned)
    gap_chars = {"-", "."}
    assert any(c in s.data for s in aln for c in gap_chars), \
        "Expected gap characters in aligned sequences"


def test_align_then_distance_then_tree(unaligned):
    # End-to-end: align → distance_matrix → fnj tree
    aln = align(unaligned)
    dm = fastphylo.distance_matrix(aln, model="k2p")
    tree = fastphylo.fnj(dm)
    assert set(tree.leaves.values()) == {s.accession for s in unaligned}
    nwk = tree.to_newick()
    assert nwk.endswith(";")
