import math
import io
import pytest
from conftest import DATA

import fastphylo
from fastphylo import read_fasta, read_stockholm, distance_matrix, DistanceMatrix


@pytest.fixture
def dna_aln():
    return read_stockholm(DATA / "dna_aligned.sto")


@pytest.fixture
def protein_aln():
    return read_fasta(DATA / "protein_aligned.fasta")


def test_dna_hamming(dna_aln):
    dm = distance_matrix(dna_aln, model="hamming")
    assert isinstance(dm, DistanceMatrix)
    assert len(dm) == 4
    # FastPhylo hamming returns raw mismatch counts (not fractions).
    # seq1 vs seq2 differ at 1 site; seq1 vs seq4 differ at 2 sites.
    assert dm[0, 1] == pytest.approx(1.0)
    assert dm[0, 3] == pytest.approx(2.0)
    # diagonal is 0
    for i in range(4):
        assert dm[i, i] == pytest.approx(0.0)


def test_dna_jc(dna_aln):
    dm = distance_matrix(dna_aln, model="jc")
    assert isinstance(dm, DistanceMatrix)
    assert len(dm) == 4
    # JC-corrected distances are slightly larger than K2P for these sequences.
    # seq1/seq2 differ at 1 of 30 sites → JC ≈ 0.034
    assert 0.03 < dm[0, 1] < 0.04
    # More diverged pair has larger distance
    assert dm[0, 3] > dm[0, 1]
    for i in range(4):
        assert dm[i, i] == pytest.approx(0.0)


def test_dna_k2p(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    assert isinstance(dm, DistanceMatrix)
    assert len(dm) == 4
    for i in range(4):
        assert dm[i, i] == pytest.approx(0.0)
    # symmetry
    for i in range(4):
        for j in range(4):
            assert dm[i, j] == pytest.approx(dm[j, i], rel=1e-6)


def test_dna_tn93(dna_aln):
    dm = distance_matrix(dna_aln, model="tn93")
    assert isinstance(dm, DistanceMatrix)
    assert len(dm) == 4
    for i in range(4):
        assert dm[i, i] == pytest.approx(0.0)


def test_default_method_dna(dna_aln):
    # Default for DNA should be k2p
    dm_default = distance_matrix(dna_aln)
    dm_k2p = distance_matrix(dna_aln, model="k2p")
    for i in range(4):
        for j in range(4):
            assert dm_default[i, j] == pytest.approx(dm_k2p[i, j], rel=1e-6)


def test_distance_matrix_names(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    assert dm.names() == ["seq1", "seq2", "seq3", "seq4"]


def test_distance_matrix_indexing(dna_aln):
    dm = distance_matrix(dna_aln, model="hamming")
    # Integer tuple indexing
    d01 = dm[0, 1]
    d10 = dm[1, 0]
    assert d01 == pytest.approx(d10, rel=1e-6)
    # Integer setitem
    dm[0, 1] = 0.5
    assert dm[0, 1] == pytest.approx(0.5)


def test_distance_matrix_string_indexing(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    # String get matches integer get
    assert dm["seq1", "seq2"] == pytest.approx(dm[0, 1], rel=1e-12)
    assert dm["seq3", "seq1"] == pytest.approx(dm[2, 0], rel=1e-12)
    # String setitem
    dm["seq1", "seq2"] = 0.99
    assert dm[0, 1] == pytest.approx(0.99)
    assert dm["seq1", "seq2"] == pytest.approx(0.99)
    # Unknown name raises KeyError
    with pytest.raises(KeyError):
        _ = dm["nosuchseq", "seq1"]


def test_distance_matrix_copy(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    dm2 = dm.copy()
    # Same values
    for i in range(len(dm)):
        for j in range(len(dm)):
            assert dm2[i, j] == pytest.approx(dm[i, j], rel=1e-12)
    # Same names
    assert dm2.names() == dm.names()
    # Independent: modifying copy does not affect original
    original_val = dm[0, 1]
    dm2[0, 1] = 999.0
    assert dm[0, 1] == pytest.approx(original_val)


def test_distance_matrix_zeros():
    names = ["alpha", "beta", "gamma"]
    dm = DistanceMatrix.zeros(names)
    assert len(dm) == 3
    assert dm.names() == names
    for i in range(3):
        for j in range(3):
            assert dm[i, j] == pytest.approx(0.0)
    # Can fill in manually via strings
    dm["alpha", "beta"] = 0.1
    dm["beta", "alpha"] = 0.1
    dm["alpha", "gamma"] = 0.3
    dm["gamma", "alpha"] = 0.3
    dm["beta", "gamma"] = 0.2
    dm["gamma", "beta"] = 0.2
    assert dm["alpha", "beta"] == pytest.approx(0.1)
    assert dm["beta", "gamma"] == pytest.approx(0.2)
    assert dm["alpha", "gamma"] == pytest.approx(0.3)


def test_distance_matrix_symmetry(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    n = len(dm)
    for i in range(n):
        for j in range(n):
            assert dm[i, j] == pytest.approx(dm[j, i], rel=1e-6)


def test_distance_matrix_to_numpy(dna_aln):
    numpy = pytest.importorskip("numpy")
    dm = distance_matrix(dna_aln, model="k2p")
    arr = dm.to_numpy()
    assert arr.shape == (4, 4)
    # Diagonal is zero
    assert numpy.allclose(numpy.diag(arr), 0.0)
    # Symmetric
    assert numpy.allclose(arr, arr.T)
    # Values match DistanceMatrix indexing
    for i in range(4):
        for j in range(4):
            assert arr[i, j] == pytest.approx(dm[i, j], rel=1e-6)


def test_k2p_saturation_nan():
    # All-T vs all-C is 100% transversions: K2P formula hits log(0) → nan.
    # FastPhylo returns nan rather than raising.
    aln_text = ">s1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n>s2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
    col = fastphylo.read_fasta(io.StringIO(aln_text))
    # Different-length sequences: wrap in an Alignment manually via the C++ layer
    names = [col[0].accession, col[1].accession]
    seqs  = [col[0].data,      col[1].data]
    from fastphylo._fastphylo import compute_dna_distances, DistMatrix
    cpp_dm = compute_dna_distances(names, seqs, model="k2p")
    from fastphylo import DistanceMatrix
    dm = DistanceMatrix(cpp_dm)
    assert math.isnan(dm[0, 1])


def test_protein_wag(protein_aln):
    dm = distance_matrix(protein_aln, model="WAG")
    assert isinstance(dm, DistanceMatrix)
    assert len(dm) == 4
    assert dm.names() == ["prot1", "prot2", "prot3", "prot4"]
    for i in range(4):
        assert dm[i, i] == pytest.approx(0.0)
    # prot1 vs prot2 differ at 1 of 41 sites → small distance
    assert 0.0 < dm[0, 1] < 0.1
    # symmetry
    for i in range(4):
        for j in range(4):
            assert dm[i, j] == pytest.approx(dm[j, i], rel=1e-6)


def test_protein_lg(protein_aln):
    dm_wag = distance_matrix(protein_aln, model="WAG")
    dm_lg  = distance_matrix(protein_aln, model="LG")
    assert isinstance(dm_lg, DistanceMatrix)
    # For highly diverged pairs (prot1 vs prot4), WAG and LG give different values
    assert abs(dm_lg[0, 3] - dm_wag[0, 3]) > 1e-4
    assert 0.0 < dm_lg[0, 1] < 0.15


def test_default_method_protein(protein_aln):
    dm_default = distance_matrix(protein_aln)
    dm_wag     = distance_matrix(protein_aln, model="WAG")
    for i in range(4):
        for j in range(4):
            assert dm_default[i, j] == pytest.approx(dm_wag[i, j], rel=1e-6)


def test_protein_unknown_model(protein_aln):
    with pytest.raises(ValueError, match="Unknown protein model"):
        distance_matrix(protein_aln, model="NOSUCHMODEL")


# ---------------------------------------------------------------------------
# Cross-validation: C++ Brent vs. pure-Python reference implementation
# ---------------------------------------------------------------------------

def _py_neg_log_likelihood(model, N, t):
    """Pure-Python negative log-likelihood using RateMatrix.get_replacement_probs."""
    import numpy as np
    P_t = model.get_replacement_probs(t)
    ll = 0.0
    for i in range(20):
        for j in range(20):
            nij = N[i, j]
            if nij > 0:
                p = max(float(P_t[i, j]), 1e-300)
                ll += nij * math.log(p)
    return -ll


def _py_brent(f, xa, xb, xtol=0.02, maxiter=20):
    """Brent bounded minimizer matching the C++ implementation."""
    mintol = 1e-11
    cg = 0.3819660112501051
    fulc = xa + cg * (xb - xa)
    nfc = fulc
    xf = fulc
    rat = e = 0.0
    x = xf
    fx = f(x)
    ffulc = fnfc = fx
    xm = 0.5 * (xa + xb)
    tol1 = xtol * abs(x) + mintol
    tol2 = 2.0 * tol1

    for _ in range(maxiter):
        if abs(x - xm) <= tol2 - 0.5 * (xb - xa):
            break
        golden = True
        p = q = r_val = 0.0
        if abs(e) > tol1:
            r_val = (x - nfc) * (fx - ffulc)
            q = (x - fulc) * (fx - fnfc)
            p = (x - fulc) * q - (x - nfc) * r_val
            q = 2.0 * (q - r_val)
            if q > 0:
                p = -p
            else:
                q = -q
            r_val, e = e, rat
            if abs(p) < abs(0.5 * q * r_val) and p > q * (xa - x) and p < q * (xb - x):
                rat = p / q
                xnew = x + rat
                if (xnew - xa) < tol2 or (xb - xnew) < tol2:
                    rat = tol1 if xm > x else -tol1
                golden = False
        if golden:
            e = xa - x if x >= xm else xb - x
            rat = cg * e
        u = x + rat if abs(rat) >= tol1 else x + (tol1 if rat > 0 else -tol1)
        fu = f(u)
        if fu <= fx:
            if u < x:
                xb = x
            else:
                xa = x
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = x, fx
            x, fx = u, fu
        else:
            if u < x:
                xa = u
            else:
                xb = u
            if fu <= fnfc or nfc == x:
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = u, fu
            elif fu <= ffulc or fulc == x or fulc == nfc:
                fulc, ffulc = u, fu
        xm = 0.5 * (xa + xb)
        tol1 = xtol * abs(x) + mintol
        tol2 = 2.0 * tol1
    return x


def _py_count_matrix(seq1, seq2):
    import numpy as np
    from fastphylo.protein import AA_ORDER
    idx = {aa: i for i, aa in enumerate(AA_ORDER)}
    idx.update({aa.lower(): i for i, aa in enumerate(AA_ORDER)})
    N = np.zeros((20, 20))
    for a, b in zip(seq1, seq2):
        i = idx.get(a, -1)
        j = idx.get(b, -1)
        if i >= 0 and j >= 0:
            N[i, j] += 1.0
    return N


def test_protein_cpp_vs_python_wag(protein_aln):
    """C++ protein distance must agree with a pure-Python Brent reference."""
    from fastphylo.protein import RateMatrix
    model = RateMatrix.instantiate("WAG")
    dm_cpp = distance_matrix(protein_aln, model="WAG")

    DELTA = 0.0001
    MAX_DIST = 3.0
    seqs = [s.data for s in protein_aln]
    n = len(seqs)
    for i in range(n):
        for j in range(i + 1, n):
            N = _py_count_matrix(seqs[i], seqs[j])
            d_py = _py_brent(
                lambda t, N=N: _py_neg_log_likelihood(model, N, t),
                DELTA, MAX_DIST,
            )
            assert dm_cpp[i, j] == pytest.approx(d_py, abs=0.021), (
                f"pair ({i},{j}): C++={dm_cpp[i,j]:.6f} vs Python={d_py:.6f}"
            )


def test_protein_cpp_vs_python_lg(protein_aln):
    """Same cross-validation for the LG model."""
    from fastphylo.protein import RateMatrix
    model = RateMatrix.instantiate("LG")
    dm_cpp = distance_matrix(protein_aln, model="LG")

    DELTA = 0.0001
    MAX_DIST = 3.0
    seqs = [s.data for s in protein_aln]
    n = len(seqs)
    for i in range(n):
        for j in range(i + 1, n):
            N = _py_count_matrix(seqs[i], seqs[j])
            d_py = _py_brent(
                lambda t, N=N: _py_neg_log_likelihood(model, N, t),
                DELTA, MAX_DIST,
            )
            assert dm_cpp[i, j] == pytest.approx(d_py, abs=0.021), (
                f"pair ({i},{j}): C++={dm_cpp[i,j]:.6f} vs Python={d_py:.6f}"
            )


def test_to_phylip(dna_aln):
    dm = distance_matrix(dna_aln, model="k2p")
    s = dm.to_phylip()
    lines = s.splitlines()
    # Header: number of taxa
    assert lines[0].strip() == "4"
    # One data line per taxon
    assert len(lines) == 5
    # Each line starts with the (padded) name
    assert lines[1].startswith("seq1")
    assert lines[2].startswith("seq2")
    # Diagonal values are 0.000000
    for i, line in enumerate(lines[1:]):
        parts = line.split()
        assert parts[0] == f"seq{i + 1}"
        assert float(parts[1 + i]) == pytest.approx(0.0)
    # Round-trips: to_phylip output can be read back by read_phylip
    import io
    from fastphylo import read_phylip, Alignment
    # Prepend a dummy sequence length to make a valid sequence phylip...
    # Actually Phylip distance format != sequence format; just verify structure.
    # Symmetry: value at (i,j) matches (j,i)
    for i in range(4):
        for j in range(4):
            row_i = lines[1 + i].split()
            row_j = lines[1 + j].split()
            assert float(row_i[1 + j]) == pytest.approx(float(row_j[1 + i]), rel=1e-6)
