# evolution

Molecular sequence analysis for Python — fast pairwise distances and
neighbour-joining trees for DNA, RNA, and protein sequences.

```python
import fastphylo

aln  = fastphylo.read("sequences.fasta")          # FASTA, Stockholm, or Phylip
dm   = fastphylo.distance_matrix(aln)             # k2p for DNA, WAG for protein
tree = fastphylo.fnj(dm)
print(tree.to_newick())
```

## Features

- **Reads** FASTA, Stockholm, and Phylip files; format auto-detected from
  extension and content
- **DNA/RNA distances**: Hamming, Jukes-Cantor, Kimura 2-parameter (default),
  Tamura-Nei 93 — computed by a fast C++ backend with SIMD (SSE2 / NEON)
- **Protein distances**: maximum-likelihood estimation under WAG (default), LG,
  JTT, Dayhoff, BLOSUM62, VT, cpREV, MtREV, RtREV, HIVb, HIVw, DCMUT,
  JTT-DCMut, and PMB, using a  Brent optimizer written in C++
- **Tree reconstruction**: Neighbour-Joining, Fast NJ, and BioNJ (with branch
  lengths) via the FastPhylo library
- **Branch length estimation**: `fit_branch_lengths` fits branch lengths to any
  tree topology by L1-minimisation against a distance matrix (requires scipy)
- **Distance matrix access**: integer or taxon-name indexing, `copy()`,
  `zeros()` factory, NumPy and Phylip export
- **Multiple alignment**: optional FAMSA integration for unaligned input
- Pure-Python sequence types with Stockholm annotation support (organism, AC,
  description); edge-set `Tree` with `to_newick()` and `merge()`

## Installation

```bash
pip install fastphylo
```

Optional extras:

```bash
pip install fastphylo[align]    # multiple-sequence alignment (pyfamsa)
pip install scipy               # branch length fitting (fit_branch_lengths)
```

## Quick start

### Aligned input → tree

```python
import fastphylo

aln  = fastphylo.read("alignment.sto")            # Stockholm, FASTA, or Phylip
dm   = fastphylo.distance_matrix(aln, model="k2p")
tree = fastphylo.fnj(dm)
print(tree.to_newick())
```

### Unaligned input → align → tree

```python
import fastphylo

seqs = fastphylo.read("sequences.fasta")          # unaligned OK
aln  = fastphylo.align(seqs)                      # requires fastphylo[align]
dm   = fastphylo.distance_matrix(aln)
tree = fastphylo.fnj(dm)
print(tree.to_newick())
```

### Protein sequences

```python
import fastphylo

aln  = fastphylo.read("proteins.fasta")
dm   = fastphylo.distance_matrix(aln, model="LG")
tree = fastphylo.bionj(dm)                        # BioNJ with branch lengths
print(tree.to_newick())
```

### Branch length fitting

NJ and FNJ return topology only (branch lengths are not computed). Use
`fit_branch_lengths` to fit branch lengths to any tree topology by minimising
the L1 deviation from the distance matrix:

```python
import fastphylo

aln  = fastphylo.read("alignment.fasta")
dm   = fastphylo.distance_matrix(aln)
tree = fastphylo.fnj(dm)                          # fast topology, no lengths

tree = fastphylo.fit_branch_lengths(tree, dm)     # requires scipy
print(tree.to_newick())                           # now includes branch lengths
```

This solves a linear program: branch lengths are chosen to minimise
`sum |path_distance(i,j) − dm[i,j]|` over all leaf pairs, subject to
non-negative branch lengths.

### Distance matrix access

```python
dm = fastphylo.distance_matrix(aln)

# Integer or name-based indexing
d = dm[0, 1]
d = dm["human", "mouse"]

# Set elements
dm["human", "mouse"] = 0.15
dm["mouse", "human"] = 0.15

# Build a matrix manually
dm = fastphylo.DistanceMatrix.zeros(["human", "mouse", "rat"])
dm["human", "mouse"] = dm["mouse", "human"] = 0.15
dm["human", "rat"]   = dm["rat",   "human"] = 0.22
dm["mouse", "rat"]   = dm["rat",   "mouse"] = 0.08

# Copy, NumPy array, Phylip string
dm2  = dm.copy()
arr  = dm.to_numpy()
text = dm.to_phylip()
```

## API overview

| Function / class | Description |
|---|---|
| `fastphylo.read(path)` | Read FASTA / Stockholm / Phylip → `SequenceCollection` or `Alignment` |
| `fastphylo.align(seqs)` | Align with FAMSA → `Alignment` |
| `fastphylo.distance_matrix(aln, model=…)` | Compute pairwise distances → `DistanceMatrix` |
| `evolution.nj(dm)` / `fnj(dm)` / `bionj(dm)` | Tree reconstruction → `Tree` |
| `fastphylo.fit_branch_lengths(tree, dm)` | L1-optimal branch lengths for a given topology |
| `DistanceMatrix.zeros(names)` | Create an all-zero matrix with taxon names |
| `DistanceMatrix.copy()` | Deep copy |
| `DistanceMatrix.to_numpy()` | Export as NumPy array |
| `DistanceMatrix.to_phylip()` | Export as Phylip-format string |
| `Tree.to_newick()` | Newick string |
| `Tree.merge(other)` | Union of two edge-set trees |

### Distance models

| Sequences | Model string | Notes |
|---|---|---|
| DNA / RNA | `"hamming"` | Raw mismatch count |
| DNA / RNA | `"jc"` | Jukes-Cantor |
| DNA / RNA | `"k2p"` *(default)* | Kimura 2-parameter |
| DNA / RNA | `"tn93"` | Tamura-Nei 93 |
| Protein | `"WAG"` *(default)* | Whelan & Goldman |
| Protein | `"LG"`, `"JTT"`, `"Dayhoff"`, … | 14 models total |

RNA is handled transparently (U → T at the C++ boundary).

## Requirements

- Python ≥ 3.12
- NumPy ≥ 1.24
- A C compiler (for the bundled FastPhylo extension, built automatically by pip)

Optional:
- pyfamsa ≥ 0.6.0 — multiple-sequence alignment (`pip install fastphylo[align]`)
- scipy — branch length fitting via `fit_branch_lengths` (`pip install scipy`)

## License

GPLv3 — see [FastPhylo](https://github.com/arvestad/FastPhylo) for the
upstream C++ library.
