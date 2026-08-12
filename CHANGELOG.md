# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed - backend replaced with libfastphylo

fastphylo-py now consumes FastPhylo's C++ core (`libfastphylo`) as a
git submodule (`extern/fastphylo`), instead of vendoring a stale,
pre-restructuring fork of it and separately reimplementing protein
maximum-likelihood distance estimation from scratch.

- Protein ML distances (`distance_matrix(..., method="ml")`, the
  default for protein alignments) are now computed by libfastphylo's
  safeguarded Newton-Raphson solver for every supported model, instead
  of this package's own from-scratch Brent's-method optimizer.
- DNA/tree distance and NJ/FNJ/BioNJ tree reconstruction now call
  FastPhylo's real, current `fillMatrix`/`computeNJTree` rather than
  an independently-drifting vendored copy.
- The public API (`distance_matrix()`, `nj()`/`fnj()`/`bionj()`,
  `Tree`, `DistanceMatrix`) is unchanged - this is a backend swap, not
  an interface change.
- Two real, disclosed behavior differences from the swap, both
  libfastphylo safeguards the old vendored/from-scratch code lacked:
  - Protein ML distances on highly-diverged pairs now saturate to a
    fixed sentinel (5.0) uniformly across models, instead of the old
    Brent search's own per-call cap.
  - K2P DNA distances on fully-saturated (100% transversion) pairs now
    return a safeguarded sentinel (10.0, with a warning) instead of
    `nan`.

### Added

- A pure-Python fallback for protein ML distance estimation
  (`FASTPHYLO_FORCE_PYTHON_ML=1`), using
  `scipy.optimize.minimize_scalar` against the existing numpy
  eigendecomposition, for environments without the compiled
  `_fastphylo` extension's solver.
- `benchmarks/`: a reproducible, Python-API-level speed comparison
  script and its results (see below).

### Performance

Measured at the `distance_matrix()` Python API level (not extrapolated
from FastPhylo's own internal C++ benchmarks), old vendored build vs.
new libfastphylo-backed build, interleaved timing, multiple dataset
sizes - see `benchmarks/results.md` for full methodology and raw data.

- **Protein ML distances: ~2.1-2.4x faster**, growing slightly with
  taxon count (10 taxa: 2.14x; 30 taxa: 2.35x; 60 taxa: 2.39x).
- **DNA/tree distances: no meaningful change** (0.97x-1.04x, i.e.
  within noise) - expected, since it's the same algorithm either way,
  now linked from a properly-built shared library instead of a
  duplicate compiled copy.

## [1.0.0] - 2026-07-29

- Real NJ and BioNJ branch lengths (previously only FNJ computed
  them); NJ/BioNJ correctness fixes.

## [0.2.0] - 2026-07-24

- `RateMatrix.instantiate()` caching and a `copy()` method for
  customizing a model without mutating the shared cached instance.
- Numerical validation framework and an expected-distance estimator.
- Capped `numpy` to `<2.3` for `manylinux2014`/`manylinux_2_17` wheel
  compatibility.

## [0.1.0] - 2026-05-08

- Initial release: DNA/protein distance computation, NJ/FNJ/BioNJ tree
  reconstruction, FASTA/Stockholm I/O, GitHub Actions wheel builds.
