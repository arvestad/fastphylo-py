# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added

- `.github/workflows/check-libfastphylo-updates.yml`: a daily scheduled
  check for new FastPhylo releases. If `extern/fastphylo`'s pinned
  commit is behind FastPhylo's latest release, it bumps the submodule,
  builds and tests fastphylo-py against the new commit, and opens a PR
  either way - the PR states whether the build passed, so cutting a
  FastPhylo release now surfaces as a GitHub notification here
  regardless of outcome, not something that has to be checked for
  manually. Requires the repo's Actions settings to allow the default
  token to create pull requests (Settings -> Actions -> General ->
  Workflow permissions) - not verified end-to-end here, since that's a
  one-time repo setting, not something a workflow run can check itself.

## [1.1.1] - 2026-08-12

### Fixed

- **Linux and macOS wheels were missing for Python 3.13** on the 1.1.0
  PyPI release: `release.yml` pinned `pypa/cibuildwheel@v2.19`, which
  predates Python 3.13's actual release (October 2024) and simply
  doesn't recognize `cp313` as a build target - `CIBW_BUILD`'s
  `"cp313-*"` silently matched nothing, no error or skip message
  anywhere in the build log. Confirmed directly: a local run under
  2.19.2 only ever planned to build `cp312-manylinux_x86_64`; the same
  run under 2.23.4 planned both `cp312-manylinux_x86_64` and
  `cp313-manylinux_x86_64`. Bumped the pin to `v2.23`. Python 3.13
  users installing `fastphylo` before this fix silently fell back to
  building from the sdist (which does work, given a C++ toolchain and
  Eigen/BLAS/LAPACK - see 1.1.0's own sdist fix below - just not the
  smooth wheel-install experience 3.12 users got).

### Added

- `release.yml` now fails loudly, before spending time on the actual
  build, if `cibuildwheel --print-build-identifiers` doesn't include
  every Python version `CIBW_BUILD` expects (`cp312`, `cp313`) -
  confirmed this catches the exact cp313 bug above: reproduced the
  check failing under the old `v2.19` pin, passing under `v2.23`. This
  is the durable fix - it catches any future silent gap (a new CPython
  release, a typo, anything) regardless of whether anyone remembers to
  keep the cibuildwheel pin current.
- `.github/dependabot.yml`: watches every pinned GitHub Action version
  across `.github/workflows/` (including `cibuildwheel`) and opens a PR
  when a newer release is available, so a stale pin like the one above
  gets surfaced for review instead of silently going stale for months.
  Doesn't replace the assertion step above - that one is the actual
  safety net; this is what makes needing it less likely in the first
  place.

## [1.1.0] - 2026-08-12

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

### Fixed

- `fastphylo.__version__` was a hardcoded literal, independent of
  `pyproject.toml`'s own version - the two had already drifted before
  anyone noticed. Now read from the installed package's own metadata.
- CI/release builds: `extern/fastphylo` (the new submodule) wasn't
  being checked out at all (`actions/checkout` doesn't fetch submodules
  by default), and libfastphylo's own new BLAS/LAPACK/Eigen3
  dependency wasn't installed in any build environment (`tests.yml`'s
  host runners, or the Linux release wheel's manylinux2014 container -
  the latter needed `epel-release` first, since a plain `yum install
  eigen3-devel` 404s on the base repos). Every build was failing
  because of this, not because of anything wrong with the backend swap
  itself. Also confirmed the sdist genuinely bundles the submodule's
  files (not just its gitlink) and installs correctly on its own.

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
