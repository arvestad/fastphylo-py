# Speed benchmark: old (vendored) vs. new (libfastphylo-backed) fastphylo-py

`fastphylo_py_integration_plan.md`, Phase 6's required deliverable: a
Python-API-level, reproducible comparison - not extrapolated from
FastPhylo's own internal C++ benchmarks or assumed from the underlying
solver swap alone.

## Method

- **Old**: commit `d9d8671` (the tip of `main` immediately before the
  libfastphylo integration) - vendored, pre-restructuring core/DNA
  code, from-scratch Brent's-method protein ML.
- **New**: commit `0bcb84a` - `extern/fastphylo` submodule, real
  FastPhylo headers, `build_ml_decomposition()`/`calculate_ml_dists()`
  (safeguarded Newton-Raphson) for protein ML.
- Each built with its own `pip install -e .` into its own Python 3.12
  venv (separate `git worktree` checkouts, so both are real, complete
  builds - not one patched in place).
- Timed via `distance_matrix(aln, model=...)`, the real top-level
  Python API a user calls - not any internal function.
- **Interleaved at the process level** (old, new, old, new, ... x5
  rounds, 15 in-process repetitions each = 75 timed calls per side per
  dataset) specifically so neither side is penalized by unrelated
  system load happening to land during one contiguous block of runs -
  the same discipline this branch's own Eigen 5.0.1 measurement used
  after a naive run order gave a misleading first result there.
- Reported number is **min of 75 samples** (standard microbenchmark
  practice - far less sensitive to one slow outlier than a mean would
  be); median is also recorded and agrees closely in every row below,
  so the numbers aren't an artifact of picking the minimum.
- Synthetic data (`generate_data.py`, fixed seed `20260812`): each
  dataset mutates one random ancestor into N taxa at a fixed per-site
  mutation rate, giving realistic (non-zero, non-saturated) pairwise
  distances at every size without needing hand-curated real data.

**Environment**: macOS 26.4.1, Darwin 25.4.0, arm64 (Apple Silicon),
Python 3.14.6 driving the benchmark orchestrator (each side's own
`distance_matrix()` call runs inside its own Python 3.12 venv/build -
see Method above), single machine, single run. Not yet run on Linux or
under CI - see the integration plan's "Remaining" section.

## Results

| dataset | N taxa | old (min) | new (min) | speedup (min) |
|---|---:|---:|---:|---:|
| DNA, k2p | 10 | 16.5 µs | 15.8 µs | 1.04x |
| DNA, k2p | 30 | 58.2 µs | 57.7 µs | 1.01x |
| DNA, k2p | 60 | 145.6 µs | 149.5 µs | 0.97x |
| DNA, k2p | 200 | 2.16 ms | 2.21 ms | 0.98x |
| Protein ML, WAG | 10 | 457 µs | 214 µs | **2.14x** |
| Protein ML, WAG | 30 | 4.33 ms | 1.85 ms | **2.35x** |
| Protein ML, WAG | 60 | 17.9 ms | 7.46 ms | **2.39x** |

Full raw data (min, median, sample counts): `results.json`.

## Interpretation

- **Protein ML: a real, consistent ~2.1-2.4x speedup**, growing
  slightly with N (more pairs amortize the one-time
  `build_ml_decomposition()` eigendecomposition cost further) - this
  is the expected larger win (Eigen + symmetrized decomposition +
  build-once API vs. the old from-scratch Brent/count-matrix code),
  and it held up under real, end-to-end measurement, not just
  component-level intuition.
- **DNA/tree path: no meaningful speedup (0.97x-1.04x, i.e. noise)** -
  also as expected going in: `fillMatrix_K2P`/`computeNJTree` are the
  same algorithm either way, just now linked from a properly-built
  shared library instead of a second, independently-compiled copy of
  the same source. No reason to expect a difference, and direct
  measurement confirms there isn't one - reported honestly rather than
  omitted or rounded up to "no change" without checking.
- The N=10/30/60 DNA rows (16-146 µs) are close to Python's own
  function-call/interpreter overhead floor - included for completeness
  and because they still show the expected old≈new pattern, but the
  N=200 row (2.2 ms) is the one actually above the noise floor, and it
  shows the same ~1x result, which is why a fourth, larger DNA tier was
  added specifically (the first pass without it wasn't trustworthy
  enough to report).

## Reproducing

```
python3 benchmarks/generate_data.py
python3 benchmarks/run_benchmarks.py \
    --old-python /path/to/old-venv/bin/python3 \
    --new-python /path/to/new-venv/bin/python3 \
    --out benchmarks/results.json
```

Building the two venvs this run used:

```
git worktree add /tmp/fastphylo-py-old d9d8671   # or any pre-integration ref
python3.12 -m venv /tmp/old-venv && /tmp/old-venv/bin/pip install pybind11 scikit-build-core
(cd /tmp/fastphylo-py-old && /tmp/old-venv/bin/pip install -e . --no-build-isolation)

python3.12 -m venv /tmp/new-venv && /tmp/new-venv/bin/pip install pybind11 scikit-build-core
pip install -e . --no-build-isolation   # from this repo, with extern/fastphylo checked out
```
