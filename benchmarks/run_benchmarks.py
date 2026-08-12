"""Reproducible Python-API-level speed comparison: old (vendored, from-scratch
Brent's-method protein ML, pre-restructuring DNA/tree code) vs. new
(libfastphylo-backed, via git submodule) fastphylo-py.

fastphylo_py_integration_plan.md, Phase 6: this is the required,
dedicated speed deliverable - measured by actually calling
distance_matrix() the way a real user would, in two separate real
interpreters/installs, not by re-quoting FastPhylo's own internal C++
benchmarks or trusting an a priori "should be faster" assumption.
Timing is interleaved (old, new, old, new, ... at the process level, on
top of each process's own internal reps) specifically so one side isn't
penalized by picking up unrelated system load that happened to occur
during a single contiguous block of runs - the same methodology used
earlier on the FastPhylo side of this engagement (see CMakeLists.txt's
Eigen 5.0.1 regression note) after a naive front-to-back run order gave
a misleading first result there.

Usage:
    python3 benchmarks/generate_data.py   # writes benchmarks/data/*.fasta once
    python3 benchmarks/run_benchmarks.py \\
        --old-python /path/to/old-venv/bin/python3 \\
        --new-python /path/to/new-venv/bin/python3 \\
        --out benchmarks/results.json

Both interpreters must have their own fastphylo installed already
(pip install -e . against the respective git ref/venv) - this script
only orchestrates and times, it doesn't build anything.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent
DATA_DIR = HERE / "data"
BENCH_ONE = HERE / "bench_one.py"

DATASETS = [
    ("dna", "k2p", "dna_n10.fasta"),
    ("dna", "k2p", "dna_n30.fasta"),
    ("dna", "k2p", "dna_n60.fasta"),
    ("dna", "k2p", "dna_n200.fasta"),
    ("protein", "WAG", "protein_n10.fasta"),
    ("protein", "WAG", "protein_n30.fasta"),
    ("protein", "WAG", "protein_n60.fasta"),
]

INTERLEAVE_ROUNDS = 5  # process-level repeats, alternating old/new
REPS_PER_PROCESS = 15  # in-process repeats (bench_one.py's own --reps)


def run_one(python_exe: str, fasta: Path, model: str) -> list[float]:
    out = subprocess.run(
        [python_exe, str(BENCH_ONE), "--fasta", str(fasta), "--model", model,
         "--reps", str(REPS_PER_PROCESS)],
        capture_output=True, text=True, check=True,
    )
    return json.loads(out.stdout)["times"]


def bench_dataset(old_python: str, new_python: str, fasta: Path, model: str) -> dict:
    old_times: list[float] = []
    new_times: list[float] = []
    # Interleaved at the process level: old, new, old, new, ...
    for _ in range(INTERLEAVE_ROUNDS):
        old_times.extend(run_one(old_python, fasta, model))
        new_times.extend(run_one(new_python, fasta, model))
    return {
        "old_min": min(old_times), "old_median": statistics.median(old_times),
        "new_min": min(new_times), "new_median": statistics.median(new_times),
        "speedup_min": min(old_times) / min(new_times),
        "speedup_median": statistics.median(old_times) / statistics.median(new_times),
        "n_old_samples": len(old_times), "n_new_samples": len(new_times),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--old-python", required=True, help="interpreter with the pre-integration fastphylo-py installed")
    p.add_argument("--new-python", required=True, help="interpreter with the libfastphylo-backed fastphylo-py installed")
    p.add_argument("--out", default=str(HERE / "results.json"))
    args = p.parse_args()

    if not DATA_DIR.exists():
        print("No benchmarks/data/ - run generate_data.py first.", file=sys.stderr)
        sys.exit(1)

    results = {}
    for kind, model, fname in DATASETS:
        fasta = DATA_DIR / fname
        key = f"{kind}:{model}:{fname}"
        print(f"Benchmarking {key} ...", file=sys.stderr)
        results[key] = bench_dataset(args.old_python, args.new_python, fasta, model)

    Path(args.out).write_text(json.dumps(results, indent=2) + "\n")

    print(f"\n{'dataset':<28}{'old (min s)':>14}{'new (min s)':>14}{'speedup (min)':>16}")
    for key, r in results.items():
        print(f"{key:<28}{r['old_min']:>14.4f}{r['new_min']:>14.4f}{r['speedup_min']:>15.2f}x")

    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
