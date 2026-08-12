"""Time a single distance_matrix() call, repeated, against whichever fastphylo
is installed in the interpreter running this script.

Invoked as a subprocess by run_benchmarks.py under two different
interpreters (old-vendored-Brent vs. new-libfastphylo-backed fastphylo-py),
so the comparison is a real, Python-API-level, apples-to-apples measurement -
not extrapolated from either side's own internal C++/pytest benchmarks.

Prints one JSON object to stdout: {"times": [seconds, ...]} - one entry per
repetition, min-of-N in the report is what actually gets compared (standard
microbenchmark practice: median/min is far less sensitive to one slow rep
from unrelated system noise than a mean would be).
"""

from __future__ import annotations

import argparse
import json
import time

import fastphylo


def read_records(fasta_path: str) -> tuple[list[str], list[str]]:
    col = fastphylo.read_fasta(fasta_path)
    names = [s.accession for s in col]
    seqs = [s.data for s in col]
    return names, seqs


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--reps", type=int, default=5)
    args = p.parse_args()

    aln = fastphylo.read_fasta(args.fasta)

    times = []
    for _ in range(args.reps):
        t0 = time.perf_counter()
        fastphylo.distance_matrix(aln, model=args.model)
        times.append(time.perf_counter() - t0)

    print(json.dumps({"version": fastphylo.__version__, "times": times}))


if __name__ == "__main__":
    main()
