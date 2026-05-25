"""
Run AliSim (IQ-TREE2) to produce a random tree and simulated alignment.

Outputs:
  <prefix>.fa        — FASTA alignment
  <prefix>.treefile  — Newick tree with branch lengths
"""
import argparse
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix",       required=True)
    ap.add_argument("--model",        required=True, help="IQ-TREE2 model name")
    ap.add_argument("--dtype",        required=True, choices=["dna", "protein"])
    ap.add_argument("--length",       type=int, required=True)
    ap.add_argument("--n-sequences",  type=int, required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.prefix), exist_ok=True)

    seqtype = "DNA" if args.dtype == "dna" else "AA"
    cmd = [
        "iqtree2",
        "--alisim", args.prefix,
        "-m", args.model,
        "--seqtype", seqtype,
        "-t", f"RANDOM{{u,{args.n_sequences}}}",
        "--length", str(args.length),
        "--out-format", "fasta",
        "--redo",
        "-quiet",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    # AliSim writes <prefix>.fa and keeps the tree as <prefix>.treefile.
    # Verify both exist.
    fa_path   = args.prefix + ".fa"
    tree_path = args.prefix + ".treefile"
    if not os.path.exists(fa_path):
        print(f"ERROR: expected FASTA not found: {fa_path}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(tree_path):
        print(f"ERROR: expected treefile not found: {tree_path}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
