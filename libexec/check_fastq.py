#!/usr/bin/env python3
import argparse
import gzip
import statistics
import sys
from pathlib import Path


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def lengths(path, limit):
    result = []
    with open_text(path) as handle:
        while len(result) < limit:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip("\r\n")
            handle.readline()
            handle.readline()
            if seq:
                result.append(len(seq))
    if not result:
        raise ValueError(f"no FASTQ reads found in {path}")
    return result


def locate(directory, sample, mate):
    root = Path(directory)
    patterns = [f"{sample}*_{mate}_*.fastq.gz", f"{sample}*_{mate}_*.fq.gz",
                f"{sample}*{mate}*.fastq.gz", f"{sample}*{mate}*.fq.gz",
                f"{sample}*{mate}*.fastq", f"{sample}*{mate}*.fq"]
    found = []
    for pattern in patterns:
        found.extend(root.glob(pattern))
    found = sorted(set(found))
    if not found:
        raise FileNotFoundError(f"cannot find {mate} FASTQ for sample {sample} under {directory}")
    return found[0]


def main():
    parser = argparse.ArgumentParser(description="Validate TSScatcher long/short library type from R1/R2 lengths")
    parser.add_argument("--fastqs", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--type", required=True, choices=["long", "short"])
    parser.add_argument("--reads", type=int, default=10000)
    parser.add_argument("--short-r1-max-length", type=int, default=35)
    parser.add_argument("--confirm-type-mismatch", action="store_true")
    args = parser.parse_args()

    r1_path = locate(args.fastqs, args.sample, "R1")
    r2_path = locate(args.fastqs, args.sample, "R2")
    r1 = lengths(r1_path, args.reads)
    r2 = lengths(r2_path, args.reads)
    r1_med = statistics.median(r1)
    r2_med = statistics.median(r2)
    inferred = "long" if r1_med == r2_med else ("short" if r1_med < r2_med and r1_med <= args.short_r1_max_length else "ambiguous")
    mismatch = inferred == "ambiguous" or inferred != args.type

    print(f"sample={args.sample}\tdeclared={args.type}\tinferred={inferred}\tR1_median={r1_med:g}\tR2_median={r2_med:g}\tR1_file={r1_path}\tR2_file={r2_path}")
    if mismatch:
        message = (f"FASTQ length check does not match declared type '{args.type}': "
                   f"R1 median={r1_med:g}, R2 median={r2_med:g}, inferred={inferred}")
        if args.confirm_type_mismatch:
            print(f"WARNING: {message}; continuing because --confirm-type-mismatch was supplied", file=sys.stderr)
        else:
            print(f"ERROR: {message}; inspect the FASTQs and rerun with --confirm-type-mismatch only if intentional", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
