#!/usr/bin/env python3
import argparse
import gzip
from pathlib import Path


def features_from_bed(path):
    result = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            result.append((f[3], f"{f[0]}:{f[1]}-{f[2]}:{f[5]}", "TSS"))
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--final-sites", required=True)
    p.add_argument("--counts", required=True, help="feature, barcode, count TSV")
    p.add_argument("--output", required=True)
    args = p.parse_args()
    features = features_from_bed(args.final_sites)
    feature_index = {x[0]: i + 1 for i, x in enumerate(features)}
    entries = []
    barcodes = set()
    with open(args.counts) as handle:
        for line in handle:
            feature, barcode, count = line.rstrip("\n").split("\t")
            if feature not in feature_index:
                raise ValueError(f"unknown final TSS ID in counts: {feature}")
            barcodes.add(barcode)
            entries.append((feature, barcode, int(count)))
    barcode_list = sorted(barcodes)
    barcode_index = {x: i + 1 for i, x in enumerate(barcode_list)}
    entries.sort(key=lambda x: (barcode_index[x[1]], feature_index[x[0]]))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with gzip.open(out / "features.tsv.gz", "wt") as handle:
        for row in features:
            handle.write("\t".join(row) + "\n")
    with gzip.open(out / "barcodes.tsv.gz", "wt") as handle:
        for barcode in barcode_list:
            handle.write(barcode + "\n")
    with gzip.open(out / "matrix.mtx.gz", "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        handle.write(f"{len(features)} {len(barcode_list)} {len(entries)}\n")
        for feature, barcode, count in entries:
            handle.write(f"{feature_index[feature]} {barcode_index[barcode]} {count}\n")


if __name__ == "__main__":
    main()
