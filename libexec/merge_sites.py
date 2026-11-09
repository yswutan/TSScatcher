#!/usr/bin/env python3
import argparse
import gzip
from collections import defaultdict


def read_manifest(path):
    records = []
    seen = set()
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        required = {"sample_id", "peak_file"}
        if not required.issubset(header):
            raise ValueError("merge manifest requires sample_id and peak_file columns")
        idx = {name: header.index(name) for name in required}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not line.strip():
                continue
            sample = fields[idx["sample_id"]]
            peak_file = fields[idx["peak_file"]]
            if sample in seen:
                raise ValueError(f"sample occurs more than once in merge manifest: {sample}")
            seen.add(sample)
            records.append((sample, peak_file))
    if not records:
        raise ValueError("merge manifest has no samples")
    return records


def load_intervals(records):
    by_key = defaultdict(list)
    for sample, path in records:
        with open(path) as handle:
            for line_no, line in enumerate(handle, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[5] not in {"+", "-"}:
                    raise ValueError(f"invalid BED6 at {path}:{line_no}")
                start, end = int(fields[1]), int(fields[2])
                if end <= start:
                    raise ValueError(f"non-positive interval at {path}:{line_no}")
                by_key[(fields[0], fields[5])].append((start, end, sample))
    return by_key


def supported_partitions(intervals, minimum):
    events = defaultdict(lambda: {"start": [], "end": []})
    for start, end, sample in intervals:
        events[start]["start"].append(sample)
        events[end]["end"].append(sample)
    active_counts = defaultdict(int)
    active = set()
    result = []
    previous = None
    for position in sorted(events):
        if previous is not None and position > previous and len(active) >= minimum:
            result.append((previous, position))
        for sample in events[position]["end"]:
            active_counts[sample] -= 1
            if active_counts[sample] <= 0:
                active.discard(sample)
        for sample in events[position]["start"]:
            active_counts[sample] += 1
            active.add(sample)
        previous = position
    return result


def merge_partitions(parts, distance):
    merged = []
    for start, end in parts:
        if not merged or start - merged[-1][1] > distance:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return merged


def main():
    parser = argparse.ArgumentParser(description="Merge strand-aware sample peaks by distinct-sample support")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-bed", required=True)
    parser.add_argument("--output-support", required=True)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--merge-distance", type=int, default=100)
    args = parser.parse_args()
    records = read_manifest(args.manifest)
    minimum = 1 if len(records) == 1 else args.min_replicates
    if minimum < 1 or args.merge_distance < 0:
        raise ValueError("thresholds must be non-negative and min-replicates at least 1")
    by_key = load_intervals(records)
    final = []
    for (chrom, strand), intervals in by_key.items():
        parts = supported_partitions(intervals, minimum)
        for start, end in merge_partitions(parts, args.merge_distance):
            support = sorted({sample for a, b, sample in intervals if a < end and b > start})
            final.append((chrom, start, end, strand, support))
    final.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    with gzip.open(args.output_bed, "wt") as bed, gzip.open(args.output_support, "wt") as support_out:
        support_out.write("final_tss_id\tsupport_samples\tsample_ids\n")
        for number, (chrom, start, end, strand, support) in enumerate(final, 1):
            feature = f"TSS{number:07d}"
            bed.write(f"{chrom}\t{start}\t{end}\t{feature}\t{len(support)}\t{strand}\n")
            support_out.write(f"{feature}\t{len(support)}\t{','.join(support)}\n")
    if not final:
        raise RuntimeError("no final sites passed the distinct-sample support threshold")


if __name__ == "__main__":
    main()
