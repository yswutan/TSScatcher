#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

input_bam= sample_id= barcodes= type= output= threads=5
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) input_bam=$2; shift 2;;
    --id) sample_id=$2; shift 2;;
    --barcodes) barcodes=$2; shift 2;;
    --type) type=$2; shift 2;;
    --output) output=$2; shift 2;;
    --threads) threads=$2; shift 2;;
    *) die "unknown preprocess option: $1";;
  esac
done
[[ "$type" == long || "$type" == short ]] || die "--type must be long or short"
[[ -n "$input_bam" && -n "$sample_id" && -n "$barcodes" && -n "$output" ]] || die "preprocess requires --bam, --id, --barcodes, --type and --output"
for cmd in samtools awk grep python3; do require_cmd "$cmd"; done
require_file "$input_bam"; require_file "$barcodes"
make_dir "$output"
work=$(mktemp -d "$output/.preprocess.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

samtools view -H "$input_bam" > "$work/header.sam"
if [[ "$type" == long ]]; then flag=(-f 64); else flag=(-F 64); fi
samtools view -@ "$threads" "${flag[@]}" "$input_bam" | LC_ALL=C grep 'xf:i:25' > "$work/body.sam" || true
cat "$work/header.sam" "$work/body.sam" | samtools view -@ "$threads" -b -o "$work/umi.bam" -
[[ $(samtools view -c "$work/umi.bam") -gt 0 ]] || die "no reads passed xf:i:25 and read-type filtering for $sample_id"
read_stream "$barcodes" | awk 'NF {print $1}' > "$work/barcodes.txt"
samtools view -@ "$threads" "$work/umi.bam" | awk -F '\t' -v barcode_file="$work/barcodes.txt" '
BEGIN {while ((getline line < barcode_file)>0) keep[line]=1}
{for(i=12;i<=NF;i++) if($i~/^CB:Z:/){bc=substr($i,6); if(bc in keep) print; break}}' > "$work/filtered.sam"
cat "$work/header.sam" "$work/filtered.sam" | samtools view -@ "$threads" -b -o "$work/filtered_cb.bam" -
[[ $(samtools view -c "$work/filtered_cb.bam") -gt 0 ]] || die "no reads matched the selected cell barcodes for $sample_id"

final="$output/${sample_id}.filtered.bam"
if [[ "$type" == long ]]; then
  samtools index -@ "$threads" "$work/filtered_cb.bam"
  python3 "$script_dir/FilterR1TSO.py" "$work/filtered_cb.bam" "$final"
else
  mv "$work/filtered_cb.bam" "$final"
fi
samtools index -@ "$threads" "$final"
require_file "$final"
note "preprocessed BAM: $final"
