#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

sample_id= reads= final_sites= output= counts_output=
while [[ $# -gt 0 ]]; do
  case "$1" in
    --id) sample_id=$2; shift 2;; --reads) reads=$2; shift 2;; --final-sites) final_sites=$2; shift 2;;
    --output) output=$2; shift 2;; --counts-output) counts_output=$2; shift 2;; *) die "unknown quantify option: $1";;
  esac
done
for v in sample_id reads final_sites output; do [[ -n "${!v}" ]] || die "missing quantify argument: $v"; done
for cmd in bedtools gzip awk sort uniq python3; do require_cmd "$cmd"; done
require_file "$reads"; require_file "$final_sites"; make_dir "$output"
tmp=$(mktemp -d "$output/.quantify.XXXXXX"); trap 'rm -rf -- "$tmp"' EXIT
read_stream "$reads" | bedtools intersect -s -wa -wb -a - -b <(read_stream "$final_sites") | awk 'BEGIN{FS=OFS="\t"} NF>=13 && $7!="" {print $11,$7}' | sort -k1,1 -k2,2 | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$1}' > "$tmp/counts.tsv"
python3 "$script_dir/matrix_from_counts.py" --final-sites "$final_sites" --counts "$tmp/counts.tsv" --output "$output"
if [[ -n "$counts_output" ]]; then cp -f "$tmp/counts.tsv" "$counts_output"; fi
note "matrix written for $sample_id: $output"
