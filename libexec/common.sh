#!/usr/bin/env bash
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[$(date '+%F %T')] $*" >&2; }
require_cmd() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
require_file() { [[ -s "$1" ]] || die "required file is missing or empty: $1"; }
make_dir() { mkdir -p "$1"; }

read_stream() {
  case "$1" in
    *.gz) gzip -cd -- "$1" ;;
    *) cat -- "$1" ;;
  esac
}

make_bigwig() {
  local input_bed=$1
  local strand=$2
  local chrom_sizes=$3
  local output_bw=$4
  local negative=${5:-no}
  local tmp_bg="${output_bw}.tmp.bedgraph"
  local tmp_sorted="${output_bw}.tmp.sorted.bedgraph"

  require_file "$input_bed"
  require_file "$chrom_sizes"
  bedtools genomecov -i "$input_bed" -g "$chrom_sizes" -bg -strand "$strand" > "$tmp_bg"
  if [[ "$negative" == yes ]]; then
    awk 'BEGIN{FS=OFS="\t"} {$4=-$4; print}' "$tmp_bg" > "$tmp_sorted"
  else
    sort -k1,1 -k2,2n "$tmp_bg" > "$tmp_sorted"
  fi
  if [[ "$negative" == yes ]]; then
    sort -k1,1 -k2,2n "$tmp_sorted" -o "$tmp_sorted"
  fi
  bedGraphToBigWig "$tmp_sorted" "$chrom_sizes" "$output_bw"
  require_file "$output_bw"
  rm -f -- "$tmp_bg" "$tmp_sorted"
}

normalize_bed6() {
  awk 'BEGIN{FS=OFS="\t"} NF>=6 {print $1,$2,$3,$4,$5,$6}' "$1"
}
