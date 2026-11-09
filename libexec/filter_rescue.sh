#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

sample_id= input= fasta= chrom_sizes= reference_tss= output_dir= threshold=0.5 rescue_distance=100
while [[ $# -gt 0 ]]; do
  case "$1" in
    --id) sample_id=$2; shift 2;;
    --input) input=$2; shift 2;;
    --fasta) fasta=$2; shift 2;;
    --chrom-sizes) chrom_sizes=$2; shift 2;;
    --reference-tss) reference_tss=$2; shift 2;;
    --output-dir) output_dir=$2; shift 2;;
    --dl-threshold) threshold=$2; shift 2;;
    --rescue-distance) rescue_distance=$2; shift 2;;
    *) die "unknown DL/rescue option: $1";;
  esac
done
[[ -n "$sample_id" && -n "$input" && -n "$fasta" && -n "$chrom_sizes" && -n "$reference_tss" && -n "$output_dir" ]] || die "DL/rescue arguments are incomplete"
for cmd in bedtools awk sort paste python3; do require_cmd "$cmd"; done
for file in "$input" "$fasta" "$chrom_sizes" "$reference_tss"; do require_file "$file"; done
make_dir "$output_dir"
work=$(mktemp -d "$output_dir/.dl_rescue.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

awk 'BEGIN{FS=OFS="\t"} NR==FNR{size[$1]=$2;next} {len=$3-$2; left=int((600-len)/2); s=$2-left; e=$3+left; if(s<0)s=0; if(e-s<600)e=s+600; if(e>size[$1]){e=size[$1];s=e-600;if(s<0)s=0} print $1,s,e,$4,$5,$6}' "$chrom_sizes" "$input" > "$work/windows.bed"
bedtools getfasta -fi "$fasta" -bed "$work/windows.bed" -s -fo "$work/windows.fa"
python3 "$script_dir/evalue_predict.py" --input_file "$work/windows.fa" --output_file "$work/probability.txt"
[[ $(wc -l < "$input") -eq $(wc -l < "$work/probability.txt") ]] || die "DL probability count does not match peak count"
paste "$input" "$work/probability.txt" > "$work/scored.tsv"
awk -v p="$threshold" 'BEGIN{FS=OFS="\t"} $7>p {print $1,$2,$3,$4,$5,$6}' "$work/scored.tsv" > "$work/kept.bed"
awk -v p="$threshold" 'BEGIN{FS=OFS="\t"} $7<=p {print $1,$2,$3,$4,$5,$6}' "$work/scored.tsv" > "$work/rejected.bed"

bedtools sort -g "$chrom_sizes" -i "$reference_tss" > "$work/reference.sorted.bed"
bedtools sort -g "$chrom_sizes" -i "$work/rejected.bed" > "$work/rejected.sorted.bed"
if [[ -s "$work/rejected.bed" ]]; then
  bedtools closest -sorted -s -d -a "$work/rejected.sorted.bed" -b "$work/reference.sorted.bed" -g "$chrom_sizes" > "$work/closest.tsv"
  awk -v d="$rescue_distance" 'BEGIN{FS=OFS="\t"} $NF>=0 && $NF<=d {print $1,$2,$3,$4,$5,$6}' "$work/closest.tsv" > "$work/rescued.bed"
else
  : > "$work/rescued.bed"
fi
cat "$work/kept.bed" "$work/rescued.bed" | sort -k1,1 -k2,2n -k3,3n > "$output_dir/${sample_id}.dl_filtered_peaks.bed"
awk 'BEGIN{FS=OFS="\t"}{print $4,1} END{print "__TSScatcher_sentinel__",0}' "$work/rescued.bed" > "$work/rescue_ids.tsv"
awk -v p="$threshold" 'BEGIN{FS=OFS="\t"} NR==FNR{if($2==1)rescued[$1]=1;next}{status=($7>p?"dl_kept":(($4 in rescued)?"rescued":"dl_rejected"));print $4,$7,status}' "$work/rescue_ids.tsv" "$work/scored.tsv" > "$work/meta.body.tsv"
{ printf 'peak_id\tprobability\tstatus\n'; cat "$work/meta.body.tsv"; } > "$work/meta.tsv"
gzip -c "$work/meta.tsv" > "$output_dir/${sample_id}.peak_metadata.tsv.gz"
require_file "$output_dir/${sample_id}.dl_filtered_peaks.bed"
