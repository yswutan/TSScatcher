#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

bam= sample_id= output= chrom_sizes= tsolike= fasta= reference_tss= threads=5 peak_count=5 peak_summit=3 dl=yes dl_threshold=0.5 rescue_distance=100
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) bam=$2; shift 2;; --id) sample_id=$2; shift 2;; --output) output=$2; shift 2;;
    --chrom-sizes) chrom_sizes=$2; shift 2;; --tsolike) tsolike=$2; shift 2;; --fasta) fasta=$2; shift 2;;
    --reference-tss) reference_tss=$2; shift 2;; --threads) threads=$2; shift 2;;
    --peak-count) peak_count=$2; shift 2;; --peak-summit) peak_summit=$2; shift 2;;
    --deep-learning) dl=$2; shift 2;; --dl-threshold) dl_threshold=$2; shift 2;;
    --rescue-distance) rescue_distance=$2; shift 2;; *) die "unknown long peak option: $1";;
  esac
done
[[ "$dl" == yes || "$dl" == no ]] || die "--deep-learning must be yes or no"
for v in bam sample_id output chrom_sizes tsolike fasta reference_tss; do [[ -n "${!v}" ]] || die "missing --${v//_/-}"; done
for cmd in samtools bedtools awk sort groupBy paraclu paraclu-cut perl python3 bedGraphToBigWig gzip; do require_cmd "$cmd"; done
for file in "$bam" "$chrom_sizes" "$tsolike"; do require_file "$file"; done
peaks="$output/peaks"; coverage="$output/coverage"; reads="$output/reads"; work_parent="$output/work"
make_dir "$peaks"; make_dir "$coverage"; make_dir "$reads"; make_dir "$work_parent"
work=$(mktemp -d "$work_parent/.long.XXXXXX")
trap 'rm -rf -- "$work"' EXIT

samtools view -@ "$threads" "$bam" | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/){print $1"\t"substr($i,6);break}}' > "$work/read_barcode.tsv"
bedtools bamtobed -cigar -i "$bam" | awk 'BEGIN{OFS="\t"}{split($4,a,"/");n=a[1];if($6=="+"){s=$2;e=$2+1}else{s=$3-1;e=$3}print $1,s,e,n,".",$6}' > "$work/read_starts.bed"
bedtools sort -g "$chrom_sizes" -i "$work/read_starts.bed" > "$work/read_starts.sorted.bed"
bedtools sort -g "$chrom_sizes" -i "$tsolike" > "$work/tsolike.sorted.bed"
bedtools closest -sorted -a "$work/read_starts.sorted.bed" -b "$work/tsolike.sorted.bed" -d -s -id -D a -g "$chrom_sizes" > "$work/tso_distance.bed"
awk 'BEGIN{FS=OFS="\t"} $NF<=-1 {print $1,$2,$3,$4,$5,$6}' "$work/tso_distance.bed" > "$work/tso_filtered.bed"
awk 'BEGIN{FS=OFS="\t"} NR==FNR{bc[$1]=$2;next}{print $1,$2,$3,$4,$5,$6,bc[$4]}' "$work/read_barcode.tsv" "$work/tso_filtered.bed" | gzip -c > "$reads/${sample_id}.reads.bed.gz"

make_bigwig "$work/tso_filtered.bed" + "$chrom_sizes" "$coverage/${sample_id}.long_R1.plus.bw" no
make_bigwig "$work/tso_filtered.bed" - "$chrom_sizes" "$coverage/${sample_id}.long_R1.minus.bw" yes

awk 'BEGIN{OFS="\t"}{print $1,$6,$2,1}' "$work/tso_filtered.bed" > "$work/paraclu.temp"
sort -k1,1 -k2,2 -k3,3n "$work/paraclu.temp" | groupBy -g 1,2,3 -c 4 -o sum > "$work/paraclu.input"
paraclu "$peak_count" "$work/paraclu.input" > "$work/paraclu.output"
paraclu-cut "$work/paraclu.output" > "$work/paraclu.cut"
awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$1"_"$3"_"$4"_"$2,$6,$2}' "$work/paraclu.cut" | sort -k1,1V -k2,2n -k6,6 > "$work/paraclu.bed"
perl "$script_dir/peak.pl" "$work/paraclu.input" "$work/paraclu.bed" > "$work/paraclu.peaks.tsv"
awk -v summit="$peak_summit" 'BEGIN{OFS="\t"}{$3=$3+1;if($7>=summit)print}' "$work/paraclu.peaks.tsv" | sort -k1,1V -k2,2n > "$peaks/${sample_id}.paraclu_summit.raw.bed"
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6}' "$peaks/${sample_id}.paraclu_summit.raw.bed" > "$peaks/${sample_id}.pre_dl_peaks.bed"

if [[ "$dl" == yes ]]; then
  "$script_dir/filter_rescue.sh" --id "$sample_id" --input "$peaks/${sample_id}.pre_dl_peaks.bed" --fasta "$fasta" --chrom-sizes "$chrom_sizes" --reference-tss "$reference_tss" --output-dir "$peaks" --dl-threshold "$dl_threshold" --rescue-distance "$rescue_distance"
fi
note "long R1 peak calling completed for $sample_id"
