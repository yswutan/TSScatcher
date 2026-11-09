#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

bam= sample_id= output= chrom_sizes= tsolike= fasta= reference_tss= rsem_ref= extended_gtf=
threads=5 peak_count=5 peak_summit=3 dl=yes dl_threshold=0.5 rescue_distance=100 extend_bp=300 shift_max_distance=1000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) bam=$2; shift 2;; --id) sample_id=$2; shift 2;; --output) output=$2; shift 2;;
    --chrom-sizes) chrom_sizes=$2; shift 2;; --tsolike) tsolike=$2; shift 2;; --fasta) fasta=$2; shift 2;;
    --reference-tss) reference_tss=$2; shift 2;; --rsem-ref) rsem_ref=$2; shift 2;; --extended-gtf) extended_gtf=$2; shift 2;;
    --threads) threads=$2; shift 2;; --peak-count) peak_count=$2; shift 2;; --peak-summit) peak_summit=$2; shift 2;;
    --deep-learning) dl=$2; shift 2;; --dl-threshold) dl_threshold=$2; shift 2;; --rescue-distance) rescue_distance=$2; shift 2;;
    --extend-bp) extend_bp=$2; shift 2;; --shift-max-distance) shift_max_distance=$2; shift 2;;
    *) die "unknown short peak option: $1";;
  esac
done
[[ "$dl" == yes || "$dl" == no ]] || die "--deep-learning must be yes or no"
for v in bam sample_id output chrom_sizes tsolike fasta reference_tss rsem_ref extended_gtf; do [[ -n "${!v}" ]] || die "missing --${v//_/-}"; done
for cmd in samtools bedtools awk sort groupBy paraclu paraclu-cut perl python3 Rscript rsem-calculate-expression gppy bedGraphToBigWig gzip; do require_cmd "$cmd"; done
for file in "$bam" "$chrom_sizes" "$tsolike" "$fasta" "$reference_tss" "$extended_gtf"; do require_file "$file"; done

peaks="$output/peaks"; coverage="$output/coverage"; reads="$output/reads"; shift_dir="$output/shift"; work_parent="$output/work"
make_dir "$peaks"; make_dir "$coverage"; make_dir "$reads"; make_dir "$shift_dir"; make_dir "$work_parent"
work=$(mktemp -d "$work_parent/.short.XXXXXX")
trap 'rm -rf -- "$work"' EXIT
samtools view -@ "$threads" "$bam" | awk '{for(i=12;i<=NF;i++) if($i~/^CB:Z:/){print $1"\t"substr($i,6);break}}' > "$work/read_barcode.tsv"

# Before-shift sites use the read five-prime coordinate, but the paired-read
# alignment strand is reversed to report the biological TSS strand.
bedtools bamtobed -i "$bam" | awk 'BEGIN{OFS="\t"}{if($6=="+"){s=$2;e=$2+1;st="-"}else{s=$3-1;e=$3;st="+"}print $1,s,e,$4,$5,st}' | sort -k1,1 -k2,2n > "$work/before_shift.bed"
awk 'BEGIN{FS=OFS="\t"} NR==FNR{bc[$1]=$2;next}{print $1,$2,$3,$4,$5,$6,bc[$4]}' "$work/read_barcode.tsv" "$work/before_shift.bed" | gzip -c > "$shift_dir/${sample_id}.before_shift.reads.bed.gz"
make_bigwig "$work/before_shift.bed" + "$chrom_sizes" "$coverage/${sample_id}.short_R1.before_shift.plus.bw" no
make_bigwig "$work/before_shift.bed" - "$chrom_sizes" "$coverage/${sample_id}.short_R1.before_shift.minus.bw" yes

samtools fastq -@ "$threads" "$bam" | awk '{if(NR%4==1)sub(/\/2$/,"",$0);print}' > "$work/reads.fq"
rsem-calculate-expression -p "$threads" --star --append-names "$work/reads.fq" "$rsem_ref" "$work/rsem"
require_file "$work/rsem.transcript.bam"
samtools sort -@ "$threads" -o "$work/transcript.sorted.bam" "$work/rsem.transcript.bam"
samtools index -@ "$threads" "$work/transcript.sorted.bam"
python3 "$script_dir/GetReadAssign.bed.py" "$work/transcript.sorted.bam" "$work/transcript.bed" ZW 0
awk -v m="$shift_max_distance" '$2<m{print $2}' "$work/transcript.bed" > "$work/distances.txt"
Rscript "$script_dir/shift_distance.R" "$work/distances.txt" "$sample_id" "$shift_max_distance" "$shift_dir"
shift_value=$(<"$shift_dir/${sample_id}.shift_distance.txt")
add_bp=$(awk -v e="$extend_bp" -v s="$shift_value" 'BEGIN{printf "%d", e-s}')
awk -v a="$add_bp" 'BEGIN{FS=OFS="\t"}{s=$2+a;e=$3+a;if(s<0){s=1;e=2}print $1,s,e,$4,$5,$6}' "$work/transcript.bed" > "$work/transcript.shifted.bed"
gppy tiv2giv -g "$extended_gtf" -i "$work/transcript.shifted.bed" | awk 'BEGIN{FS=OFS="\t"}{if($6=="+"){$3=$2+1}else if($6=="-"){$2=$3-1}print $1,$2,$3,$4,$5,$6}' | sort -k1,1 -k2,2n > "$work/after_shift.bed"
awk 'BEGIN{FS=OFS="\t"} NR==FNR{bc[$1]=$2;next}{print $1,$2,$3,$4,$5,$6,bc[$4]}' "$work/read_barcode.tsv" "$work/after_shift.bed" > "$work/after_shift.with_barcode.bed"
gzip -c "$work/after_shift.with_barcode.bed" > "$shift_dir/${sample_id}.after_shift.reads.bed.gz"
cp -f "$shift_dir/${sample_id}.after_shift.reads.bed.gz" "$reads/${sample_id}.reads.bed.gz"
make_bigwig "$work/after_shift.bed" + "$chrom_sizes" "$coverage/${sample_id}.short_R1.after_shift.plus.bw" no
make_bigwig "$work/after_shift.bed" - "$chrom_sizes" "$coverage/${sample_id}.short_R1.after_shift.minus.bw" yes

awk 'BEGIN{OFS="\t"}{print $1,$6,$2,1}' "$work/after_shift.bed" > "$work/paraclu.temp"
sort -k1,1 -k2,2 -k3,3n "$work/paraclu.temp" | groupBy -g 1,2,3 -c 4 -o sum > "$work/paraclu.input"
paraclu "$peak_count" "$work/paraclu.input" > "$work/paraclu.output"
paraclu-cut "$work/paraclu.output" > "$work/paraclu.cut"
awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$1"_"$3"_"$4"_"$2,$6,$2}' "$work/paraclu.cut" | sort -k1,1V -k2,2n -k6,6 > "$work/paraclu.bed"
perl "$script_dir/peak.pl" "$work/paraclu.input" "$work/paraclu.bed" > "$work/paraclu.peaks.tsv"
awk -v summit="$peak_summit" 'BEGIN{OFS="\t"}{$3=$3+1;if($7>=summit)print}' "$work/paraclu.peaks.tsv" | sort -k1,1V -k2,2n > "$peaks/${sample_id}.paraclu_summit.raw.bed"
bedtools sort -g "$chrom_sizes" -i "$peaks/${sample_id}.paraclu_summit.raw.bed" > "$work/summit.sorted.bed"
bedtools sort -g "$chrom_sizes" -i "$tsolike" > "$work/tsolike.sorted.bed"
bedtools closest -sorted -s -d -a "$work/summit.sorted.bed" -b "$work/tsolike.sorted.bed" -g "$chrom_sizes" > "$work/summit.tso_distance.bed"
awk 'BEGIN{FS=OFS="\t"} $NF>=100{print $1,$2,$3,$4,$5,$6}' "$work/summit.tso_distance.bed" | sort -k1,1 -k2,2n > "$peaks/${sample_id}.pre_dl_peaks.bed"

if [[ "$dl" == yes ]]; then
  "$script_dir/filter_rescue.sh" --id "$sample_id" --input "$peaks/${sample_id}.pre_dl_peaks.bed" --fasta "$fasta" --chrom-sizes "$chrom_sizes" --reference-tss "$reference_tss" --output-dir "$peaks" --dl-threshold "$dl_threshold" --rescue-distance "$rescue_distance"
fi
note "short R1 peak calling completed for $sample_id; shift=$shift_value"
