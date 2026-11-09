#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

fasta= gtf= output= genome_name=TSScatcher threads=10 extend_bp=300 cellranger_bin=cellranger
while [[ $# -gt 0 ]]; do
  case "$1" in
    --fasta) fasta=$2; shift 2;;
    --gtf) gtf=$2; shift 2;;
    --output) output=$2; shift 2;;
    --genome-name) genome_name=$2; shift 2;;
    --threads) threads=$2; shift 2;;
    --extend-bp) extend_bp=$2; shift 2;;
    --cellranger) cellranger_bin=$2; shift 2;;
    *) die "unknown mkref option: $1";;
  esac
done
[[ -n "$fasta" && -n "$gtf" && -n "$output" ]] || die "mkref requires --fasta, --gtf and --output"
for cmd in samtools Rscript rsem-prepare-reference "$cellranger_bin"; do require_cmd "$cmd"; done
require_file "$fasta"; require_file "$gtf"
fasta=$(realpath "$fasta"); gtf=$(realpath "$gtf"); output=$(realpath -m "$output")
make_dir "$output/genome"; make_dir "$output/rsem"; make_dir "$output/cellranger"
cp -f -- "$fasta" "$output/genome/genome.fa"
cp -f -- "$gtf" "$output/genome/genes.gtf"
samtools faidx "$output/genome/genome.fa"
cut -f1,2 "$output/genome/genome.fa.fai" > "$output/genome/chrom.sizes"
Rscript "$script_dir/extend_gtf.R" "$output/genome/genes.gtf" "$output/genome/genes.extend${extend_bp}.gtf" "$extend_bp"
awk 'BEGIN{FS=OFS="\t"} $3=="transcript" && ($7=="+" || $7=="-") {t=($7=="+"?$4:$5); s=t-1; if(s<0)s=0; k=$1 FS s FS t FS $7; if(!seen[k]++) print $1,s,t,$1":"s"-"t":"$7,".",$7}' "$output/genome/genes.gtf" | sort -k1,1 -k2,2n > "$output/genome/reference_tss.bed"
require_file "$output/genome/genes.extend${extend_bp}.gtf"
require_file "$output/genome/reference_tss.bed"
rsem-prepare-reference -p "$threads" --gtf "$output/genome/genes.gtf" --star "$output/genome/genome.fa" "$output/rsem/reference"
[[ ! -e "$output/cellranger/$genome_name" ]] || die "Cell Ranger reference target already exists: $output/cellranger/$genome_name"
(cd "$output/cellranger" && "$cellranger_bin" mkref --genome="$genome_name" --fasta="$output/genome/genome.fa" --genes="$output/genome/genes.gtf" --nthreads="$threads")
[[ -d "$output/cellranger/$genome_name" ]] || die "Cell Ranger reference directory was not created"
note "reference bundle created at $output"
