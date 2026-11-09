#!/usr/bin/env bash
set -euo pipefail
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${script_dir}/common.sh"

sample_id= fastqs= fastq_sample= type= transcriptome= output= threads=10 memory=64 cellranger_bin=cellranger confirm_mismatch=no type_check_reads=10000 short_r1_max_length=35
while [[ $# -gt 0 ]]; do
  case "$1" in
    --id) sample_id=$2; shift 2;; --fastqs) fastqs=$2; shift 2;; --sample) fastq_sample=$2; shift 2;;
    --type) type=$2; shift 2;; --transcriptome) transcriptome=$2; shift 2;; --output) output=$2; shift 2;;
    --threads) threads=$2; shift 2;; --memory) memory=$2; shift 2;; --cellranger) cellranger_bin=$2; shift 2;;
    --confirm-type-mismatch) confirm_mismatch=yes; shift;; --type-check-reads) type_check_reads=$2; shift 2;;
    --short-r1-max-length) short_r1_max_length=$2; shift 2;;
    *) die "unknown map option: $1";;
  esac
done
[[ "$type" == long || "$type" == short ]] || die "--type must be long or short"
for v in sample_id fastqs fastq_sample transcriptome output; do [[ -n "${!v}" ]] || die "missing map argument: $v"; done
require_cmd "$cellranger_bin"; require_cmd python3; [[ -d "$fastqs" ]] || die "FASTQ directory not found: $fastqs"; [[ -d "$transcriptome" ]] || die "Cell Ranger reference not found: $transcriptome"
fastqs=$(realpath "$fastqs"); transcriptome=$(realpath "$transcriptome"); output=$(realpath -m "$output")
check_args=(--fastqs "$fastqs" --sample "$fastq_sample" --type "$type" --reads "$type_check_reads" --short-r1-max-length "$short_r1_max_length")
[[ "$confirm_mismatch" == yes ]] && check_args+=(--confirm-type-mismatch)
python3 "$script_dir/check_fastq.py" "${check_args[@]}"
chemistry=SC5P-PE; [[ "$type" == short ]] && chemistry=SC5P-R2
make_dir "$output"; cd "$output"
"$cellranger_bin" count --id="$sample_id" --fastqs="$fastqs" --sample="$fastq_sample" --transcriptome="$transcriptome" --localcores="$threads" --localmem="$memory" --nosecondary --chemistry="$chemistry"
require_file "$output/$sample_id/outs/possorted_genome_bam.bam"
note "Cell Ranger mapping completed for $sample_id ($chemistry)"
