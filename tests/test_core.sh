#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/a.bed" <<'EOF'
chr1	100	110	a1	.	+
chr1	102	109	a1_duplicate	.	+
chr1	500	510	a2	.	-
EOF
cat > "$tmp/b.bed" <<'EOF'
chr1	105	115	b1	.	+
chr1	505	515	b2	.	-
EOF
cat > "$tmp/manifest.tsv" <<EOF
sample_id	peak_file
A	$tmp/a.bed
B	$tmp/b.bed
EOF
"$root/bin/tsscatcher" merge --manifest "$tmp/manifest.tsv" --output-bed "$tmp/final.bed.gz" --output-support "$tmp/support.tsv.gz" --min-replicates 2 --merge-distance 100
[[ $(gzip -cd "$tmp/final.bed.gz" | wc -l) -eq 2 ]]
gzip -cd "$tmp/final.bed.gz" | awk 'BEGIN{FS="\t"} $5!=2{exit 1}'

cat > "$tmp/counts.tsv" <<'EOF'
TSS0000001	A#BC1	2
TSS0000002	A#BC1	1
EOF
python3 "$root/libexec/matrix_from_counts.py" --final-sites "$tmp/final.bed.gz" --counts "$tmp/counts.tsv" --output "$tmp/matrix"
[[ -s "$tmp/matrix/matrix.mtx.gz" ]]
[[ -s "$tmp/matrix/barcodes.tsv.gz" ]]
[[ -s "$tmp/matrix/features.tsv.gz" ]]

mkdir -p "$tmp/fastq"
{
  printf '@r1\n'; printf 'A%.0s' {1..26}; printf '\n+\n'; printf 'I%.0s' {1..26}; printf '\n'
} | gzip -c > "$tmp/fastq/S1_R1_001.fastq.gz"
{
  printf '@r1\n'; printf 'A%.0s' {1..90}; printf '\n+\n'; printf 'I%.0s' {1..90}; printf '\n'
} | gzip -c > "$tmp/fastq/S1_R2_001.fastq.gz"
"$root/bin/tsscatcher" check-fastq --fastqs "$tmp/fastq" --sample S1 --type short --reads 10 --short-r1-max-length 35 >/dev/null
if "$root/bin/tsscatcher" check-fastq --fastqs "$tmp/fastq" --sample S1 --type long --reads 10 --short-r1-max-length 35 >/dev/null 2>&1; then
  echo "type mismatch was not rejected" >&2
  exit 1
fi
mkdir -p "$tmp/transcriptome" "$tmp/map_output"
cat > "$tmp/mock_cellranger" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
id=
for arg in "$@"; do case "$arg" in --id=*) id=${arg#--id=};; esac; done
[[ -n "$id" ]]
mkdir -p "$id/outs"
printf 'mock-bam\n' > "$id/outs/possorted_genome_bam.bam"
EOF
chmod +x "$tmp/mock_cellranger"
"$root/bin/tsscatcher" map --id S1 --fastqs "$tmp/fastq" --sample S1 --type short --transcriptome "$tmp/transcriptome" --output "$tmp/map_output" --cellranger "$tmp/mock_cellranger" >/dev/null
[[ -s "$tmp/map_output/S1/outs/possorted_genome_bam.bam" ]]
echo CORE_TESTS_OK
