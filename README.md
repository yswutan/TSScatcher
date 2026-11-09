# TSScatcher

TSScatcher calls and quantifies strand-aware transcription start site (TSS)
intervals from 10x Genomics 5′ single-cell RNA-seq. It supports one library type
per sample:

- `long`: equal-length R1/R2 FASTQs, mapped with Cell Ranger `SC5P-PE`;
- `short`: a shorter R1 (normally no more than 35 bp), mapped with `SC5P-R2` and
  corrected by transcriptome-based shift estimation.

Long and short samples may be included in the same project. Their sample-level
peaks are combined into one strand-aware set of final sites.

## Requirements

Linux is required. Create the software environment with:

```bash
micromamba create -f env/environment.yml
micromamba activate tsscatcher
```

Cell Ranger is distributed separately by 10x Genomics. The original TSScatcher
workflow used Cell Ranger 7.2.0. Install it before running TSScatcher and ensure
`cellranger` is on `PATH`, or pass its executable using `--cellranger`.

The runtime never installs packages. If a dependency is absent, the command stops
and reports it before the affected step.

Add TSScatcher to `PATH`:

```bash
export PATH=/path/to/TSScatcher/bin:$PATH
```

## Prepare a reference

```bash
tsscatcher mkref \
  --fasta /reference/genome.fa \
  --gtf /reference/genes.gtf \
  --genome-name GRCh38 \
  --extend-bp 300 \
  --threads 20 \
  --output /reference/tsscatcher_GRCh38
```

This builds the Cell Ranger and RSEM/STAR references, chromosome sizes,
strand-aware reference TSS BED, and the GTF extended at transcript 5′ ends.
Genome-scale reference indexes are not stored in this Git repository.

## Sample sheet

The tab-separated file has exactly these first four columns:

```text
sample_id    fastq_dir       fastq_sample    type
sample01     /data/run_01    sample01        long
sample02     /data/run_02    sample02        short
```

Each sample occurs once and has one type. A sample is never counted twice toward
the replicate threshold.

## Run a project

```bash
tsscatcher run \
  --samples samples.tsv \
  --reference /reference/tsscatcher_GRCh38 \
  --output /analysis/project01 \
  --threads 20 \
  --memory 128 \
  --deep-learning yes \
  --dl-threshold 0.5 \
  --rescue-distance 100 \
  --min-replicates 2 \
  --merge-distance 100 \
  --shift-max-distance 1000
```

Single-sample direct input is also supported:

```bash
tsscatcher run \
  --id sample01 \
  --fastqs /data/run_01 \
  --sample sample01 \
  --type long \
  --reference /reference/tsscatcher_GRCh38 \
  --output /analysis/sample01
```

Before Cell Ranger runs, TSScatcher samples reads from both R1 and R2. Equal
median lengths are expected for `long`; unequal lengths with R1 shorter than R2
and no more than 35 bp are expected for `short`. A mismatch stops the run. After
manual inspection, an intentional mismatch can be accepted with:

```bash
--confirm-type-mismatch
```

## Principal parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `--peak-count` | 5 | Paraclu minimum count |
| `--peak-summit` | 3 | Minimum summit support |
| `--deep-learning` | `yes` | Enable DL filtering and reference rescue |
| `--dl-threshold` | 0.5 | DL probability threshold |
| `--rescue-distance` | 100 | Maximum same-strand distance to a reference TSS |
| `--extend-bp` | 300 | Transcript 5′ extension used by short R1 shifting |
| `--shift-max-distance` | 1000 | Maximum distance included in median shift estimation |
| `--min-replicates` | 2 | Required distinct samples for multi-sample projects |
| `--merge-distance` | 100 | Distance used to join retained TSS partitions |

For a one-sample project, the replicate threshold is automatically one.

## Outputs

See [docs/OUTPUTS.md](docs/OUTPUTS.md). Important outputs include:

- original Paraclu summit files;
- TSO-filtered BED6 peaks before DL;
- when DL is enabled, `dl_filtered_peaks.bed`, containing DL-kept plus
  reference-rescued peaks;
- two strand-specific bigWigs for long R1;
- before/after-shift strand-specific bigWigs, shift plot, shift value, and
  read-level BEDs for short R1;
- a single `final_sites.bed.gz` across all long and short samples;
- Cell Ranger-style `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`
  directories for every sample and for the complete project.

Minus-strand bigWigs use negative coverage for mirrored browser display.
Temporary coverage bedGraph files are removed only after a non-empty bigWig has
been created successfully. Read-level BED files with read names and barcodes are
retained for final-site quantification.

## Algorithm notes

The exact confirmed workflow is documented in
[docs/ALGORITHM.md](docs/ALGORITHM.md). Every genomic comparison, rescue,
replicate calculation, merge and quantification is strand-aware.

## Test

The dependency-light core tests exercise FASTQ type validation, distinct-sample
merging and sparse-matrix generation:

```bash
bash tests/test_core.sh
```

Full end-to-end validation additionally requires Cell Ranger, RSEM/STAR, Paraclu,
the supplied model and a real reference/sample dataset.

## License

No redistribution license has been selected yet. Add the intended license before
publishing the repository publicly.
