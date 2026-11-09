# Output layout

Each sample is processed as exactly one library type (`long` or `short`). Long
libraries use Cell Ranger chemistry `SC5P-PE`; short libraries use `SC5P-R2`.

```text
results/
  samples/<sample_id>/
    peaks/
      <sample_id>.paraclu_summit.raw.bed
      <sample_id>.pre_dl_peaks.bed
      <sample_id>.dl_filtered_peaks.bed       # DL-kept plus reference-rescued
      <sample_id>.peak_metadata.tsv.gz
    coverage/
      <sample_id>.long_R1.plus.bw
      <sample_id>.long_R1.minus.bw
      # short libraries instead produce before/after-shift plus/minus bigWigs
    shift/                                    # short libraries only
      <sample_id>.shift_distance.pdf
      <sample_id>.shift_distance.txt
      <sample_id>.before_shift.reads.bed.gz
      <sample_id>.after_shift.reads.bed.gz
    matrix/
      matrix.mtx.gz
      barcodes.tsv.gz
      features.tsv.gz
  aggregate/
    final_sites.bed.gz
    final_sites.support.tsv.gz
    matrix/
      matrix.mtx.gz
      barcodes.tsv.gz
      features.tsv.gz
```

The minus-strand bigWigs contain negative coverage for mirrored genome-browser
display. Read names and cell barcodes are retained in read-level BED files, not
in bedGraph files; standard bedGraph column 4 remains numeric coverage.

Coverage bedGraph files are temporary. After `bedGraphToBigWig` succeeds and the
resulting bigWig is verified to be non-empty, the bedGraph is deleted. It is
retained only when conversion fails, so the failed conversion can be diagnosed.
Read-level BED files used for quantification are not deleted with these temporary
coverage files.

`pre_dl_peaks.bed` is BED6 after TSO-like filtering and before deep-learning
filtering. When deep learning is enabled, `dl_filtered_peaks.bed` is BED6 and
contains both model-kept and reference-rescued peaks. No separate DL-only BED is
written. When deep learning is disabled, neither `dl_filtered_peaks.bed` nor a
rescue output is written, and `pre_dl_peaks.bed` enters multi-sample merging.
