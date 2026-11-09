# Confirmed workflow

1. Inspect sampled reads from both R1 and R2 FASTQs.
2. Validate the declared type. Equal R1/R2 lengths are expected for `long`;
   unequal lengths with the shorter R1 (normally no more than 35 bp) are expected
   for `short`. A mismatch stops unless `--confirm-type-mismatch` is supplied.
3. Map FASTQs with Cell Ranger (`SC5P-PE` for long, `SC5P-R2` for short).
4. Preprocess the Cell Ranger BAM and selected cell barcodes.
5. Call Paraclu peaks per sample and retain the original summit file unchanged.
6. Apply strand-aware TSO-like filtering to form the BED6 pre-DL peak set.
7. Optionally apply the deep-learning filter, then rescue rejected peaks that are
   on the same strand and within `--rescue-distance` of a TSS derived from the
   user-supplied reference GTF.
8. Combine long and short samples into one strand-aware merge. A sample contributes
   at most once to a partition. With multiple samples, retain partitions supported
   by at least `--min-replicates` distinct sample IDs (default 2); with one sample,
   the threshold is automatically 1. Merge retained partitions within
   `--merge-distance` (default 100 bp).
9. Recompute distinct sample support on the final merged intervals.
10. Quantify every sample against the same final intervals, writing per-sample and
    aggregate Cell Ranger-style sparse matrices.

For short libraries, shift estimation uses the median of distances below
`--shift-max-distance` (default 1000 bp). The reference transcript extension is
controlled by `--extend-bp` (default 300 bp).
