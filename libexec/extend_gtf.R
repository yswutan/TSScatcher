options(scipen = 999)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(rtracklayer))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: extend_gtf.R <input.gtf> <output.gtf> <extend_bp>")
input_gtf <- args[1]
output_gtf <- args[2]
extend_bp <- as.integer(args[3])
if (!file.exists(input_gtf)) stop("input GTF not found")
if (extend_bp < 0) stop("extend_bp must be non-negative")

gg <- import(input_gtf)
ggf <- as.data.frame(gg)
required <- c("seqnames", "strand", "gene_id", "type")
if (!all(required %in% colnames(ggf))) stop("GTF must contain gene_id and feature type metadata")
genenames <- unique(ggf[, c("seqnames", "strand", "gene_id")])
pieces <- vector("list", nrow(genenames))

for (i in seq_len(nrow(genenames))) {
  geneinfo <- gg[gg$gene_id %in% genenames$gene_id[i] &
                 as.character(seqnames(gg)) %in% as.character(genenames$seqnames[i]) &
                 as.character(strand(gg)) %in% as.character(genenames$strand[i])]
  tx <- geneinfo[geneinfo$type == "transcript"]
  if (length(tx) == 0) {
    pieces[[i]] <- geneinfo
    next
  }
  if (as.character(strand(geneinfo)[1]) == "+") {
    tss <- start(tx)
    for (site in tss) {
      idx <- which(start(geneinfo) == site)
      old_start <- start(geneinfo)[idx]
      delta <- pmin(extend_bp, old_start - 1L)
      width(geneinfo)[idx] <- width(geneinfo)[idx] + delta
      start(geneinfo)[idx] <- old_start - delta
    }
  } else if (as.character(strand(geneinfo)[1]) == "-") {
    tss <- end(tx)
    for (site in tss) {
      idx <- which(end(geneinfo) == site)
      end(geneinfo)[idx] <- end(geneinfo)[idx] + extend_bp
    }
  }
  pieces[[i]] <- geneinfo
  if (i %% 1000 == 0) message(i)
}
extended <- do.call(c, pieces)
export(extended, output_gtf)
