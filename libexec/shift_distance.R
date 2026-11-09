options(warn = 1)
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: shift_distance.R <distances.txt> <sample_id> <max_distance> <output_dir>")
}
distance_file <- args[1]
sample_id <- args[2]
max_distance <- as.numeric(args[3])
output_dir <- args[4]

values <- scan(distance_file, what = numeric(), quiet = TRUE)
values <- values[is.finite(values) & values >= 0 & values < max_distance]
if (length(values) == 0) stop("no valid shift distances remain after filtering")
shift_value <- as.integer(round(median(values)))
density_obj <- density(values)
plot_data <- data.frame(distance = values)
p <- ggplot(plot_data, aes(x = distance)) +
  geom_density(color = "#E69F00") +
  geom_vline(xintercept = shift_value, color = "blue", linetype = "dashed") +
  annotate("label", x = shift_value, y = max(density_obj$y) / 2,
           label = shift_value, vjust = 2) +
  coord_cartesian(xlim = c(0, max_distance)) +
  theme_bw() + xlab("Distance") + ylab("Density") +
  ggtitle(paste0(sample_id, ": median distance to reference TSS"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(output_dir, paste0(sample_id, ".shift_distance.pdf")), p, width = 4, height = 4)
write.table(shift_value, file.path(output_dir, paste0(sample_id, ".shift_distance.txt")),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
