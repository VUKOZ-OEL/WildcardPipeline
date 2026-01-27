library(lidR)
library(raster)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript MaxHeight98.R <input_folder> <output_folder> <output_prefix> <crs_string>")
}
input_folder  <- args[1]
output_folder <- args[2]
output_prefix <- args[3]
crs_string    <- args[4]

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

cat("Step: Calculating Top Canopy Height (98th percentile of Z > 1m non-ground points)...\n")

las_files <- list.files(input_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS/LAZ files found.")

cat("Found", length(las_files), "LAS/LAZ file(s) in", input_folder, "\n")

cat("Reading as LAScatalog (chunked processing for large files)...\n")
# IMPORTANT: read the whole folder, not a single file
ctg <- readLAScatalog(input_folder)

opt_chunk_size(ctg)   <- 1000
opt_chunk_buffer(ctg) <- 10

# IMPORTANT: avoid read-time filtering that can empty whole chunks (state ∅)
opt_filter(ctg) <- ""

p98_canopy <- function(z, cls) {
  zz <- z[cls != 2L & z > 1]              # non-ground AND above 1 m
  if (length(zz) == 0L) return(NA_real_)  # no canopy points in this cell
  as.numeric(stats::quantile(zz, 0.98, na.rm = TRUE, names = FALSE))
}

cat("Computing 98th percentile canopy height raster...\n")
r1 <- grid_metrics(ctg, ~p98_canopy(Z, Classification), res = 1)

crs(r1) <- crs_string

# Histogram only if at least one finite cell exists
vals <- values(r1)
finite_n <- sum(is.finite(vals))
if (finite_n > 0) {
  png(
    filename = file.path(output_folder, paste0(output_prefix, "_height98_hist.png")),
    width    = 1000,
    height   = 800
  )
  raster::hist(r1, main = "Max height histogram", xlab = "Height (m)", ylab = "Frequency")
  dev.off()
} else {
  cat("Histogram skipped: raster has no finite values (no canopy points >1m).\n")
}

cat("Extracting max value of the raster...\n")
max_top_canopy_height <- suppressWarnings(max(vals, na.rm = TRUE))
if (!is.finite(max_top_canopy_height)) max_top_canopy_height <- NA_real_

cat("Global Max Top Canopy Height (98th percentile):", max_top_canopy_height, "\n")

cat("Writing final output...\n")
output_file <- file.path(output_folder, paste0(output_prefix, "_maxHeight.txt"))
write(max_top_canopy_height, file = output_file, ncolumns = 1)

cat("Processing done. Moving forward with analysis.\n")
