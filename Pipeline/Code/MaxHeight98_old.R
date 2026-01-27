library(lidR)
library(raster)

# Args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript calculate_max_height.R <input_folder> <output_folder> <output_prefix>")
}
input_folder  <- args[1]
output_folder <- args[2]
output_prefix <- args[3]

cat("Step: Calculating Top Canopy Height (98th percentile of Z > 1m non-ground points)...\n")

# List LAS/LAZ files (non-empty)
las_files <- list.files(input_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS files found.")

cat("Reading and merging LAS files...\n")
las_combined <- do.call(rbind, lapply(las_files, function(file) {
  las <- tryCatch(readLAS(file), error = function(e) NULL)
  if (!is.null(las) && !is.empty(las)) return(las) else return(NULL)
}))
if (is.null(las_combined) || is.empty(las_combined)) stop("No valid point cloud data after merging.")

# Filter to non-ground points above 1 meters
cat("Filtering nonground points...")
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

cat("Writing raster")
r1 <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
crs(r1) <- "EPSG:5514"

# ---------- 98th percentile height histogram (all cells) ----------
png(
  filename = file.path(output_folder, paste0(output_prefix, "_height98_hist.png")),
  width    = 1000,
  height   = 800
)
raster::hist(
  r1,
  main = "Max height histogram",
  xlab = "Height (m)",
  ylab = "Frequency"
)
dev.off()

# Extract max value from raster
cat("Extracting max value of a raster")
max_top_canopy_height <- max(values(r1), na.rm = TRUE)
cat("Global Max Top Canopy Height (98th percentile):", max_top_canopy_height, "\n")

cat("Writing final output")
# Write result to output file
output_file <- file.path(output_folder, paste0(output_prefix, "_maxHeight.txt"))
write(max_top_canopy_height, file = output_file, ncolumns = 1)

cat("Processing done. Moving forward with analysis")
