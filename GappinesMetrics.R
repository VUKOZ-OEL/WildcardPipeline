# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript GappinesMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix>")
}

las_folder <- args[1]
output_folder <- args[2]
max_height_file <- args[3]
base_name <- args[4]

# Read global max height
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global)) stop("Invalid max height value in the file.")

# Output raster folder
raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)

# Load LAS files
las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS/LAZ files found.")

# Combine LAS files
cat("Reading and merging LAS files...\n")
las_combined <- do.call(rbind, lapply(las_files, function(file) {
  las <- tryCatch(readLAS(file), error = function(e) NULL)
  if (!is.null(las) && !is.empty(las)) return(las) else return(NULL)
}))
if (is.null(las_combined) || is.empty(las_combined)) stop("No valid point cloud data after merging.")

# Filter to non-ground points above 1 meter
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

# Compute binary presence rasters (1 = gap, 0 = canopy)
# Deep Gap: cell has no points > 2m
r_deep_gap_bin <- grid_metrics(las_nonground, ~as.numeric(length(Z[Z > 2]) == 0), res = 1)

# Gap: threshold = 25% of site max height
gap_thresh <- 0.25 * max_height_global
r_gap_bin <- grid_metrics(las_nonground, ~as.numeric(length(Z[Z > gap_thresh]) == 0), res = 1)

# Cover: inverse of deep gap
r_cover_bin <- calc(r_deep_gap_bin, fun = function(x) ifelse(is.na(x), NA, 1 - x))

# Set CRS and save rasters
crs_code <- "EPSG:5514"
crs(r_cover_bin) <- crs_code
crs(r_deep_gap_bin) <- crs_code
crs(r_gap_bin) <- crs_code

writeRaster(r_cover_bin,      file.path(raster_output_path, paste0(base_name, "_CoverFraction.tif")),     overwrite = TRUE)
writeRaster(r_deep_gap_bin,   file.path(raster_output_path, paste0(base_name, "_DeepGapFraction.tif")),   overwrite = TRUE)
writeRaster(r_gap_bin,        file.path(raster_output_path, paste0(base_name, "_GapFraction.tif")),       overwrite = TRUE)

# Summary metrics (as % of area)
cover_vals <- values(r_cover_bin)
deep_gap_vals <- values(r_deep_gap_bin)
gap_vals <- values(r_gap_bin)

metrics <- data.frame(
  CoverFraction_percent    = round(100 * sum(cover_vals == 1, na.rm = TRUE) / sum(!is.na(cover_vals)), 2),
  DeepGapFraction_percent  = round(100 * sum(deep_gap_vals == 1, na.rm = TRUE) / sum(!is.na(deep_gap_vals)), 2),
  GapFraction_percent      = round(100 * sum(gap_vals == 1, na.rm = TRUE) / sum(!is.na(gap_vals)), 2)
)

write.table(metrics, file = file.path(output_folder, paste0(base_name, "_gappines_metrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "        ", quote = FALSE)

# Save PNGs
png(file.path(output_folder, paste0(base_name, "_CoverFraction.png")), width = 800, height = 600)
plot(r_cover_bin, main = "Canopy Cover Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
dev.off()

png(file.path(output_folder, paste0(base_name, "_DeepGapFraction.png")), width = 800, height = 600)
plot(r_deep_gap_bin, main = "Deep Gap Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
dev.off()

png(file.path(output_folder, paste0(base_name, "_GapFraction.png")), width = 800, height = 600)
plot(r_gap_bin, main = "Gap Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
dev.off()

cat("Gap and cover metrics processing complete. Results saved to:", output_folder, "\n")
