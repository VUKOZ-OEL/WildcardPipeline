# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript HeightVariabilityMetrics.R <input_folder> <output_folder> <output_prefix> <max_height_file>")
}

input_folder <- args[1]
output_folder <- args[2]
base_name <- args[3]
max_height_file <- args[4]

# Read global max height
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global)) stop("Invalid max height value in the file.")

# Create output directory
raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)

# Load LAS files
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
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

# Compute metrics
cat("Metric 1\n")
r1 <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
crs(r1) <- "EPSG:5514"
writeRaster(r1, file.path(raster_output_path, paste0(base_name, "_TopCanopyHeight.tif")), overwrite = TRUE)

cat("Metric 2\n")
r2 <- grid_metrics(las_nonground, ~mean(Z, na.rm = TRUE), res = 1)
crs(r2) <- "EPSG:5514"
writeRaster(r2, file.path(raster_output_path, paste0(base_name, "_MeanCanopyHeight.tif")), overwrite = TRUE)

cat("Metric 3\n")
cv_metric <- function(z) {
  m <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA) else return(s / m)
}
r3 <- grid_metrics(las_nonground, ~cv_metric(Z), res = 1)

cat("Metric 4\n")
voxel_fun <- function(z, dz = 1) {
  z_bin <- unique(floor(z / dz))
  if (length(z_bin) == 0) return(NA)
  return(max(z_bin) - min(z_bin) + 1)
}
r4 <- grid_metrics(las_nonground, voxel_fun(Z), res = 1)
crs(r4) <- "EPSG:5514"
writeRaster(r4, file.path(raster_output_path, paste0(base_name, "_CanopyHeightVoxelRange.tif")), overwrite = TRUE)

# Compute rumple index from surface
cat("Metric 5\n")
surface_raster <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
rumple <- rumple_index(surface_raster)

cat("Metric 6\n")
# Site-level max canopy height
max_canopy_height_site <- quantile(las_nonground$Z, 0.98, na.rm = TRUE)

cat("Saving metrics to text file\n")
# Save summary metrics
metrics <- data.frame(
  MaxCanopyHeight = max_canopy_height_site,
  RumpleIndex = rumple,
  AvgTopCanopyHeight = mean(values(r1), na.rm = TRUE),
  AvgMeanCanopyHeight = mean(values(r2), na.rm = TRUE),
  AvgCanopyHeightCoV = mean(values(r3), na.rm = TRUE),
  AvgCanopyHeightVoxelRange = mean(values(r4), na.rm = TRUE)
)

cat("Write table\n")
write.table(metrics, file = file.path(output_folder, paste0(base_name, "_HeightVariabilityMetrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# PNGs
cat("Make pngs\n")
png(file.path(output_folder, paste0(base_name, "_TopCanopyHeight.png")), width = 800, height = 600)
plot(r1, main = "Top Canopy Height", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_MeanCanopyHeight.png")), width = 800, height = 600)
plot(r2, main = "Mean Canopy Height", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_CanopyHeightCoV.png")), width = 800, height = 600)
plot(r3, main = "Canopy Height Coefficient of Variation", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_CanopyHeightVoxelRange.png")), width = 800, height = 600)
plot(r4, main = "Canopy Height Voxel Range", xlab = "X", ylab = "Y")
dev.off()

cat("Height variability metrics processing complete. Results saved to:", output_folder, "\n")
