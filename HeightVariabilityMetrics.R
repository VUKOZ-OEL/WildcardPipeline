# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(sf))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
cat("Amount of arguments is: ", length(args), "\n")
cat("Arguments are: ", args[1], args[2], args[3], args[4])
if (length(args) != 6) {
  stop("Usage: Rscript HeightVariabilityMetrics.R <input_folder> <output_folder> <output_prefix> <max_height_file> <nOfElements> <LOCAL|EPSG:code>")
}

input_folder <- args[1]
output_folder <- args[2]
base_name <- args[3]
max_height_file <- args[4]
nOfElements <- args[5]
coordinates <- trimws(args[6])

if (is.na(coordinates) || !nzchar(coordinates)) {
  stop("coordinates must be LOCAL or a valid CRS such as EPSG:3035.")
}

local_mode <- identical(toupper(coordinates), "LOCAL")
if (local_mode) {
  coordinates <- "LOCAL"
  configured_crs <- sf::NA_crs_
} else {
  configured_crs <- suppressWarnings(sf::st_crs(coordinates))
  if (is.na(configured_crs)) {
    stop("Invalid coordinate system. Use LOCAL or a valid CRS such as EPSG:3035; got: ", coordinates)
  }
}

crs_missing <- function(x) is.null(x) || is.na(x)
crs_text_for_raster <- function(x) {
  if (crs_missing(x)) return(coordinates)
  input <- tryCatch(x$input, error = function(e) NA_character_)
  wkt <- tryCatch(x$wkt, error = function(e) NA_character_)
  if (!is.na(input) && nzchar(input)) return(input)
  if (!is.na(wkt) && nzchar(wkt)) return(wkt)
  coordinates
}
set_output_crs <- function(r, source_crs) {
  if (local_mode) {
    raster::crs(r) <- NA_character_
  } else {
    raster::crs(r) <- crs_text_for_raster(source_crs)
  }
  r
}

if (!dir.exists(input_folder))  stop("input_folder not found: ", input_folder)
#if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
if (!file.exists(max_height_file)) stop("max_height_file not found: ", max_height_file)

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

las_crs_obj <- tryCatch(lidR::st_crs(las_combined), error = function(e) NULL)
if (local_mode) {
  lidR::st_crs(las_combined) <- sf::NA_crs_
  las_crs_obj <- sf::NA_crs_
  cat("Coordinate mode: LOCAL; output rasters will have no CRS.\n")
} else if (crs_missing(las_crs_obj)) {
  lidR::st_crs(las_combined) <- configured_crs
  las_crs_obj <- lidR::st_crs(las_combined)
  cat("LAS had no CRS; using configured fallback: ", coordinates, "\n", sep = "")
}

# Filter to non-ground points above 1 meters
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

# Compute metrics
cat("Metric 1\n")
r1 <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
r1 <- set_output_crs(r1, las_crs_obj)
writeRaster(r1, file.path(raster_output_path, paste0(base_name, "_TopCanopyHeight.tif")), overwrite = TRUE)

# ---------- TopCanopyHeight value histogram ----------
png(
  file.path(output_folder, paste0(base_name, "_TopCanopyHeight_hist.png")),
  width  = 1000,
  height = 800
)
raster::hist(
  r1,
  main = "Top Canopy Height Histogram",
  xlab  = "Height (m)",
  ylab  = "Frequency"
)
dev.off()
# -------------------------------------------


cat("Metric 2\n")
r2 <- grid_metrics(las_nonground, ~mean(Z, na.rm = TRUE), res = 1)
r2 <- set_output_crs(r2, las_crs_obj)
writeRaster(r2, file.path(raster_output_path, paste0(base_name, "_MeanCanopyHeight.tif")), overwrite = TRUE)

cat("Metric 3\n")
cv_metric <- function(z) {
  m <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  if (is.na(m) || m == 0) return(NA) else return(s / m)
}
r3 <- grid_metrics(las_nonground, ~cv_metric(Z), res = 1)
r3 <- set_output_crs(r3, las_crs_obj)
writeRaster(r3, file.path(raster_output_path, paste0(base_name, "_CanopyHeightCoV.tif")), overwrite = TRUE)

#cat("Metric 4\n")
#voxel_fun <- function(z, dz = 1, nOfElements = 10) {
#  z <- z[is.finite(z)]
#  if (length(z) == 0) return(NA_real_)
#  bins <- floor(z / dz)
#  bin_counts <- table(bins)
#  valid_bins <- as.integer(names(bin_counts)[bin_counts >= nOfElements])
#  if (length(valid_bins) == 0) return(NA_real_)
#  max(valid_bins) - min(valid_bins) + 1
#}
#r4 <- grid_metrics(las_nonground, ~voxel_fun(Z, dz = 1, nOfElements = nOfElements), res = 1)
#r4 <- set_output_crs(r4, las_crs_obj)
#writeRaster(r4, file.path(raster_output_path, paste0(base_name, "_CanopyHeightVoxelRange.tif")), overwrite = TRUE)

# Compute rumple index from surface
cat("Metric 5\n")
chm <- grid_canopy(las_nonground, res = 1, algorithm = p2r(subcircle = 0.2))
rumple <- rumple_index(chm, na.rm=TRUE)

cat("Metric 6\n")
# Site-level max canopy height
max_canopy_height_site <- quantile(las_nonground$Z, 0.98, na.rm = TRUE)

cat("Metric 7\n")
# Site-level top canopy rugosity
top_canopy_rugosity <- sd(values(r1), na.rm = TRUE)

cat("Saving metrics to text file\n")
# Save summary metrics
metrics <- data.frame(
  MaxCanopyHeight = max_height_global,
  thPercentileMaxHeight = max_canopy_height_site,
  RumpleIndex = rumple,
  TopCanopyRugosity = top_canopy_rugosity,
  AvgTopCanopyHeight = mean(values(r1), na.rm = TRUE),
  AvgMeanCanopyHeight = mean(values(r2), na.rm = TRUE),
  AvgCanopyHeightCoV = mean(values(r3), na.rm = TRUE)
  #AvgCanopyHeightVoxelRange = mean(values(r4), na.rm = TRUE)
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

#png(file.path(output_folder, paste0(base_name, "_CanopyHeightVoxelRange.png")), width = 800, height = 600)
#plot(r4, main = "Canopy Height Voxel Range", xlab = "X", ylab = "Y")
#dev.off()

cat("Height variability metrics processing complete. Results saved to:", output_folder, "\n")
