# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(sf))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript CanopyDensityMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix> <LOCAL|EPSG:code>")
}

input_folder <- args[1]
output_folder <- args[2]
max_height_file <- args[3]
base_name <- args[4]
coordinates <- trimws(args[5])

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

# Filter to non-ground points above 1m
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

cat("Canopy density metrics calculation \n")
cat("Metric 1 calculation \n")
# Metric 1: Canopy Density (points per column)
# Canopy density = number of nonground points above 1 m per m²
r_density <- grid_density(las_nonground, res = 1)
r_density <- set_output_crs(r_density, las_crs_obj)
writeRaster(r_density, file.path(raster_output_path, paste0(base_name, "_CanopyDensity.tif")), overwrite = TRUE)

cat("Metric 3 calculation \n")
voxels <- voxel_metrics(las_nonground, ~list(density = length(Z)), res = 1)
voxel_df <- as.data.frame(voxels)
column_rugosity <- voxel_df %>% group_by(X, Y) %>% summarise(rugosity = sd(density, na.rm = TRUE), .groups = "drop")
rugosity_raster <- rasterFromXYZ(column_rugosity)
rugosity_raster <- set_output_crs(rugosity_raster, las_crs_obj)
writeRaster(rugosity_raster, file.path(raster_output_path, paste0(base_name, "_CanopyRugosity.tif")), overwrite = TRUE)
global_rugosity <- mean(column_rugosity$rugosity, na.rm = TRUE)

cat("Metric 4 calculation \n")
# Metric 4: Clumping Index
calculate_clumping_index <- function(chm, max_height_global, height_bins = seq(0, max_height_global, by = 1)) {
  gap_fractions <- c()
  log_gap_fractions <- c()
  for (h in height_bins) {
    gap_cells <- sum(values(chm) < h, na.rm = TRUE)
    total_cells <- sum(!is.na(values(chm)))
    gap_fraction <- ifelse(total_cells > 0, gap_cells / total_cells, NA)
    if (!is.na(gap_fraction) && gap_fraction > 0) {
      gap_fractions <- c(gap_fractions, gap_fraction)
      log_gap_fractions <- c(log_gap_fractions, log(gap_fraction))
    }
  }
  if (length(gap_fractions) > 0 && length(log_gap_fractions) > 0) {
    avg_gap_fraction <- mean(gap_fractions, na.rm = TRUE)
    avg_log_gap_fraction <- mean(log_gap_fractions, na.rm = TRUE)
    return(log(avg_gap_fraction) / avg_log_gap_fraction)
  } else {
    return(NA)
  }
}

chm <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
clumping_val <- calculate_clumping_index(chm, max_height_global)

cat("Save metric txt \n")
# Save summary metrics
metrics <- data.frame(
  AvgClumpingIndex = clumping_val,
  AvgCanopyDensity = mean(values(r_density), na.rm = TRUE),
  CanopyRugosity = global_rugosity
)

write.table(metrics, file = file.path(output_folder, paste0(base_name, "_CanopyDensityMetrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

cat("Creating PNGs \n")
# PNGs
png(file.path(output_folder, paste0(base_name, "_CanopyDensity.png")), width = 800, height = 600)
plot(r_density, main = "Canopy Density", xlab = "X", ylab = "Y")
dev.off()

#png(file.path(output_folder, paste0(base_name, "_CanopyPorosity.png")), width = 800, height = 600)
#plot(r_porosity, main = "Canopy Porosity", xlab = "X", ylab = "Y")
#dev.off()

png(file.path(output_folder, paste0(base_name, "_CanopyRugosity.png")), width = 800, height = 600)
plot(rugosity_raster, main = "Canopy Rugosity", xlab = "X", ylab = "Y")
dev.off()

cat("Canopy density metrics processing complete. Results saved to:", output_folder, "\n")
