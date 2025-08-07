# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript CanopyDensityMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix>")
}

input_folder <- args[1]
output_folder <- args[2]
max_height_file <- args[3]
base_name <- args[4]

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

# Filter to non-ground points above 1m
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

cat("Metric 1 calculation \n")
# Metric 1: Canopy Density (points per column)
r_density <- grid_density(las_nonground, res = 1)
crs(r_density) <- "EPSG:5514"
writeRaster(r_density, file.path(raster_output_path, paste0(base_name, "_CanopyDensity.tif")), overwrite = TRUE)

cat("Metric 2 calculation \n")
# Metric 2: Canopy Porosity (ratio of ground to all points in column, but from nonground points context)
las_all <- filter_poi(las_combined, Z > 2)  # Above 2m points (both ground & non-ground)
porosity_metric <- function(z, classification) {
  total <- length(z)
  ground <- sum(classification == 2, na.rm = TRUE)
  if (total == 0) return(NA) else return(ground / total)
}
r_porosity <- grid_metrics(las_all, porosity_metric(Z, Classification), res = 1)
crs(r_porosity) <- "EPSG:5514"
writeRaster(r_porosity, file.path(raster_output_path, paste0(base_name, "_CanopyPorosity.tif")), overwrite = TRUE)

cat("Metric 3 calculation \n")
# Metric 3: Top Canopy Rugosity = SD of per-column 98th percentile height
r_top_height <- grid_metrics(las_nonground, ~quantile(Z, 0.98, na.rm = TRUE), res = 1)
crs(r_top_height) <- "EPSG:5514"
writeRaster(r_top_height, file.path(raster_output_path, paste0(base_name, "_TopCanopyHeight.tif")), overwrite = TRUE)
top_rugosity <- sd(values(r_top_height), na.rm = TRUE)

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

chm <- grid_metrics(las_nonground, ~mean(Z, na.rm = TRUE), res = 1)
clumping_val <- calculate_clumping_index(chm, max_height_global)

cat("Save metric txt \n")
# Save summary metrics
metrics <- data.frame(
  AvgClumpingIndex = clumping_val,
  AvgCanopyPorosity = mean(values(r_porosity), na.rm = TRUE),
  AvgCanopyDensity = mean(values(r_density), na.rm = TRUE),
  TopCanopyRugosity = top_rugosity
)

write.table(metrics, file = file.path(output_folder, paste0(base_name, "_CanopyDensityMetrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

cat("Creating PNGs \n")
# PNGs
png(file.path(output_folder, paste0(base_name, "_CanopyDensity.png")), width = 800, height = 600)
plot(r_density, main = "Canopy Density", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_CanopyPorosity.png")), width = 800, height = 600)
plot(r_porosity, main = "Canopy Porosity", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_TopCanopyHeight.png")), width = 800, height = 600)
plot(r_top_height, main = "Top Canopy Height", xlab = "X", ylab = "Y")
dev.off()

cat("Canopy density metrics processing complete. Results saved to:", output_folder, "\n")
