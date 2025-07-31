# Load required libraries
library(lidR)
library(raster)
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript preprocess_single_las.R <las_folder_path> <processed_folder_path> <results_folder_path>")
}

las_folder_path       <- args[1]
processed_folder_path <- args[2]
results_folder_path   <- args[3]

dir.create(processed_folder_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_folder_path, showWarnings = FALSE, recursive = TRUE)

# List LAS/LAZ files and take the first one
las_files <- list.files(las_folder_path, pattern = "\\.(las|laz)$", full.names = TRUE)
if (length(las_files) == 0) {
  stop("No LAS or LAZ files found in the specified folder: ", las_folder_path)
}
las_file <- las_files[1]
cat("Processing file:", las_file, "\n")

# Read the LAS file safely
las <- tryCatch(readLAS(las_file), error = function(e) NULL)

if (is.null(las) || npoints(las) == 0) {
  stop("Failed to read LAS file or file contains no points: ", las_file)
}

# Step 1: Generate DTM
cat("Step 1: Building DTM...\n")
# Classify ground points
# Step 1: Store original ground points
original_ground <- filter_poi(las, Classification == 2)
n_original <- npoints(original_ground)
cat("Original ground points:", n_original, "\n")

# Step 2: Try to reclassify ground
las_reclassified <- classify_ground(las, algorithm = csf())
new_ground <- filter_poi(las_reclassified, Classification == 2)
n_new <- npoints(new_ground)
cat("Newly classified ground points:", n_new, "\n")

# Step 3: Decide which to use
if (n_new > n_original) {
  cat("Using newly classified ground points.\n")
  las_ground <- las_reclassified
}else {
  cat("Using old classification of ground points")
  las_ground <- original_ground
}

ground_points <- filter_poi(las_ground, Classification == 2)
cat("Final ground point count:", npoints(ground_points), "\n")

# Safety check
if (npoints(ground_points) < 3) {
  stop("Too few ground points after filtering.")
}

# Step 4: Try DTM generation with knnidw; fallback to tin if needed, third is failsafe if the plot was already preprocessed
cat("Generating DTM using knnidw...\n")
dtm <- tryCatch({
  cat("Generating DTM using knnidw...\n")
  grid_terrain(ground_points, res = 1, algorithm = knnidw(k = 10, p = 2))
}, error = function(e_knn) {
  cat("knnidw failed:", conditionMessage(e_knn), "\n")
  cat("Falling back to TIN interpolation...\n")

  # Deduplicate
  xy_keys <- paste(ground_points@data$X, ground_points@data$Y)
  keep_idx <- !duplicated(xy_keys)
  ground_points_dedup <- filter_poi(ground_points, keep_idx)
  cat("Ground points after removing duplicates:", npoints(ground_points_dedup), "\n")

  # Try TIN or fallback to flat DTM
  tryCatch({
    if (npoints(ground_points_dedup) < 3) {
      stop("Too few unique ground points for TIN interpolation.")
    }
    grid_terrain(ground_points_dedup, res = 1, algorithm = tin())
  }, error = function(e_tin) {
    cat("TIN also failed:", conditionMessage(e_tin), "\n")
  })
})

print(dtm)
print(crs(dtm))

# Write the DTM
cat("Creating dtm raster...")
dtm_tif <- file.path(results_folder_path, "DTM.tif")
raster::writeRaster(dtm, dtm_tif, overwrite = TRUE)
cat("DTM saved to", dtm_tif, "\n")

# Save a PNG preview of the DTM
png(file.path(results_folder_path, "DTM.png"), width = 1000, height = 800)
plot(dtm, main = "Digital Terrain Model (DTM)")
dev.off()

# Step 2: Normalize heights using the DTM
cat("Step 2: Normalizing point cloud...\n")
las_norm <- normalize_height(las, dtm)
# Check for valid normalized cloud
if (npoints(las_norm) == 0) {
  stop("Normalization resulted in empty point cloud for: ", las_file)
}

# Clamp negative Z values to zero
las_norm$Z[las_norm$Z < 0] <- 0

# Write the normalized LAS to the processed folder
out_las <- file.path(processed_folder_path, basename(las_file))
writeLAS(las_norm, out_las)
cat("Normalized LAS saved to", out_las, "\n")

# Step 3: Generate CHM
cat("Step 3: Building CHM...\n")
chm <- grid_canopy(las_norm, res = 1, algorithm = p2r(subcircle = 0.2))

# Write the CHM to disk
chm_tif <- file.path(results_folder_path, "CHM.tif")
raster::writeRaster(chm, chm_tif, overwrite = TRUE)
cat("CHM saved to", chm_tif, "\n")
png(file.path(results_folder_path, "CHM.png"), width = 1000, height = 800)
plot(chm, main = "Canopy Height Model (CHM)")
dev.off()

cat("Processing complete.\n")
cat("Processed LAS file saved to:", out_las, "\n")
cat("DTM and CHM outputs saved to:", results_folder_path, "\n") 