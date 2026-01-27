suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))

# ----------------------------
# Args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop("Usage: Rscript CanopyVolumeMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix> <nOfElements> <coordinates>")
}

las_folder      <- args[1]
output_folder   <- args[2]
max_height_file <- args[3]
base_name       <- args[4]
nOfElements     <- args[5]
coordinates     <- args[6]

cat("nOfElements (point-count threshold per voxel) =", nOfElements, "\n")

if (!dir.exists(las_folder)) stop("input_folder not found: ", las_folder)
if (!file.exists(max_height_file)) stop("max_height_file not found: ", max_height_file)

dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Read global max height
# ----------------------------
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global) || max_height_global <= 0) stop("Invalid max height value in the file.")
total_layers_global <- ceiling(max_height_global)

# ----------------------------
# Load and merge LAS/LAZ
# ----------------------------
las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS/LAZ files found in: ", las_folder)

cat("Reading and merging LAS files...\n")
las_list <- lapply(las_files, function(f) {
  tryCatch({
    las <- readLAS(f)
    if (!is.null(las) && !is.empty(las)) las else NULL
  }, error = function(e) NULL)
})
las_list <- Filter(Negate(is.null), las_list)
if (length(las_list) == 0) stop("No valid point cloud data after reading LAS/LAZ files.")

las_combined <- do.call(rbind, las_list)
if (is.null(las_combined) || is.empty(las_combined)) stop("No valid point cloud data after merging.")

# ----------------------------
# Filter: non-ground, above 1 m
# ----------------------------
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)
if (is.null(las_nonground) || is.empty(las_nonground)) stop("No non-ground points above 1 m after filtering.")

rm(las_combined, las_list)
gc(FALSE)

# ----------------------------
# Voxel metrics
#   - V1 = number of points in each 1m voxel
#   - all_voxels=TRUE keeps empty voxels (V1 will be NA for those)
# ----------------------------
cat("Computing voxel point counts...\n")
voxels <- voxel_metrics(las_nonground, ~length(Z), res = 1)#, all_voxels=TRUE)
voxel_dt <- as.data.table(voxels)

# Safety: keep only rows with defined coordinates
voxel_dt <- voxel_dt[is.finite(X) & is.finite(Y) & is.finite(Z)]

# "Full/filled voxel" definition: point-count threshold
full_voxels <- voxel_dt[!is.na(V1) & V1 >= nOfElements]

cat("Voxel V1 summary (points per 1m voxel):\n")
print(summary(voxel_dt$V1))

cat("Example voxels (first 20 rows): X Y Z V1\n")
print(voxel_dt[, .(X, Y, Z, V1)][1:min(20, .N)])

# If nothing meets the threshold, write empty outputs but do not crash
if (nrow(full_voxels) == 0) {
  warning("No voxels met the threshold (V1 >= nOfElements). Outputs will be empty/NA.")
}

# ----------------------------
# Occupancy per (X,Y) column:
# FullVoxelCount = number of Z-voxels (1m layers) with >= nOfElements points
# ----------------------------
occupancy <- full_voxels[, .(FullVoxelCount = .N), by = .(X, Y)]

# Relative occupancy vs global maximum height (in 1m layers)
occupancy[, RelativeTotalOccupancy := pmin((FullVoxelCount / total_layers_global) * 100, 100)]

# Local maximum height proxy per column based on voxel grid:
# Use max voxel Z center; convert to approximate number of 1m layers via ceiling().
max_local <- voxel_dt[, .(MaxLocalHeight = suppressWarnings(max(Z, na.rm = TRUE))), by = .(X, Y)]
occupancy <- merge(occupancy, max_local, by = c("X", "Y"), all.x = TRUE)

occupancy[, LocalLayers := {
  h <- MaxLocalHeight
  ifelse(is.finite(h) & h > 0, ceiling(h), NA_integer_)
}]

occupancy[, RelativeLocalOccupancy := {
  denom <- LocalLayers
  ifelse(!is.na(denom) & denom > 0, pmin((FullVoxelCount / denom) * 100, 100), NA_real_)
}]

# ----------------------------
# Output rasters
# ----------------------------
cat("Writing rasters and PNGs...\n")
crs_code <- coordinates
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

r1 <- rasterFromXYZ(occupancy[, .(X, Y, FullVoxelCount)])
crs(r1) <- crs_code
writeRaster(r1, file.path(output_folder, paste0(base_name, "_CanopyOccupancy.tif")), overwrite = TRUE)

#r2 <- rasterFromXYZ(occupancy[, .(X, Y, RelativeTotalOccupancy)])
#crs(r2) <- crs_code
#writeRaster(r2, file.path(output_folder, paste0(base_name, "_RelativeTotalOccupancy.tif")), overwrite = TRUE)

r3 <- rasterFromXYZ(occupancy[, .(X, Y, RelativeLocalOccupancy)])
crs(r3) <- crs_code
writeRaster(r3, file.path(output_folder, paste0(base_name, "_RelativeLocalOccupancy.tif")), overwrite = TRUE)

png(file.path(output_folder, paste0(base_name, "_CanopyOccupancy.png")), width = 800, height = 600)
plot(r1, main = "Canopy Occupancy (Full voxels per 1m column)", xlab = "X", ylab = "Y")
dev.off()

#png(file.path(output_folder, paste0(base_name, "_RelativeTotalOccupancy.png")), width = 800, height = 600)
#plot(r2, main = "Relative Total Occupancy (%)", xlab = "X", ylab = "Y", zlim = c(0, 100))
#dev.off()

png(file.path(output_folder, paste0(base_name, "_RelativeLocalOccupancy.png")), width = 800, height = 600)
plot(r3, main = "Relative Local Occupancy (%)", xlab = "X", ylab = "Y", zlim = c(0, 100))
dev.off()

# ----------------------------
# Summary metrics
# ----------------------------
summary_metrics <- data.frame(
  AvgCanopyOccupancy        = mean(occupancy$FullVoxelCount, na.rm = TRUE),
  #AvgRelativeTotalOccupancy = mean(occupancy$RelativeTotalOccupancy, na.rm = TRUE),
  AvgRelativeLocalOccupancy = mean(occupancy$RelativeLocalOccupancy, na.rm = TRUE)
)

write.table(
  summary_metrics,
  file = file.path(output_folder, paste0(base_name, "_CanopyVolumeMetrics.txt")),
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  quote = FALSE
)

cat("Canopy volume metrics processing complete. Results saved to:", output_folder, "\n")
