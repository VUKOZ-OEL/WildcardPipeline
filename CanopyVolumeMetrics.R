# Load required packages
suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript CanopyVolumeMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix> <voxel_size>")
}

las_folder <- args[1]
output_folder <- args[2]
max_height_file <- args[3]
base_name <- args[4]
voxel_size <- args[5]

# Read global max height
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global)) stop("Invalid max height value in the file.")

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

# Filter to non-ground points above 1 meters
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

# Compute voxel metrics on non-ground points
voxels <- voxel_metrics(las_nonground, ~length(Z), res = 1)
voxel_dt <- as.data.table(voxels)

# Parse voxel size
voxel_size_numeric <- as.numeric(voxel_size)
# Align to 1m² columns and bin Z axis to 1m layers
voxel_dt[, colX := floor(X)]
voxel_dt[, colY := floor(Y)]
voxel_dt[, colZbin := floor(Z)]

# Count unique filled voxels per (X,Y,Z)
filled_voxels <- voxel_dt[V1 >= 10, .N, by = .(colX, colY, colZbin)]

# Now count number of Z-layers per (X,Y)
local_density_per_column <- filled_voxels[, .N, by = .(colX, colY)]
setnames(local_density_per_column, "N", "FilledVoxels1m3")

# Compute occupancy metrics
full_voxels <- voxel_dt[V1 >= 10, ]
occupancy <- full_voxels[, .N, by = .(X, Y)]
setnames(occupancy, "N", "FullVoxelCount")

occupancy[, RelativeTotalOccupancy := (FullVoxelCount / max_height_global) * 100]
max_local <- voxel_dt[, .(MaxLocalHeight = max(Z, na.rm = TRUE)), by = .(X, Y)]
occupancy <- merge(occupancy, max_local, by = c("X", "Y"), all.x = TRUE)
occupancy[, RelativeLocalOccupancy := pmin((FullVoxelCount / MaxLocalHeight) * 100, 100)]

# Create and save rasters
crs_code <- "EPSG:5514"
raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)

r1 <- rasterFromXYZ(occupancy[, .(X, Y, FullVoxelCount)])
crs(r1) <- crs_code
writeRaster(r1, file.path(raster_output_path, paste0(base_name, "_CanopyOccupancy.tif")), overwrite = TRUE)

r2 <- rasterFromXYZ(occupancy[, .(X, Y, RelativeTotalOccupancy)])
crs(r2) <- crs_code
writeRaster(r2, file.path(raster_output_path, paste0(base_name, "_RelativeTotalOccupancy.tif")), overwrite = TRUE)

r3 <- rasterFromXYZ(occupancy[, .(X, Y, RelativeLocalOccupancy)])
crs(r3) <- crs_code
writeRaster(r3, file.path(raster_output_path, paste0(base_name, "_RelativeLocalOccupancy.tif")), overwrite = TRUE)

# Save PNGs for rasters
png(file.path(output_folder, paste0(base_name, "_CanopyOccupancy.png")), width = 800, height = 600)
plot(r1, main = "Canopy Occupancy (Filled Voxels per Column)", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_RelativeTotalOccupancy.png")), width = 800, height = 600)
plot(r2, main = "Relative Total Occupancy (%)", xlab = "X", ylab = "Y", zlim = c(0, 100))
dev.off()

png(file.path(output_folder, paste0(base_name, "_RelativeLocalOccupancy.png")), width = 800, height = 600)
plot(r3, main = "Relative Local Occupancy (%)", xlab = "X", ylab = "Y", zlim = c(0, 100))
dev.off()

# Save summary metrics
write.table(data.frame(
  AvgCanopyOccupancy = mean(occupancy$FullVoxelCount, na.rm = TRUE),
  AvgRelativeTotalOccupancy = mean(occupancy$RelativeTotalOccupancy, na.rm = TRUE),
  AvgRelativeLocalOccupancy = mean(occupancy$RelativeLocalOccupancy, na.rm = TRUE)
), file = file.path(output_folder, paste0(base_name, "_CanopyVolumeMetrics.txt")),
row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

cat("Canopy volume metrics processing complete. Results saved to:", output_folder, "\n")
