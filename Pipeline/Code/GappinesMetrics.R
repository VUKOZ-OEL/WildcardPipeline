suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(sf))
suppressMessages(library(terra))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("Usage: Rscript GappinesMetrics.R <input_folder> <output_folder> <max_height_file> <output_prefix> <polygon_file> <numberOfElements> <coordinates>")
}

las_folder       <- args[1]
output_folder    <- args[2]
max_height_file  <- args[3]
base_name        <- args[4]
polygon_file     <- args[5]
nOfElements      <- args[6]
coordinates      <- args[7]

cat("Output folder: ", output_folder, "\n")
cat("Polygon file: ", polygon_file, "\n")
cat("Number of elements: ", nOfElements, "\n")
cat("Coordinate system: ", coordinates, "\n")

# Read global max height
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global)) stop("Invalid max height value in the file.")

# Output raster folder
raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)

crs_code <- coordinates
cat("Loading las file... \n")
# --- Load LAS files --------------------------------------------------------- #
las_files <- list.files(las_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS/LAZ files found.")
  
las_crs <- crs_code
if (is.na(las_crs) || las_crs == "") {
  stop("LAS has no CRS in its header. Set it explicitly before computing rasters.")
}

cat("Loading polygon... \n")
# --- Read & prep polygon (sf) ---------------------------------------------- #
poly <- st_read(polygon_file, quiet = TRUE)
if (is.null(st_crs(poly))) stop("Polygon has no CRS. Please define one.")
poly_crs <- st_crs(poly)
if (is.na(poly_crs)) stop("Polygon has no CRS (.prj missing/invalid).")

cat("Reproject polygon to las crs \n")
# Reproject polygon to LAS CRS
poly <- st_transform(poly, crs = las_crs)
poly <- st_make_valid(poly)
poly <- suppressWarnings(st_union(poly))  # dissolve to one geometry

cat("Reading and merging LAS files \n")
# Now actually read & merge all LAS tiles
las_list <- lapply(las_files, function(file) {
  las <- tryCatch(readLAS(file), error = function(e) NULL)
  if (is.null(las) || npoints(las) == 0) return(NULL)
  return(las)
})

las_combined <- do.call(rbind, las_list)
if (is.null(las_combined) || npoints(las_combined) == 0) stop("No points in LAS files.")

# --- Non-ground & height filter -------------------------------------------- #
las_nonground <- filter_poi(las_combined, Classification != 2 & Z > 1)

# --- Metrics over full LAS extent ------------------------------------------ #
deep_gap_metric <- function(z) {
  n_above_1 <- length(z[z > 1])
  if (is.na(n_above_1) || n_above_1 < nOfElements) return(1) else return(0)
}

# Base rasters (1 m) produced from the NON-CLIPPED LAS (full tiles)
r_deep_gap_bin <- grid_metrics(las_nonground, ~deep_gap_metric(Z), res = 1)

gap_thresh <- 0.25 * max_height_global
r_gap_bin  <- grid_metrics(las_nonground, ~as.numeric(length(Z[Z > gap_thresh]) < nOfElements), res = 1)

cat("Cliping rasters with polygons")
# --- Clip rasters with polygon using terra --------------------------------- #

# Convert raster::RasterLayer -> terra::SpatRaster
s_deep_gap_bin <- rast(r_deep_gap_bin)
s_gap_bin      <- rast(r_gap_bin)

# Ensure polygon is a terra vector in the same CRS
poly_v <- vect(poly)  # terra will read CRS from sf object

# Treat NA (no returns) as gaps EVERYWHERE first
s_deep_gap_bin_fill <- ifel(is.na(s_deep_gap_bin), 1, s_deep_gap_bin)
s_gap_bin_fill      <- ifel(is.na(s_gap_bin),      1, s_gap_bin)

# Now mask with polygon: outside polygon -> NA, inside keeps 0/1 (including filled 1)
s_deep_gap_bin_clip <- mask(s_deep_gap_bin_fill, poly_v)
s_gap_bin_clip      <- mask(s_gap_bin_fill,      poly_v)

# Cover = inverse of deep gap inside polygon; keep NA outside polygon
#s_cover_bin <- ifel(is.na(s_deep_gap_bin_clip), NA, 1 - s_deep_gap_bin_clip)

# Convert back to RasterLayer for the rest of your pipeline
r_deep_gap_bin <- raster(s_deep_gap_bin_clip)
r_gap_bin      <- raster(s_gap_bin_clip)
#r_cover_bin    <- raster(s_cover_bin)

# --- Save rasters ----------------------------------------------------------- #
#writeRaster(r_cover_bin,    file.path(raster_output_path, paste0(base_name, "_CoverFraction.tif")),   overwrite = TRUE)
writeRaster(r_deep_gap_bin, file.path(raster_output_path, paste0(base_name, "_DeepGapFraction.tif")), overwrite = TRUE)
writeRaster(r_gap_bin,      file.path(raster_output_path, paste0(base_name, "_GapFraction.tif")),     overwrite = TRUE)

# --- Summary metrics (percent of polygon area only) ------------------------ #
#cover_vals    <- values(r_cover_bin)
deep_gap_vals <- values(r_deep_gap_bin)
gap_vals      <- values(r_gap_bin)

metrics <- data.frame(
  #CoverFraction_percent    = round(100 * sum(cover_vals == 1,    na.rm = TRUE) / sum(!is.na(cover_vals)),    2),
  DeepGapFraction_percent  = round(100 * sum(deep_gap_vals == 1, na.rm = TRUE) / sum(!is.na(deep_gap_vals)), 2),
  GapFraction_percent      = round(100 * sum(gap_vals == 1,      na.rm = TRUE) / sum(!is.na(gap_vals)),      2)
)

write.table(metrics, file = file.path(output_folder, paste0(base_name, "_gappines_metrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "        ", quote = FALSE)

# --- PNGs ------------------------------------------------------------------ #
#png(file.path(output_folder, paste0(base_name, "_CoverFraction.png")), width = 800, height = 600)
#plot(r_cover_bin, main = "Canopy Cover Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
#dev.off()

png(file.path(output_folder, paste0(base_name, "_DeepGapFraction.png")), width = 800, height = 600)
plot(r_deep_gap_bin, main = "Deep Gap Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
dev.off()

png(file.path(output_folder, paste0(base_name, "_GapFraction.png")), width = 800, height = 600)
plot(r_gap_bin, main = "Gap Fraction (Binary)", xlab = "X", ylab = "Y", zlim = c(0, 1))
dev.off()

cat("Gap and cover metrics processing complete. Results saved to:", output_folder, "\n")
