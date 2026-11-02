suppressMessages(library(lidR))
suppressMessages(library(raster))
suppressMessages(library(data.table))
suppressMessages(library(sf))
suppressMessages(library(terra))

# =============================================================================
# Gapiness metrics script
#
# This script calculates two families of gap metrics:
#
#   1. VOXEL-BASED GAP METRICS
#      A 1 m x 1 m raster cell is treated as one vertical column.
#      Inside each column, the Z dimension is divided into 1 m vertical voxels.
#
#      A voxel is considered "full" if it contains at least nOfElements points.
#
#      A column is a gap if there is NO full voxel above the selected threshold.
#
#   2. POINT-COUNT GAP METRICS
#      A 1 m x 1 m raster cell is treated as one vertical column.
#
#      A column is a gap if it contains fewer than nOfElements total points
#      above the selected threshold.
#
# Thresholds:
#
#   DeepGapFraction:
#      threshold = 1 m
#
#   GapFraction:
#      threshold = 0.25 * max_height_global
#
# No-return behavior:
#
#   Cells with no returns inside the polygon are treated as gap = 1.
#   Cells outside the polygon are masked to NA and excluded from percentages.
#
# Optional cover metrics:
#
#   CoverFraction = 1 - GapFraction
#
#   Fully functional cover-fraction code is included below, but commented out.
# =============================================================================


# =============================================================================
# Read command line arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  stop(
    "Usage: Rscript GappinesMetrics.R ",
    "<input_folder> <output_folder> <max_height_file> ",
    "<output_prefix> <polygon_file> <numberOfElements> <LOCAL|EPSG:code>"
  )
}

las_folder       <- args[1]
output_folder    <- args[2]
max_height_file  <- args[3]
base_name        <- args[4]
polygon_file     <- args[5]

# Command line arguments are read as text.
# Convert numberOfElements to integer because it is used as a point-count threshold.
nOfElements <- as.integer(args[6])

if (is.na(nOfElements) || nOfElements < 1) {
  stop("numberOfElements must be a positive integer.")
}

coordinates <- trimws(args[7])

if (is.na(coordinates) || !nzchar(coordinates)) {
  stop("Coordinate system argument is empty. Use LOCAL or a valid CRS such as EPSG:3035.")
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

cat("Output folder: ", output_folder, "\n")
cat("Polygon file: ", polygon_file, "\n")
cat("Number of elements: ", nOfElements, "\n")
cat("Coordinate system: ", coordinates, "\n")
cat(
  "Coordinate mode: ",
  if (local_mode) "LOCAL metric Cartesian coordinates (no reprojection)" else "georeferenced CRS",
  "\n",
  sep = ""
)


# =============================================================================
# Read global max height
# =============================================================================

max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))

if (is.na(max_height_global)) {
  stop("Invalid max height value in the file.")
}

gap_thresh <- 0.25 * max_height_global

cat("Global maximum height: ", max_height_global, "\n")
cat("GapFraction threshold, 25% of max height: ", gap_thresh, "\n")


# =============================================================================
# Prepare output folder
# =============================================================================

raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# Coordinate system
# =============================================================================

# Fallback CRS in georeferenced mode. It is deliberately NA in LOCAL mode.
crs_code <- configured_crs

format_crs_for_log <- function(crs_obj) {
  if (is.null(crs_obj) || is.na(crs_obj)) return("NA")
  epsg <- tryCatch(crs_obj$epsg, error = function(e) NA_integer_)
  input <- tryCatch(crs_obj$input, error = function(e) NA_character_)
  wkt <- tryCatch(crs_obj$wkt, error = function(e) NA_character_)

  parts <- character(0)
  if (!is.na(epsg)) parts <- c(parts, paste0("EPSG:", epsg))
  if (!is.na(input) && nzchar(input)) parts <- c(parts, paste0("input=", input))
  if (!length(parts) && !is.na(wkt) && nzchar(wkt)) {
    parts <- c(parts, paste0("wkt=", substr(wkt, 1, 160)))
  }
  if (!length(parts)) return("NA")
  paste(parts, collapse = "; ")
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


# =============================================================================
# Load LAS / LAZ files
# =============================================================================

cat("Loading LAS/LAZ files...\n")

las_files <- list.files(
  las_folder,
  pattern = "\\.(las|laz)$",
  full.names = TRUE
)

# Ignore zero-byte files.
las_files <- las_files[file.info(las_files)$size > 0]

if (length(las_files) == 0) {
  stop("No valid LAS/LAZ files found.")
}


# =============================================================================
# Read and merge LAS files
# =============================================================================

cat("Reading and merging LAS files...\n")

las_list <- lapply(las_files, function(file) {
  las <- tryCatch(
    readLAS(file),
    error = function(e) {
      warning("Failed to read LAS/LAZ file: ", file)
      return(NULL)
    }
  )

  if (is.null(las) || npoints(las) == 0) {
    return(NULL)
  }

  return(las)
})

# Remove failed or empty reads before merging.
las_list <- Filter(Negate(is.null), las_list)

if (!length(las_list)) {
  stop("No points found in LAS/LAZ files.")
}

las_combined <- do.call(rbind, las_list)

if (is.null(las_combined) || npoints(las_combined) == 0) {
  stop("No points found after merging LAS/LAZ files.")
}


# =============================================================================
# Check / set LAS CRS
# =============================================================================

# Get CRS from merged LAS. In LOCAL mode the raw values are used with an
# undefined CRS. In georeferenced mode, missing LAS metadata receive the CRS
# explicitly supplied for this run.
las_crs_obj <- tryCatch(lidR::st_crs(las_combined), error = function(e) NULL)

if (local_mode) {
  cat("LOCAL mode; LAS CRS metadata ignored: ", format_crs_for_log(las_crs_obj), "\n")
  lidR::st_crs(las_combined) <- sf::NA_crs_
  las_crs_obj <- sf::NA_crs_
} else if (is.null(las_crs_obj) || is.na(las_crs_obj)) {
  lidR::st_crs(las_combined) <- crs_code
  las_crs_obj <- lidR::st_crs(las_combined)
  cat("LAS had no CRS, setting fallback: ", coordinates, "\n")
} else {
  cat("LAS CRS found: ", format_crs_for_log(las_crs_obj), "\n")
}


# =============================================================================
# Read and prepare polygon
# =============================================================================

cat("Loading polygon...\n")

poly <- sf::st_read(polygon_file, quiet = TRUE)

poly_crs_obj <- sf::st_crs(poly)

if (local_mode) {
  cat("LOCAL mode; polygon CRS metadata ignored: ", format_crs_for_log(poly_crs_obj), "\n")
  suppressWarnings(sf::st_crs(poly) <- sf::NA_crs_)
} else if (is.null(poly_crs_obj) || is.na(poly_crs_obj)) {
  # The polygon and LAS are expected to be numerically aligned. Attach the
  # effective LAS CRS as missing metadata without moving any coordinates.
  suppressWarnings(sf::st_crs(poly) <- las_crs_obj)
  poly_crs_obj <- sf::st_crs(poly)
  cat(
    "Polygon had no CRS, assigning effective LAS CRS: ",
    format_crs_for_log(las_crs_obj),
    "\n"
  )
} else {
  cat("Polygon CRS before transform: ", format_crs_for_log(poly_crs_obj), "\n")
}

if (!local_mode) {
  cat("LAS CRS used for polygon transform: ", format_crs_for_log(las_crs_obj), "\n")
  cat("Reprojecting polygon to LAS CRS...\n")
} else {
  cat("LOCAL mode: polygon coordinates will not be transformed.\n")
}

# Repair invalid geometries and drop Z/M dimensions. Reproject only in
# georeferenced mode; LOCAL X/Y values must remain numerically unchanged.
poly <- sf::st_make_valid(poly)
suppressWarnings({
  poly <- sf::st_zm(poly, drop = TRUE, what = "ZM")
})
if (!local_mode) {
  poly <- sf::st_transform(poly, crs = las_crs_obj)
}

# Dissolve all polygon features into one geometry.
poly <- suppressWarnings(sf::st_union(poly))

# -----------------------------------------------------------------------------
# Remove a 1 m border strip from the calculation area.
#
# The original polygon is preserved in `poly`. A second polygon (`poly_calc`) is
# created by buffering 1 m toward the polygon interior. Only cells inside this
# eroded polygon are retained in the final gap rasters and summary percentages.
#
# In georeferenced mode the polygon has been transformed to the LAS CRS. In
# LOCAL mode its raw X/Y values are unchanged. In both cases, the pipeline
# requires horizontal coordinate units to be metres, so dist = -1 means 1 m.
# -----------------------------------------------------------------------------
inner_buffer_m <- 1

if (isTRUE(sf::st_is_longlat(poly))) {
  stop(
    "Cannot apply a 1 m inward polygon buffer because the LAS CRS is geographic ",
    "(longitude/latitude). Please use a projected CRS with metre units."
  )
}

cat("Applying inward polygon buffer: ", inner_buffer_m, " m\n", sep = "")

poly_calc <- suppressWarnings(
  sf::st_buffer(poly, dist = -inner_buffer_m)
)
poly_calc <- sf::st_make_valid(poly_calc)

# A negative buffer can completely remove very small or narrow polygons.
if (length(poly_calc) == 0 || all(sf::st_is_empty(poly_calc))) {
  stop(
    "The polygon becomes empty after applying the 1 m inward buffer. ",
    "The site is too small/narrow for this buffer distance."
  )
}


# =============================================================================
# Non-ground and height filter
# =============================================================================
#
# This preserves the original behavior:
#
#   - Ground points are excluded: Classification != 2
#   - Points at or below 1 m are excluded: Z > 1
#
# Both the voxel-based and point-count metrics below are therefore calculated
# only from non-ground points above 1 m.
# =============================================================================

cat("Filtering non-ground points above 1 m...\n")

las_nonground <- filter_poi(
  las_combined,
  Classification != 2 & Z > 1
)

if (is.null(las_nonground) || npoints(las_nonground) == 0) {
  warning("No non-ground points above 1 m remain after filtering.")
}


# =============================================================================
# VOXEL-BASED GAP METRIC FUNCTIONS
# =============================================================================
#
# grid_metrics(..., res = 1) creates 1 m x 1 m horizontal raster cells.
#
# Each raster cell is treated as one vertical column.
#
# For each column:
#
#   1. Keep only points above the selected height threshold.
#   2. Divide Z values into 1 m vertical bins using floor(Z).
#   3. Count how many points fall inside each vertical voxel.
#   4. Check whether at least one voxel contains nOfElements or more points.
#
# Output:
#
#   0 = not a gap
#   1 = gap
#
# A column is not a gap if it has at least one full 1 x 1 x 1 m voxel
# above the selected threshold.
# =============================================================================

has_full_voxel_above_threshold <- function(z, threshold, min_points) {
  # Remove NA, NaN, Inf and -Inf values.
  z <- z[is.finite(z)]
  
  # If there are no valid points, there is no full voxel.
  if (!length(z)) {
    return(FALSE)
  }
  
  # Keep only points above the relevant height threshold.
  z <- z[z > threshold]
  
  # If no points remain above the threshold, the column cannot be full.
  if (!length(z)) {
    return(FALSE)
  }
  
  # Create 1 m vertical voxel IDs.
  #
  # Examples:
  #   Z in [1, 2) -> voxel_id = 1
  #   Z in [2, 3) -> voxel_id = 2
  #   Z in [3, 4) -> voxel_id = 3
  voxel_id <- as.integer(floor(z))
  
  # Count how many points are inside each vertical voxel.
  #
  # The offset makes voxel IDs start at 1 for tabulate().
  voxel_counts <- tabulate(voxel_id - min(voxel_id) + 1L)
  
  # A column contains a full voxel if at least one voxel has enough points.
  any(voxel_counts >= min_points)
}


voxel_gap_metric <- function(z, threshold, min_points) {
  # Return binary gap value:
  #
  #   0 = not a gap
  #   1 = gap
  #
  # A column is a gap if it contains no full voxel above the threshold.
  if (has_full_voxel_above_threshold(z, threshold, min_points)) {
    return(0)
  } else {
    return(1)
  }
}


# =============================================================================
# POINT-COUNT GAP METRIC FUNCTION
# =============================================================================
#
# This reproduces the older/simple metric definition.
#
# For each 1 m x 1 m raster cell:
#
#   1. Keep only points above the selected height threshold.
#   2. Count total number of points above that threshold.
#   3. Compare the count with nOfElements.
#
# Output:
#
#   0 = not a gap
#   1 = gap
#
# A column is not a gap if it has at least nOfElements total points
# above the selected threshold.
#
# Unlike the voxel method, this does NOT require points to be concentrated
# inside the same 1 m vertical voxel.
# =============================================================================

point_count_gap_metric <- function(z, threshold, min_points) {
  # Remove NA, NaN, Inf and -Inf values.
  z <- z[is.finite(z)]
  
  # Count all points above the selected threshold.
  n_above_threshold <- sum(z > threshold)
  
  # Return binary gap value:
  #
  #   0 = not a gap
  #   1 = gap
  if (n_above_threshold < min_points) {
    return(1)
  } else {
    return(0)
  }
}


# =============================================================================
# Calculate gap rasters over full LAS extent
# =============================================================================
#
# Base rasters are produced from the NON-CLIPPED LAS, preserving the original
# behavior. They are filled and masked by the polygon later.
# =============================================================================

cat("Calculating voxel-based gap rasters...\n")

# -----------------------------------------------------------------------------
# DeepGapFraction_voxels
#
# A 1 m x 1 m column is a deep gap if it has no full 1 x 1 x 1 m voxel
# above 1 m height.
# -----------------------------------------------------------------------------

r_deep_gap_voxels <- grid_metrics(
  las_nonground,
  ~voxel_gap_metric(
    Z,
    threshold = 1,
    min_points = nOfElements
  ),
  res = 1
)


# -----------------------------------------------------------------------------
# GapFraction_voxels
#
# A 1 m x 1 m column is a gap if it has no full 1 x 1 x 1 m voxel
# above 25% of the global maximum height.
# -----------------------------------------------------------------------------

r_gap_voxels <- grid_metrics(
  las_nonground,
  ~voxel_gap_metric(
    Z,
    threshold = gap_thresh,
    min_points = nOfElements
  ),
  res = 1
)


cat("Calculating point-count gap rasters...\n")

# -----------------------------------------------------------------------------
# DeepGapFraction_points
#
# A 1 m x 1 m column is a deep gap if it has fewer than nOfElements total points
# above 1 m height.
# -----------------------------------------------------------------------------

r_deep_gap_points <- grid_metrics(
  las_nonground,
  ~point_count_gap_metric(
    Z,
    threshold = 1,
    min_points = nOfElements
  ),
  res = 1
)


# -----------------------------------------------------------------------------
# GapFraction_points
#
# A 1 m x 1 m column is a gap if it has fewer than nOfElements total points
# above 25% of the global maximum height.
# -----------------------------------------------------------------------------

r_gap_points <- grid_metrics(
  las_nonground,
  ~point_count_gap_metric(
    Z,
    threshold = gap_thresh,
    min_points = nOfElements
  ),
  res = 1
)

# Make the output CRS policy explicit before converting to terra objects.
r_deep_gap_voxels <- set_output_crs(r_deep_gap_voxels, las_crs_obj)
r_gap_voxels      <- set_output_crs(r_gap_voxels,      las_crs_obj)
r_deep_gap_points <- set_output_crs(r_deep_gap_points, las_crs_obj)
r_gap_points      <- set_output_crs(r_gap_points,      las_crs_obj)


# =============================================================================
# Clip rasters with polygon using terra
# =============================================================================

cat("Filling no-return cells and clipping rasters with polygon...\n")

# Convert the inward-buffered calculation polygon to a terra vector.
# terra keeps the CRS from the sf object.
# The original polygon boundary is therefore excluded by a 1 m border strip.
poly_v <- vect(poly_calc)


# -----------------------------------------------------------------------------
# Helper function:
#
#   1. Convert raster::RasterLayer to terra::SpatRaster.
#   2. Convert NA cells to gap = 1.
#   3. Mask with polygon.
#   4. Convert back to raster::RasterLayer.
#
# This preserves the intended behavior:
#
#   NA inside polygon  -> gap = 1
#   outside polygon    -> NA
# -----------------------------------------------------------------------------

fill_no_returns_and_mask <- function(r, polygon_vect) {
  s <- rast(r)
  
  # Treat no-return cells as gaps before polygon masking.
  s_fill <- ifel(
    is.na(s),
    1,
    s
  )
  
  # Keep only cells inside the polygon.
  # Outside the polygon becomes NA and is excluded from percentages.
  s_clip <- mask(
    s_fill,
    polygon_vect
  )
  
  return(raster(s_clip))
}


# Apply the same fill-and-mask behavior to all gap rasters.
r_deep_gap_voxels <- fill_no_returns_and_mask(r_deep_gap_voxels, poly_v)
r_gap_voxels      <- fill_no_returns_and_mask(r_gap_voxels,      poly_v)
r_deep_gap_points <- fill_no_returns_and_mask(r_deep_gap_points, poly_v)
r_gap_points      <- fill_no_returns_and_mask(r_gap_points,      poly_v)


# =============================================================================
# Optional cover rasters
# =============================================================================
#
# CoverFraction is defined here as:
#
#   cover = 1 - gap
#
# These lines are fully functional but commented out.
#
# If you want cover outputs, uncomment this whole section and the related
# writeRaster(), percentage, and PNG sections below.
# =============================================================================

# r_cover_deep_voxels <- 1 - r_deep_gap_voxels
# r_cover_gap_voxels  <- 1 - r_gap_voxels
# r_cover_deep_points <- 1 - r_deep_gap_points
# r_cover_gap_points  <- 1 - r_gap_points


# =============================================================================
# Save output rasters
# =============================================================================

cat("Saving gap rasters...\n")

writeRaster(
  r_deep_gap_voxels,
  file.path(raster_output_path, paste0(base_name, "_DeepGapFraction_voxels.tif")),
  overwrite = TRUE
)

writeRaster(
  r_gap_voxels,
  file.path(raster_output_path, paste0(base_name, "_GapFraction_voxels.tif")),
  overwrite = TRUE
)

writeRaster(
  r_deep_gap_points,
  file.path(raster_output_path, paste0(base_name, "_DeepGapFraction_points.tif")),
  overwrite = TRUE
)

writeRaster(
  r_gap_points,
  file.path(raster_output_path, paste0(base_name, "_GapFraction_points.tif")),
  overwrite = TRUE
)


# -----------------------------------------------------------------------------
# Optional cover raster outputs
#
# Uncomment if cover rasters are needed.
# -----------------------------------------------------------------------------

# cat("Saving optional cover rasters...\n")
#
# writeRaster(
#   r_cover_deep_voxels,
#   file.path(raster_output_path, paste0(base_name, "_CoverFraction_deep_voxels.tif")),
#   overwrite = TRUE
# )
#
# writeRaster(
#   r_cover_gap_voxels,
#   file.path(raster_output_path, paste0(base_name, "_CoverFraction_gap_voxels.tif")),
#   overwrite = TRUE
# )
#
# writeRaster(
#   r_cover_deep_points,
#   file.path(raster_output_path, paste0(base_name, "_CoverFraction_deep_points.tif")),
#   overwrite = TRUE
# )
#
# writeRaster(
#   r_cover_gap_points,
#   file.path(raster_output_path, paste0(base_name, "_CoverFraction_gap_points.tif")),
#   overwrite = TRUE
# )


# =============================================================================
# Summary metrics
# =============================================================================
#
# Final percentages are calculated only inside the 1 m inward-buffered polygon.
#
# Cells outside the inward-buffered polygon are NA and excluded from the denominator.
#
# For each binary raster:
#
#   percent = 100 * number_of_cells_equal_to_1 / number_of_valid_polygon_cells
# =============================================================================

binary_percent <- function(r) {
  vals <- values(r)
  
  valid_cells <- sum(!is.na(vals))
  
  if (valid_cells == 0) {
    return(NA_real_)
  }
  
  round(
    100 * sum(vals == 1, na.rm = TRUE) / valid_cells,
    2
  )
}


metrics <- data.frame(
  DeepGapFraction_voxels_percent = binary_percent(r_deep_gap_voxels),
  GapFraction_voxels_percent     = binary_percent(r_gap_voxels),
  DeepGapFraction_points_percent = binary_percent(r_deep_gap_points),
  GapFraction_points_percent     = binary_percent(r_gap_points)
  
  # Optional cover percentages.
  # Uncomment these lines if the optional cover rasters above are enabled.
  #
  # CoverFraction_deep_voxels_percent = binary_percent(r_cover_deep_voxels),
  # CoverFraction_gap_voxels_percent  = binary_percent(r_cover_gap_voxels),
  # CoverFraction_deep_points_percent = binary_percent(r_cover_deep_points),
  # CoverFraction_gap_points_percent  = binary_percent(r_cover_gap_points)
)

write.table(
  metrics,
  file = file.path(output_folder, paste0(base_name, "_gappines_metrics.txt")),
  row.names = FALSE,
  col.names = TRUE,
  sep = "        ",
  quote = FALSE
)


# =============================================================================
# PNG previews
# =============================================================================

cat("Saving PNG previews...\n")

plot_binary_raster_png <- function(r, file_name, plot_title) {
  png(
    file.path(output_folder, file_name),
    width = 800,
    height = 600
  )
  
  plot(
    r,
    main = plot_title,
    xlab = "X",
    ylab = "Y",
    zlim = c(0, 1)
  )
  
  dev.off()
}


# Gap PNGs: voxel method
plot_binary_raster_png(
  r_deep_gap_voxels,
  paste0(base_name, "_DeepGapFraction_voxels.png"),
  "Deep Gap Fraction - Voxels (Binary)"
)

plot_binary_raster_png(
  r_gap_voxels,
  paste0(base_name, "_GapFraction_voxels.png"),
  "Gap Fraction - Voxels (Binary)"
)


# Gap PNGs: point-count method
plot_binary_raster_png(
  r_deep_gap_points,
  paste0(base_name, "_DeepGapFraction_points.png"),
  "Deep Gap Fraction - Points (Binary)"
)

plot_binary_raster_png(
  r_gap_points,
  paste0(base_name, "_GapFraction_points.png"),
  "Gap Fraction - Points (Binary)"
)


# -----------------------------------------------------------------------------
# Optional cover PNGs
#
# Uncomment if cover rasters are enabled.
# -----------------------------------------------------------------------------

# plot_binary_raster_png(
#   r_cover_deep_voxels,
#   paste0(base_name, "_CoverFraction_deep_voxels.png"),
#   "Cover Fraction from Deep Gap - Voxels (Binary)"
# )
#
# plot_binary_raster_png(
#   r_cover_gap_voxels,
#   paste0(base_name, "_CoverFraction_gap_voxels.png"),
#   "Cover Fraction from Gap - Voxels (Binary)"
# )
#
# plot_binary_raster_png(
#   r_cover_deep_points,
#   paste0(base_name, "_CoverFraction_deep_points.png"),
#   "Cover Fraction from Deep Gap - Points (Binary)"
# )
#
# plot_binary_raster_png(
#   r_cover_gap_points,
#   paste0(base_name, "_CoverFraction_gap_points.png"),
#   "Cover Fraction from Gap - Points (Binary)"
# )


cat("Gap metrics processing complete. Results saved to:", output_folder, "\n")
