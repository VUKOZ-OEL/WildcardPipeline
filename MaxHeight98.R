library(lidR)
library(raster)
library(sf)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: Rscript MaxHeight98.R <input_folder> <output_folder> <output_prefix> <LOCAL|EPSG:code>")
}
input_folder  <- args[1]
output_folder <- args[2]
output_prefix <- args[3]
crs_string    <- trimws(args[4])

if (is.na(crs_string) || !nzchar(crs_string)) {
  stop("crs_string must be LOCAL or a valid CRS such as EPSG:3035.")
}

local_mode <- identical(toupper(crs_string), "LOCAL")
if (local_mode) {
  crs_string <- "LOCAL"
  configured_crs <- sf::NA_crs_
} else {
  configured_crs <- suppressWarnings(sf::st_crs(crs_string))
  if (is.na(configured_crs)) {
    stop("Invalid coordinate system. Use LOCAL or a valid CRS such as EPSG:3035; got: ", crs_string)
  }
}

crs_missing <- function(x) is.null(x) || is.na(x)
crs_text_for_raster <- function(x) {
  if (crs_missing(x)) return(crs_string)
  input <- tryCatch(x$input, error = function(e) NA_character_)
  wkt <- tryCatch(x$wkt, error = function(e) NA_character_)
  if (!is.na(input) && nzchar(input)) return(input)
  if (!is.na(wkt) && nzchar(wkt)) return(wkt)
  crs_string
}
set_output_crs <- function(r, source_crs) {
  if (local_mode) {
    raster::crs(r) <- NA_character_
  } else {
    raster::crs(r) <- crs_text_for_raster(source_crs)
  }
  r
}

if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

cat("Step: Calculating Top Canopy Height (98th percentile of Z > 1m non-ground points)...\n")

las_files <- list.files(input_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]
if (length(las_files) == 0) stop("No valid LAS/LAZ files found.")

cat("Found", length(las_files), "LAS/LAZ file(s) in", input_folder, "\n")

cat("Reading as LAScatalog (chunked processing for large files)...\n")
# IMPORTANT: read the whole folder, not a single file
ctg <- readLAScatalog(input_folder)

ctg_crs_obj <- tryCatch(lidR::st_crs(ctg), error = function(e) NULL)
if (local_mode) {
  lidR::st_crs(ctg) <- sf::NA_crs_
  ctg_crs_obj <- sf::NA_crs_
  cat("Coordinate mode: LOCAL; raster coordinates remain unreferenced.\n")
} else if (crs_missing(ctg_crs_obj)) {
  lidR::st_crs(ctg) <- configured_crs
  ctg_crs_obj <- lidR::st_crs(ctg)
  cat("LAS catalog had no CRS; using configured fallback: ", crs_string, "\n", sep = "")
}

opt_chunk_size(ctg)   <- 1000
opt_chunk_buffer(ctg) <- 10

# IMPORTANT: avoid read-time filtering that can empty whole chunks (state ∅)
opt_filter(ctg) <- ""

p98_canopy <- function(z, cls) {
  zz <- z[cls != 2L & z > 1]              # non-ground AND above 1 m
  if (length(zz) == 0L) return(NA_real_)  # no canopy points in this cell
  as.numeric(stats::quantile(zz, 0.98, na.rm = TRUE, names = FALSE))
}

cat("Computing 98th percentile canopy height raster...\n")
r1 <- grid_metrics(ctg, ~p98_canopy(Z, Classification), res = 1)

r1 <- set_output_crs(r1, ctg_crs_obj)

# Histogram only if at least one finite cell exists
vals <- values(r1)
finite_n <- sum(is.finite(vals))
if (finite_n > 0) {
  png(
    filename = file.path(output_folder, paste0(output_prefix, "_height98_hist.png")),
    width    = 1000,
    height   = 800
  )
  raster::hist(r1, main = "Max height histogram", xlab = "Height (m)", ylab = "Frequency")
  dev.off()
} else {
  cat("Histogram skipped: raster has no finite values (no canopy points >1m).\n")
}

cat("Extracting 98th percentile value of finite raster cells...\n")
finite_vals <- vals[is.finite(vals)]
p98_raster_top_canopy_height <- if (length(finite_vals) > 0L) {
  as.numeric(stats::quantile(finite_vals, probs = 0.98, na.rm = TRUE, names = FALSE))
} else {
  NA_real_
}

cat("Global 98th Percentile Top Canopy Height across raster cells:", p98_raster_top_canopy_height, "\n")

cat("Writing final output...\n")
output_file <- file.path(output_folder, paste0(output_prefix, "_maxHeight.txt"))
write(p98_raster_top_canopy_height, file = output_file, ncolumns = 1)

cat("Processing done. Moving forward with analysis.\n")
