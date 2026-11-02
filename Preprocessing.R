suppressPackageStartupMessages({
  library(lidR)
  library(raster)
  library(terra)
  library(ggplot2)
  library(sf)      # for polygon handling
})

terra::terraOptions(overwrite = TRUE)

# ---------- Pin native/OpenMP threads ----------
Sys.setenv(
  OMP_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS      = "1",
  BLIS_NUM_THREADS     = "1",
  OMP_WAIT_POLICY      = "PASSIVE"
)
if (utils::packageVersion("lidR") >= "2.1.0") {
  try(lidR::set_lidr_threads(1L), silent = TRUE)
}

# ---------- Force sequential chunk processing if {future} exists ----------
.old_plan <- NULL
if (requireNamespace("future", quietly = TRUE)) {
  .old_plan <- future::plan()
  future::plan("sequential")
}
.onexit_future <- function() {
  if (!is.null(.old_plan)) try(future::plan(.old_plan), silent = TRUE)
}
on.exit(.onexit_future(), add = TRUE)

# ---------- Memory logger (RSS + R heap) ----------
.memlog <- data.frame(step = character(), rss_MiB = numeric(), r_MiB = numeric(), stringsAsFactors = FALSE)
.as_mb  <- function(x) if (is.na(x)) "NA" else sprintf("%.2f", x)
.mem_used <- function() {
  rss <- NA_real_
  if (file.exists("/proc/self/status")) {
    ln <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
    if (length(ln)) {
      kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", ln[1])))
      if (!is.na(kb)) rss <- kb / 1024
    }
  }
  r <- NA_real_
  if (requireNamespace("lobstr", quietly = TRUE)) r <- as.numeric(lobstr::mem_used()) / (1024^2)
  list(rss = rss, r = r)
}
.mem_snap <- function(step) {
  step <- gsub("[\t\r\n]", " ", as.character(step))
  gc()
  u <- tryCatch(.mem_used(), error = function(e) list(rss = NA_real_, r = NA_real_))
  .memlog <<- rbind(.memlog, data.frame(step = step, rss_MiB = u$rss, r_MiB = u$r))
  cat(sprintf("[MEM]\tstep=%s\trss_MiB=%s\tr_MiB=%s\n", step, .as_mb(u$rss), .as_mb(u$r)))
  flush(stdout())
}
dump_memlog <- function() {
  if (exists(".memlog", inherits = FALSE) && nrow(.memlog)) {
    cat("\n===== MEM SUMMARY =====\n")
    write.table(.memlog, file = stdout(), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("=======================\n")
    flush(stdout())
  }
}
on.exit(dump_memlog(), add = TRUE)

# also dump last warnings at end
on.exit({
  ws <- warnings()
  if (!is.null(ws)) {
    cat("\n===== LAST WARNINGS (max 50) =====\n")
    print(utils::tail(ws, 50))
    cat("===================================\n")
  }
}, add = TRUE)

# Set unique chunk-raster output template per product
set_raster_chunk_output <- function(ctg, out_dir, prefix) {
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  lidR::opt_output_files(ctg) <- file.path(out_dir, paste0(prefix, "_{XLEFT}_{YBOTTOM}_{ID}"))
  ctg
}


# ---------- Args ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript preprocess_single_las.R <las_folder_path> <processed_folder_path> <results_folder_path> <polygon_file> <LOCAL|EPSG:code> [max_normalized_height]")
}
las_folder_path       <- args[1]
processed_folder_path <- args[2]
results_folder_path   <- args[3]
coordinates           <- trimws(args[5])
if (is.na(coordinates) || !nzchar(coordinates)) {
  stop("Missing coordinates argument.")
}

# Supported modes:
#   LOCAL       - aligned local Cartesian coordinates in metres; no CRS is used
#   EPSG:<code> - georeferenced data; this CRS is used as a fallback when input
#                 metadata are missing
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

max_normalized_height_arg <- if (length(args) >= 6L) trimws(args[6]) else ""
max_normalized_height <- NA_real_
if (nzchar(max_normalized_height_arg)) {
  max_normalized_height <- suppressWarnings(as.numeric(max_normalized_height_arg))
  if (length(max_normalized_height) != 1L ||
      !is.finite(max_normalized_height) ||
      max_normalized_height < 0) {
    stop(
      "max_normalized_height must be empty or a finite non-negative number; got: ",
      max_normalized_height_arg
    )
  }
}

cat("Coordinates: ", coordinates, "\n")
cat(
  "Coordinate mode: ",
  if (local_mode) "LOCAL metric Cartesian coordinates (no reprojection)" else "georeferenced CRS",
  "\n",
  sep = ""
)
if (is.na(max_normalized_height)) {
  cat("Max normalized height: disabled\n")
} else {
  cat("Max normalized height: ", max_normalized_height, " m above ground\n", sep = "")
}

# Optional polygon argument
polygon_file <- if (length(args) >= 4) args[4] else NA_character_
use_polygon  <- !is.na(polygon_file) && nzchar(polygon_file)
if (use_polygon && !file.exists(polygon_file)) {
  stop("Polygon file does not exist: ", polygon_file)
}

# Fallback CRS in georeferenced mode. It is deliberately NA in LOCAL mode.
crs_code <- configured_crs

dir.create(processed_folder_path, showWarnings = FALSE, recursive = TRUE)
dir.create(results_folder_path,   showWarnings = FALSE, recursive = TRUE)
.mem_snap("Setup complete")

# ---------- Pick first LAS/LAZ ----------
las_files <- list.files(las_folder_path, pattern = "\\.(las|laz)$", full.names = TRUE)
if (!length(las_files)) stop("No LAS/LAZ in: ", las_folder_path)
las_file  <- las_files[1]
base_name <- basename(las_file)
stem      <- tools::file_path_sans_ext(base_name)
ext       <- tolower(tools::file_ext(base_name))
write_laz <- (ext == "laz")
.mem_snap(paste("Processing file:", las_file))

# ---------- Work dirs + raster tmpdir ----------
work_dir  <- file.path(results_folder_path, ".work", stem)
RESUME    <- as.integer(Sys.getenv("RESUME", "1"))
if (RESUME == 0L && dir.exists(work_dir)) unlink(work_dir, recursive = TRUE, force = TRUE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
class_dir <- file.path(work_dir, "class")
norm_dir  <- file.path(work_dir, "norm")
dir.create(class_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(norm_dir,  showWarnings = FALSE, recursive = TRUE)

run_tag <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", as.integer(Sys.getpid()))
tmpdir  <- file.path(work_dir, paste0("raster_tmp_", run_tag))
dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)
raster::rasterOptions(tmpdir = tmpdir, overwrite = TRUE, progress = "text")
.mem_snap(paste("Raster tmpdir:", tmpdir))

cleanup_temp_outputs <- function() {
  if (exists("work_dir") && dir.exists(work_dir)) {
    try(unlink(work_dir, recursive = TRUE, force = TRUE), silent = TRUE)
  }
  if (exists("work_dir")) {
    work_parent <- dirname(work_dir)
    if (dir.exists(work_parent) && length(list.files(work_parent, all.files = FALSE)) == 0L) {
      try(unlink(work_parent, recursive = TRUE, force = TRUE), silent = TRUE)
    }
  }
  if (exists("processed_folder_path") && exists("stem")) {
    leftovers <- Sys.glob(file.path(processed_folder_path, paste0(stem, "_merged*")))
    if (length(leftovers)) {
      try(unlink(leftovers, recursive = TRUE, force = TRUE), silent = TRUE)
    }
  }
  if (exists("results_folder_path")) {
    norm_dtm_tmp <- file.path(results_folder_path, "DTM_unmasked_for_normalization.tif")
    if (file.exists(norm_dtm_tmp)) {
      try(unlink(norm_dtm_tmp, force = TRUE), silent = TRUE)
    }
  }
}
on.exit(cleanup_temp_outputs(), add = TRUE)

build_plot_mask <- function(template_raster, polygon_path, local_mode, fallback_crs) {
  if (is.na(polygon_path) || !nzchar(polygon_path) || !file.exists(polygon_path)) return(NULL)
  poly <- sf::st_read(polygon_path, quiet = TRUE)
  tr_crs <- tryCatch(raster::crs(template_raster), error = function(e) NA)

  if (local_mode) {
    # LOCAL means that raw X/Y values are already aligned. Ignore any attached
    # CRS label and never transform the coordinates.
    suppressWarnings(sf::st_crs(poly) <- sf::NA_crs_)
  } else if (is.na(sf::st_crs(poly))) {
    # A missing .prj is a metadata omission. The polygon is already aligned
    # with the point cloud, so prefer the template/LAS CRS when available.
    if (!is.na(tr_crs) && nzchar(as.character(tr_crs))) {
      suppressWarnings(sf::st_crs(poly) <- sf::st_crs(tr_crs))
    } else {
      suppressWarnings(sf::st_crs(poly) <- fallback_crs)
    }
  }

  poly <- sf::st_make_valid(poly)
  suppressWarnings({
    poly <- sf::st_zm(poly, drop = TRUE, what = "ZM")
  })

  if (!local_mode) {
    if (!is.na(tr_crs) && nzchar(as.character(tr_crs))) {
      poly <- sf::st_transform(poly, crs = tr_crs)
    } else {
      poly <- sf::st_transform(poly, crs = fallback_crs)
    }
  }

  poly <- suppressWarnings(sf::st_union(poly))
  poly_sp <- methods::as(poly, "Spatial")
  raster::rasterize(poly_sp, template_raster, field = 1, background = NA)
}


format_crs_for_log <- function(crs_obj) {
  if (is.null(crs_obj) || is.na(crs_obj)) return("NA")
  epsg <- tryCatch(crs_obj$epsg, error = function(e) NA_integer_)
  input <- tryCatch(crs_obj$input, error = function(e) NA_character_)
  wkt <- tryCatch(crs_obj$wkt, error = function(e) NA_character_)

  parts <- character(0)
  if (!is.na(epsg)) parts <- c(parts, paste0("EPSG:", epsg))
  if (!is.na(input) && nzchar(input)) parts <- c(parts, paste0("input=", input))
  if (!length(parts) && !is.na(wkt) && nzchar(wkt)) parts <- c(parts, paste0("wkt=", substr(wkt, 1, 160)))
  if (!length(parts)) return("NA")
  paste(parts, collapse = "; ")
}

log_bbox <- function(label, xmin, xmax, ymin, ymax) {
  .mem_snap(sprintf(
    "%s bbox: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f",
    label, xmin, xmax, ymin, ymax
  ))
}


# ---------- Clip LAS to polygon at the beginning (if requested) ----------
if (use_polygon) {
  .mem_snap("Polygon clip: reading LAS and polygon")
  
  las_full <- readLAS(las_file)
  if (is.null(las_full) || lidR::is.empty(las_full)) {
    stop("Input LAS is empty before clipping: ", las_file)
  }
  
  # Resolve the LAS and polygon coordinate handling for this run.
  las_crs_obj <- tryCatch(lidR::st_crs(las_full), error = function(e) NULL)
  poly <- sf::st_read(polygon_file, quiet = TRUE)
  poly_crs_obj <- sf::st_crs(poly)

  if (local_mode) {
    .mem_snap(paste(
      "Polygon clip: LOCAL mode; LAS CRS metadata ignored:",
      format_crs_for_log(las_crs_obj)
    ))
    .mem_snap(paste(
      "Polygon clip: LOCAL mode; polygon CRS metadata ignored:",
      format_crs_for_log(poly_crs_obj)
    ))
    lidR::st_crs(las_full) <- sf::NA_crs_
    suppressWarnings(sf::st_crs(poly) <- sf::NA_crs_)
    las_crs_obj <- sf::NA_crs_
  } else {
    if (is.null(las_crs_obj) || is.na(las_crs_obj)) {
      lidR::st_crs(las_full) <- crs_code
      las_crs_obj <- lidR::st_crs(las_full)
      .mem_snap(paste("Polygon clip: LAS had no CRS, setting fallback:", coordinates))
    } else {
      .mem_snap(paste("Polygon clip: LAS CRS found:", format_crs_for_log(las_crs_obj)))
    }

    if (is.null(poly_crs_obj) || is.na(poly_crs_obj)) {
      # The polygon and LAS are expected to be numerically aligned. When the
      # .prj is absent, attach the effective LAS CRS without moving coordinates.
      suppressWarnings(sf::st_crs(poly) <- las_crs_obj)
      poly_crs_obj <- sf::st_crs(poly)
      .mem_snap(paste(
        "Polygon clip: polygon had no CRS, assigning effective LAS CRS:",
        format_crs_for_log(las_crs_obj)
      ))
    } else {
      .mem_snap(paste("Polygon clip: polygon CRS before transform:", format_crs_for_log(poly_crs_obj)))
    }
    .mem_snap(paste("Polygon clip: LAS CRS used for clipping:", format_crs_for_log(las_crs_obj)))
  }

  # Clean geometry. Reproject only in georeferenced mode; LOCAL coordinates
  # must remain numerically unchanged.
  poly <- sf::st_make_valid(poly)
  suppressWarnings({
    poly <- sf::st_zm(poly, drop = TRUE, what = "ZM")  # drop Z/M if present
  })
  if (!local_mode) {
    poly <- sf::st_transform(poly, crs = las_crs_obj)
  }
  poly <- suppressWarnings(sf::st_union(poly))
  
  las_xr <- range(las_full@data$X, na.rm = TRUE)
  las_yr <- range(las_full@data$Y, na.rm = TRUE)
  poly_bb <- sf::st_bbox(poly)
  log_bbox("Polygon clip: LAS", las_xr[1], las_xr[2], las_yr[1], las_yr[2])
  log_bbox("Polygon clip: polygon in working coordinates", poly_bb["xmin"], poly_bb["xmax"], poly_bb["ymin"], poly_bb["ymax"])
  
  .mem_snap("Polygon clip: running clip_roi on LAS")
  las_clip <- clip_roi(las_full, poly)
  
  # clip_roi(LAS, MULTIPOLYGON) may return list of LAS -> merge
  if (is.list(las_clip)) {
    las_clip <- do.call(rbind, las_clip)
  }
  
  if (is.null(las_clip) || lidR::is.empty(las_clip)) {
    stop("Polygon clipping produced no points. Check overlap between LAS and polygon.")
  }
  
  clipped_file <- file.path(work_dir, paste0(stem, "_clipped.", ext))
  writeLAS(las_clip, clipped_file)
  .mem_snap(paste("Polygon clip: clipped LAS written:", clipped_file))
  
  rm(las_full, las_clip); gc()
  
  # From here on, the entire pipeline uses the clipped LAS only
  las_file <- clipped_file
  .mem_snap(paste("Polygon clip: using clipped LAS as input:", las_file))
}

# ---------- Tunables (env overrides) ----------
CLASS_CHUNK_M  <- as.numeric(Sys.getenv("CLASS_CHUNK_M",  "300"))
CLASS_BUFFER_M <- as.numeric(Sys.getenv("CLASS_BUFFER_M", "60"))    #30
NORM_CHUNK_M   <- as.numeric(Sys.getenv("NORM_CHUNK_M",   "40"))    #80
NORM_BUFFER_M  <- as.numeric(Sys.getenv("NORM_BUFFER_M",  "10"))    #2
CHM_CHUNK_M    <- as.numeric(Sys.getenv("CHM_CHUNK_M",    "300"))
CHM_BUFFER_M   <- as.numeric(Sys.getenv("CHM_BUFFER_M",   "30"))
NORM_LAZ       <- as.integer(Sys.getenv("NORM_LAZ", "0"))       # 0=LAS (faster), 1=LAZ
DTM_FILL_NA    <- as.integer(Sys.getenv("DTM_FILL_NA", "1"))    # 1=fill isolated NAs
COVERAGE_OK    <- as.numeric(Sys.getenv("COVERAGE_OK", "0.98")) # min coverage to reuse/accept
FORCE_REMERGE  <- as.integer(Sys.getenv("FORCE_REMERGE", "0"))  # manual override for final merge

.mem_snap(sprintf(
  "Chunk params: CLASS=%dm/%dm  NORM=%dm/%dm  CHM=%dm/%dm  NORM_LAZ=%d  FORCE_REMERGE=%d",
  CLASS_CHUNK_M, CLASS_BUFFER_M, NORM_CHUNK_M, NORM_BUFFER_M, CHM_CHUNK_M, CHM_BUFFER_M, NORM_LAZ, FORCE_REMERGE
))

# ---------- Resume helpers ----------
ensure_clean_or_resume <- function(dir, expected_row) {
  meta <- file.path(dir, "_meta.tsv")
  have <- length(list.files(dir, pattern="\\.(las|laz)$", full.names=TRUE)) > 0L
  ok <- FALSE
  if (file.exists(meta)) {
    old <- tryCatch(read.delim(meta, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)
    ok <- !is.null(old) && identical(colnames(old), colnames(expected_row)) &&
      all(as.character(old[1,]) == as.character(expected_row[1,]))
  }
  if (have && !ok) {
    unlink(dir, recursive=TRUE, force=TRUE)
    dir.create(dir, recursive=TRUE, showWarnings=FALSE)
  }
  write.table(expected_row, meta, sep="\t", row.names=FALSE, quote=FALSE)
}
las_header_extent <- function(path) {
  h <- lidR::readLASheader(path)
  c(xmin=h@PHB[["Min X"]], xmax=h@PHB[["Max X"]], ymin=h@PHB[["Min Y"]], ymax=h@PHB[["Max Y"]])
}
union_extent_from_headers <- function(files) {
  if (!length(files)) return(rep(NA_real_, 4))
  e <- las_header_extent(files[[1]])
  for (f in files[-1]) {
    ef <- las_header_extent(f)
    e["xmin"] <- min(e["xmin"], ef["xmin"])
    e["xmax"] <- max(e["xmax"], ef["xmax"])
    e["ymin"] <- min(e["ymin"], ef["ymin"])
    e["ymax"] <- max(e["ymax"], ef["ymax"])
  }
  unname(e)
}
extent_area <- function(ext) (ext[2]-ext[1]) * (ext[4]-ext[3])

ensure_clean_or_resume(class_dir, data.frame(kind="class", chunk=CLASS_CHUNK_M, buffer=CLASS_BUFFER_M))
ensure_clean_or_resume(norm_dir,  data.frame(kind="norm",  chunk=NORM_CHUNK_M,  buffer=NORM_BUFFER_M))

# ---------- Helper: optional DTM NA fill ----------
fill_dtm_if_needed <- function(dtm_raster, label) {
  if (DTM_FILL_NA == 1L) {
    .mem_snap(paste("Step 1:", label, "progressive fill start"))
    wins <- c(3, 5, 9, 15, 31, 61, 121)
    for (w in wins) {
      na_before <- raster::cellStats(is.na(dtm_raster), stat = "sum")
      if (is.na(na_before) || na_before == 0) break
      dtm_filled <- raster::focal(
        dtm_raster, w = matrix(1, w, w), pad = TRUE,
        fun = function(x, ...) {
          center <- x[ceiling(length(x) / 2)]
          if (!is.na(center)) return(center)
          vals <- x[!is.na(x)]
          if (!length(vals)) return(NA_real_)
          mean(vals)
        }
      )
      dtm_raster <- raster::cover(dtm_raster, dtm_filled)
      na_after <- raster::cellStats(is.na(dtm_raster), stat = "sum")
      .mem_snap(sprintf("%s fill window %dx%d: %d -> %d NA cells", label, w, w, as.integer(na_before), as.integer(na_after)))
      if (is.na(na_after) || na_after == 0 || identical(na_after, na_before)) next
    }
    na_cells <- raster::cellStats(is.na(dtm_raster), stat = "sum")
    .mem_snap(sprintf("%s NA cells after progressive fill: %d", label, as.integer(na_cells)))
    .mem_snap(paste("Step 1:", label, "progressive fill done"))
  } else {
    na_cells <- raster::cellStats(is.na(dtm_raster), stat = "sum")
    .mem_snap(sprintf("%s NA cells (no fill): %d", label, as.integer(na_cells)))
  }
  dtm_raster
}


extend_dtm_to_las_extent <- function(dtm_raster, las_path, pad) {
  src_ext <- las_header_extent(las_path)
  dtm_res <- suppressWarnings(raster::res(dtm_raster))
  pad_use <- suppressWarnings(max(c(pad, dtm_res * 2), na.rm = TRUE))
  if (!is.finite(pad_use)) pad_use <- 0

  target_ext <- raster::extent(
    src_ext[1] - pad_use,
    src_ext[2] + pad_use,
    src_ext[3] - pad_use,
    src_ext[4] + pad_use
  )

  current_ext <- raster::extent(dtm_raster)
  current_vec <- c(current_ext[1], current_ext[2], current_ext[3], current_ext[4])
  target_vec  <- c(target_ext[1],  target_ext[2],  target_ext[3],  target_ext[4])

  .mem_snap(sprintf(
    "Step 1: DTM source LAS extent for normalization: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f pad=%.3f",
    src_ext[1], src_ext[2], src_ext[3], src_ext[4], pad_use
  ))
  .mem_snap(sprintf(
    "Step 1: DTM extent before normalization extension: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f",
    current_vec[1], current_vec[2], current_vec[3], current_vec[4]
  ))

  needs_extension <- target_vec[1] < current_vec[1] ||
    target_vec[2] > current_vec[2] ||
    target_vec[3] < current_vec[3] ||
    target_vec[4] > current_vec[4]

  if (isTRUE(needs_extension)) {
    dtm_raster <- raster::extend(dtm_raster, target_ext, value = NA)
    extended_ext <- raster::extent(dtm_raster)
    extended_vec <- c(extended_ext[1], extended_ext[2], extended_ext[3], extended_ext[4])
    .mem_snap(sprintf(
      "Step 1: DTM extended for normalization support: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f",
      extended_vec[1], extended_vec[2], extended_vec[3], extended_vec[4]
    ))
  } else {
    .mem_snap("Step 1: DTM already covers clipped LAS extent for normalization")
  }

  dtm_raster
}

# =========================================================
# STEP 1: Classify ground (CSF) & build DTM — CHUNKED
# =========================================================
ctg <- readLAScatalog(las_file)
opt_chunk_size(ctg)        <- CLASS_CHUNK_M
opt_chunk_buffer(ctg)      <- CLASS_BUFFER_M
opt_progress(ctg)          <- TRUE
opt_laz_compression(ctg)   <- write_laz
opt_output_files(ctg)      <- file.path(class_dir, "{XLEFT}_{YBOTTOM}_{ID}")
opt_select(ctg)            <- "xyz"
.mem_snap("Step 1: catalog ready for CSF")

class_files <- list.files(class_dir, pattern="\\.(las|laz)$", full.names=TRUE)
if (length(class_files) > 0L) {
  .mem_snap("Step 1: reuse existing classified chunks (compat OK)")
  ctg_class <- readLAScatalog(class_dir)
} else {
  .mem_snap("Step 1: classify_ground(CSF) start")
  ctg_class <- classify_ground(ctg, algorithm = csf())
  gc()
  .mem_snap("Step 1: classify_ground(CSF) done")
}

# DTM from classified catalog (ground only), minimal attrs
opt_filter(ctg_class) <- "-keep_class 2"
opt_select(ctg_class) <- "xyz"
.mem_snap("Step 1: DTM prep (filter=ground, select=xyz)")

ctg_class <- set_raster_chunk_output(ctg_class, file.path(work_dir, "dtm1m_tiles"), "DTM1m")
.mem_snap("Step 1: grid_terrain(knnidw) start")
dtm <- grid_terrain(ctg_class, res = 1, algorithm = knnidw(k = 10, p = 2))
gc()
.mem_snap("Step 1: grid_terrain(knnidw) done")

if (local_mode) {
  raster::crs(dtm) <- NA_character_
  .mem_snap("Step 1: LOCAL mode; DTM CRS left undefined")
}

# Extend the unmasked normalization DTM to the full clipped/source LAS extent.
# grid_terrain() builds terrain support from ground points, so its raster extent can be
# slightly smaller than the clipped LAS point extent. Extending before filling prevents
# normalize_height() from seeing edge points outside the DTM raster.
dtm <- extend_dtm_to_las_extent(dtm, las_file, pad = NORM_BUFFER_M)

# Optional: small-hole and edge fill for the 1 m DTM used later
# This also fills NA cells introduced by the normalization-support extension above.
dtm <- fill_dtm_if_needed(dtm, "DTM")

# Build plot mask for output rasters, but DO NOT mask the DTM used for normalization.
# normalize_height() needs terrain support slightly beyond the exact plot footprint
# because chunks use buffers and edge points would otherwise fall on NA terrain.
plot_mask <- NULL
if (use_polygon) {
  .mem_snap("Step 1: building plot mask for terrain outputs")
  plot_mask <- build_plot_mask(dtm, polygon_file, local_mode, crs_code)
}

# Save the unmasked DTM and reopen file-backed. This unmasked DTM is used for normalization.
dtm_norm_tif <- file.path(results_folder_path, "DTM_unmasked_for_normalization.tif")
if (file.exists(dtm_norm_tif)) file.remove(dtm_norm_tif)
writeRaster(dtm, dtm_norm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
dtm <- raster(dtm_norm_tif)
.mem_snap(paste("Step 1: unmasked DTM for normalization written", dtm_norm_tif))

# Prepare masked DTM output for users.
dtm_out <- dtm
if (!is.null(plot_mask)) {
  dtm_out <- raster::mask(dtm_out, plot_mask)
  .mem_snap("Step 1: DTM output masked to plot footprint")
}

dtm_tif <- file.path(results_folder_path, "DTM.tif")
if (file.exists(dtm_tif)) file.remove(dtm_tif)
writeRaster(dtm_out, dtm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 1: DTM written", dtm_tif))

png(file.path(results_folder_path, "DTM.png"), width = 1000, height = 800)
plot(dtm_out, main = "Digital Terrain Model (DTM)")
dev.off()
.mem_snap("Step 1: DTM preview PNG")

# ---------- Additional 0.5 m DTM (derived from masked 1 m DTM output) ----------
.mem_snap("Step 1: DTM_05 from DTM output resample start")
dtm05_template <- raster::raster(
  raster::extent(dtm_out),
  res = 0.5,
  crs = raster::crs(dtm_out)
)
dtm05 <- raster::resample(dtm_out, dtm05_template, method = "bilinear")

if (!is.null(plot_mask)) {
  plot_mask_05 <- raster::resample(plot_mask, dtm05, method = "ngb")
  dtm05 <- raster::mask(dtm05, plot_mask_05)
} else {
  dtm05 <- raster::mask(dtm05, raster::resample(dtm_out, dtm05, method = "ngb"))
}

gc()
.mem_snap("Step 1: DTM_05 from DTM output resample done")

# Save 0.5 m DTM (no further use downstream)
dtm05_tif <- file.path(results_folder_path, "DTM_05.tif")
if (file.exists(dtm05_tif)) file.remove(dtm05_tif)
writeRaster(dtm05, dtm05_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 1: DTM_05 written", dtm05_tif))

# DTM_05 preview as PNG
dtm05_png <- file.path(results_folder_path, "DTM_05.png")
png(dtm05_png, width = 1000, height = 800)
plot(dtm05, main = "Digital Terrain Model 0.5 m (DTM_05)")
dev.off()
.mem_snap(paste("Step 1: DTM_05 preview PNG", dtm05_png))

# =========================================================
# STEP 2: Normalize — CHUNKED with sentinel+coverage → MERGE
# =========================================================
norm_done  <- file.path(norm_dir, "_DONE")
norm_files <- list.files(norm_dir, pattern="\\.(las|laz)$", full.names=TRUE)

need_normalize <- TRUE
if (length(norm_files) > 0L && file.exists(norm_done)) {
  src_ext <- las_header_extent(las_file)
  dst_ext <- union_extent_from_headers(norm_files)
  ar <- extent_area(dst_ext) / extent_area(src_ext)
  if (is.finite(ar) && ar >= COVERAGE_OK) {
    .mem_snap(sprintf("Step 2: reuse normalized chunks (coverage %.1f%%)", 100*ar))
    need_normalize <- FALSE
  } else {
    .mem_snap(sprintf("Step 2: normalized chunks incomplete (coverage %.1f%%) -> recomputing", 100*ar))
    unlink(norm_dir, recursive=TRUE, force=TRUE)
    dir.create(norm_dir, recursive=TRUE, showWarnings=FALSE)
  }
} else if (length(norm_files) > 0L && !file.exists(norm_done)) {
  .mem_snap("Step 2: normalized chunks present without sentinel -> recomputing")
  unlink(norm_dir, recursive=TRUE, force=TRUE)
  dir.create(norm_dir, recursive=TRUE, showWarnings=FALSE)
}

if (need_normalize) {
  .mem_snap("Step 2: normalize_height (chunked) start")
  ctg_norm <- readLAScatalog(las_file)
  opt_chunk_size(ctg_norm)      <- NORM_CHUNK_M
  opt_chunk_buffer(ctg_norm)    <- NORM_BUFFER_M
  opt_progress(ctg_norm)        <- TRUE
  opt_output_files(ctg_norm)    <- file.path(norm_dir, "{XLEFT}_{YBOTTOM}_{ID}")
  opt_select(ctg_norm)          <- "*"   # preserve attributes
  opt_laz_compression(ctg_norm) <- (NORM_LAZ == 1L)
  ctg_norm <- catalog_map(ctg_norm, function(l, dtm) {
    if (lidR::is.empty(l)) return(NULL)
    xr <- range(l@data$X, na.rm = TRUE)
    yr <- range(l@data$Y, na.rm = TRUE)
    .mem_snap(sprintf("Step 2 chunk bbox: xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f", xr[1], xr[2], yr[1], yr[2]))

    dtm_res <- suppressWarnings(raster::res(dtm))
    pad <- max(c(NORM_BUFFER_M, dtm_res * 2), na.rm = TRUE)
    if (!is.finite(pad)) pad <- NORM_BUFFER_M

    dtm_extent <- raster::extent(
      xr[1] - pad,
      xr[2] + pad,
      yr[1] - pad,
      yr[2] + pad
    )
    dtm_sub <- tryCatch(raster::crop(dtm, dtm_extent), error = function(e) NULL)
    if (is.null(dtm_sub) || raster::ncell(dtm_sub) == 0) {
      stop(sprintf("DTM crop is empty for padded chunk bbox xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f pad=%.3f", xr[1], xr[2], yr[1], yr[2], pad))
    }
    dtm_vals <- raster::getValues(dtm_sub)
    if (!length(dtm_vals) || all(is.na(dtm_vals))) {
      stop(sprintf("DTM crop is all NA for padded chunk bbox xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f pad=%.3f", xr[1], xr[2], yr[1], yr[2], pad))
    }

    z_dtm <- tryCatch(raster::extract(dtm, cbind(l@data$X, l@data$Y)), error = function(e) NULL)
    if (is.null(z_dtm) || !length(z_dtm)) {
      stop(sprintf("DTM point extraction failed for chunk bbox xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f", xr[1], xr[2], yr[1], yr[2]))
    }
    point_coverage <- mean(is.finite(z_dtm))
    if (!is.finite(point_coverage) || point_coverage < COVERAGE_OK) {
      stop(sprintf(
        "DTM point coverage %.1f%% for chunk bbox xmin=%.3f xmax=%.3f ymin=%.3f ymax=%.3f (< %.1f%%).",
        100 * point_coverage, xr[1], xr[2], yr[1], yr[2], 100 * COVERAGE_OK
      ))
    }

    l <- normalize_height(l, dtm)
    if (!is.null(l$Z)) {
      bad <- l$Z < 0
      if (any(bad, na.rm=TRUE)) l$Z[bad] <- 0
    }
    l
  }, dtm = dtm)
  gc()

  # Validate normalized chunk coverage against the clipped/source LAS BEFORE
  # writing the sentinel. This catches silent tile loss during catalog_map().
  norm_files_now <- list.files(norm_dir, pattern="\\.(las|laz)$", full.names=TRUE)
  if (!length(norm_files_now)) {
    stop("Normalization produced no output tiles in: ", norm_dir)
  }
  src_ext_now <- las_header_extent(las_file)
  dst_ext_now <- union_extent_from_headers(norm_files_now)
  ar_now <- extent_area(dst_ext_now) / extent_area(src_ext_now)
  if (!is.finite(ar_now) || ar_now < COVERAGE_OK) {
    stop(sprintf(
      "Normalized chunk coverage only %.1f%% of source (< %.1f%%).",
      100 * ar_now, 100 * COVERAGE_OK
    ))
  }

  # create sentinel AFTER successful validation
  file.create(norm_done)
  .mem_snap("Step 2: normalize_height (chunked) done + sentinel written")
}

# Free DTM before retile
rm(dtm); gc()
.mem_snap("Step 2: freed DTM before retile")

# ---------- MERGE normalized chunks → single site-wide file (HARDENED) ----------
final_out <- file.path(processed_folder_path, base_name)

# Optional manual override to force rebuild of the final merged file
if (FORCE_REMERGE == 1L && file.exists(final_out)) {
  .mem_snap("Step 2: FORCE_REMERGE=1 -> removing existing final_out")
  try(file.remove(final_out), silent = TRUE)
}

# Normalize tile set must exist
norm_files <- list.files(norm_dir, pattern="\\.(las|laz)$", full.names=TRUE)
if (!length(norm_files)) stop("No normalized tiles found in: ", norm_dir)

# Compute both source and normalized extents.
# Use the clipped/source LAS as the authoritative site extent so missing normalized
# edge tiles are detected instead of being treated as the expected footprint.
src_ext   <- las_header_extent(las_file)
norm_ext  <- union_extent_from_headers(norm_files)
width     <- src_ext[2] - src_ext[1]
height    <- src_ext[4] - src_ext[3]
site_area <- extent_area(src_ext)

# Decide whether we need to (re)merge
need_remerge <- TRUE
if (file.exists(final_out)) {
  fin_ext <- las_header_extent(final_out)
  ar <- extent_area(fin_ext) / site_area
  if (is.finite(ar) && ar >= COVERAGE_OK) {
    .mem_snap(sprintf("Step 2: final normalized file seems complete (coverage %.1f%%) — keeping", 100*ar))
    need_remerge <- FALSE
  } else {
    .mem_snap(sprintf("Step 2: final normalized file incomplete (coverage %.1f%%) — remerging", 100*ar))
    try(file.remove(final_out), silent = TRUE)
  }
}

if (need_remerge) {
  # Clean any stale merged_* candidates from prior runs
  old_merged <- list.files(
    processed_folder_path,
    pattern = paste0("^", stem, "_merged_.*\\.(las|laz)$"),
    full.names = TRUE
  )
  if (length(old_merged)) try(file.remove(old_merged), silent = TRUE)
  
  # We'll try a few times with progressively larger tiles if needed
  attempt <- 1L
  max_attempts <- 3L
  merged_candidates <- character(0)
  
  repeat {
    # Build a catalog from the normalized tiles
    ctg_norm2 <- readLAScatalog(norm_dir)
    opt_progress(ctg_norm2)        <- TRUE
    opt_laz_compression(ctg_norm2) <- write_laz          # match input compression for the final file
    opt_select(ctg_norm2)          <- "*"                # keep all attributes
    opt_chunk_buffer(ctg_norm2)    <- 0                  # output has no buffer
    
    # Critical bits: size > site span and align grid on the site's lower-left corner
    fudge <- c(1, 10, 100)[attempt]                      # grow margin if needed
    tile_size <- max(width, height) + fudge
    opt_chunk_size(ctg_norm2)       <- tile_size
    # If function is available, align; otherwise we still rely on a bigger tile
    try({ opt_chunk_alignment(ctg_norm2) <- c(norm_ext[1], norm_ext[3]) }, silent = TRUE)
    
    tmp_prefix <- file.path(processed_folder_path, paste0(stem, "_merged"))
    opt_output_files(ctg_norm2)     <- paste0(tmp_prefix, "_{XLEFT}_{YBOTTOM}")
    
    .mem_snap(sprintf(
      "Step 2: catalog_retile attempt %d (tile=%.3f, align=%.3f/%.3f)",
      attempt, tile_size, src_ext[1], src_ext[3]
    ))
    invisible(catalog_retile(ctg_norm2))
    gc()
    
    merged_candidates <- list.files(
      processed_folder_path,
      pattern = paste0("^", stem, "_merged_.*\\.(las|laz)$"),
      full.names = TRUE
    )
    
    if (length(merged_candidates) == 1L) break
    if (attempt >= max_attempts) {
      stop(sprintf(
        "Merge produced %d tiles after %d attempts. Increase tile size or inspect inputs.",
        length(merged_candidates), max_attempts
      ))
    }
    
    # Clean up and retry with a larger tile (and same alignment)
    if (length(merged_candidates)) try(file.remove(merged_candidates), silent = TRUE)
    attempt <- attempt + 1L
  }
  
  # Move the single merged candidate to the final destination name
  if (file.exists(final_out)) file.remove(final_out)
  file.rename(merged_candidates[[1]], final_out)
  
  # Final coverage validation against the tile set union extent
  fin_ext <- las_header_extent(final_out)
  ar2 <- extent_area(fin_ext) / site_area
  if (!is.finite(ar2) || ar2 < COVERAGE_OK) {
    stop(sprintf(
      "Final normalized file covers only %.1f%% of the site (< %.1f%%).",
      100*ar2, 100*COVERAGE_OK
    ))
  }
  .mem_snap(paste("Step 2: normalized single file ready:", final_out))
} else {
  .mem_snap(paste("Step 2: final normalized file already complete:", final_out))
}

# =========================================================
# STEP 2B: Apply optional hard height ceiling, then classify/remove noise
# =========================================================
.mem_snap("Step 2B: read normalized LAS for height/noise filtering")

las_clean <- readLAS(final_out, select = "*")
if (is.null(las_clean) || lidR::is.empty(las_clean)) {
  stop("Normalized LAS is empty before height/noise filtering: ", final_out)
}

if (!is.na(max_normalized_height)) {
  .mem_snap(sprintf(
    "Step 2B: hard normalized-height cutoff start (Z <= %.6f m)",
    max_normalized_height
  ))

  n_before_height <- lidR::npoints(las_clean)
  las_clean <- filter_poi(las_clean, is.finite(Z) & Z <= max_normalized_height)
  n_after_height <- lidR::npoints(las_clean)

  if (is.null(las_clean) || lidR::is.empty(las_clean)) {
    stop(sprintf(
      "Hard normalized-height cutoff Z <= %.6f m removed every point from: %s",
      max_normalized_height,
      final_out
    ))
  }

  .mem_snap(sprintf(
    "Step 2B: height cutoff removed %d points above %.6f m; %d points kept",
    n_before_height - n_after_height,
    max_normalized_height,
    n_after_height
  ))
} else {
  .mem_snap("Step 2B: hard normalized-height cutoff disabled")
}

.mem_snap("Step 2B: classify_noise(SOR) start")

las_clean <- classify_noise(
  las_clean,
  algorithm = sor(k = 15, m = 5)
)

.mem_snap("Step 2B: filter LASNOISE start")

n_before <- lidR::npoints(las_clean)
las_clean <- filter_poi(las_clean, Classification != LASNOISE)
n_after <- lidR::npoints(las_clean)

.mem_snap(sprintf(
  "Step 2B: noise removed: %d points removed; %d points kept",
  n_before - n_after,
  n_after
))

tmp_clean_out <- file.path(
  processed_folder_path,
  paste0(stem, "_noise_cleaned.", ext)
)

if (file.exists(tmp_clean_out)) file.remove(tmp_clean_out)
writeLAS(las_clean, tmp_clean_out)

rm(las_clean)
gc()

if (file.exists(final_out)) file.remove(final_out)
file.rename(tmp_clean_out, final_out)

.mem_snap(paste("Step 2B: cleaned normalized LAS written:", final_out))

# =========================================================
# STEP 3: CHM from merged normalized LAS (single-file mode)
# =========================================================
.mem_snap("Step 3: read merged normalized LAS for CHM")
las_norm <- readLAS(final_out, select = "xyz")
if (is.null(las_norm) || lidR::is.empty(las_norm)) {
  stop("Merged normalized LAS is empty: ", final_out)
}

.mem_snap("Step 3: rasterize_canopy(p2r) start")
chm <- rasterize_canopy(las_norm, res = 1, algorithm = p2r())
gc()
.mem_snap("Step 3: rasterize_canopy(p2r) done")

dtm_template <- raster(dtm_tif)
chm <- raster::raster(chm)

if (!raster::compareRaster(
  chm, dtm_template,
  extent = TRUE, rowcol = TRUE, crs = FALSE, res = TRUE, orig = TRUE,
  stopiffalse = FALSE
)) {
  .mem_snap("Step 3: aligning CHM to DTM grid")
  chm <- raster::extend(chm, raster::extent(dtm_template))
  if (!raster::compareRaster(
    chm, dtm_template,
    extent = TRUE, rowcol = TRUE, crs = FALSE, res = TRUE, orig = TRUE,
    stopiffalse = FALSE
  )) {
    chm <- raster::resample(chm, dtm_template, method = "bilinear")
  }
}

raster::crs(chm) <- raster::crs(dtm_template)
chm <- raster::mask(chm, dtm_template)

chm_vals <- raster::getValues(chm)
if (!any(is.finite(chm_vals))) {
  stop("CHM is all NA after rasterization/alignment.")
}

chm_tif <- file.path(results_folder_path, "CHM.tif")
if (file.exists(chm_tif)) file.remove(chm_tif)
writeRaster(chm, chm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 3: CHM written", chm_tif))

png(file.path(results_folder_path, "CHM.png"), width = 1000, height = 800)
plot(chm, main = "Canopy Height Model (CHM)")
dev.off()
.mem_snap("Step 3: CHM preview PNG")

# ---------- CHM height histogram ----------
chm_hist_png <- file.path(results_folder_path, "CHM_hist.png")
tryCatch({
  vals <- raster::getValues(chm)
  vals <- vals[is.finite(vals)]

  if (length(vals) > 0) {
    png(chm_hist_png, width = 1000, height = 800)
    hist(
      vals,
      main = "CHM Height Histogram",
      xlab = "Height (m)",
      ylab = "Frequency"
    )
    dev.off()
    .mem_snap("Step 3: CHM histogram PNG")
  } else {
    .mem_snap("Step 3: CHM histogram skipped (no finite values)")
  }
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  .mem_snap(paste("Step 3: CHM histogram skipped due to error:", conditionMessage(e)))
})

# ---------- Additional 0.5 m CHM aligned to DTM_05 ----------
.mem_snap("Step 3: rasterize_canopy(p2r) 0.5m start")
chm05 <- rasterize_canopy(las_norm, res = 0.5, algorithm = p2r())
gc()
.mem_snap("Step 3: rasterize_canopy(p2r) 0.5m done")

dtm05_chm_template <- raster(dtm05_tif)
chm05 <- raster::raster(chm05)

if (!raster::compareRaster(
  chm05, dtm05_chm_template,
  extent = TRUE, rowcol = TRUE, crs = FALSE, res = TRUE, orig = TRUE,
  stopiffalse = FALSE
)) {
  .mem_snap("Step 3: aligning 0.5m CHM to DTM_05 grid")
  chm05 <- raster::extend(chm05, raster::extent(dtm05_chm_template))
  if (!raster::compareRaster(
    chm05, dtm05_chm_template,
    extent = TRUE, rowcol = TRUE, crs = FALSE, res = TRUE, orig = TRUE,
    stopiffalse = FALSE
  )) {
    chm05 <- raster::resample(chm05, dtm05_chm_template, method = "bilinear")
  }
}

raster::crs(chm05) <- raster::crs(dtm05_chm_template)
chm05 <- raster::mask(chm05, dtm05_chm_template)

chm05_vals <- raster::getValues(chm05)
if (!any(is.finite(chm05_vals))) {
  stop("CHM_05 is all NA after rasterization/alignment.")
}

chm05_tif <- file.path(results_folder_path, "CHM_05.tif")
if (file.exists(chm05_tif)) file.remove(chm05_tif)
writeRaster(chm05, chm05_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 3: CHM_05 written", chm05_tif))

png(file.path(results_folder_path, "CHM_05.png"), width = 1000, height = 800)
plot(chm05, main = "Canopy Height Model 0.5 m (CHM_05)")
dev.off()
.mem_snap("Step 3: CHM_05 preview PNG")

.mem_snap("Processing complete")
