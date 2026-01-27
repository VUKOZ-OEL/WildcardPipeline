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
if (length(args) < 4) {
  stop("Usage: Rscript preprocess_single_las.R <las_folder_path> <processed_folder_path> <results_folder_path> <polygon_file> <coordinates>")
}
las_folder_path       <- args[1]
processed_folder_path <- args[2]
results_folder_path   <- args[3]
coordinates           <- args[5]

cat("Corrdinates: ", coordinates, "\n")

# Optional polygon argument
polygon_file <- if (length(args) >= 4) args[4] else NA_character_
use_polygon  <- !is.na(polygon_file) && nzchar(polygon_file)
if (use_polygon && !file.exists(polygon_file)) {
  stop("Polygon file does not exist: ", polygon_file)
}

# Fallback CRS if LAS has none (same behaviour as your other script)
crs_code <- coordinates

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

# ---------- Clip LAS to polygon at the beginning (if requested) ----------
if (use_polygon) {
  .mem_snap("Polygon clip: reading LAS and polygon")
  
  las_full <- readLAS(las_file)
  if (is.null(las_full) || lidR::is.empty(las_full)) {
    stop("Input LAS is empty before clipping: ", las_file)
  }
  
  # Get or set CRS on the LAS using lidR::st_crs
  las_crs_obj <- tryCatch(lidR::st_crs(las_full), error = function(e) NULL)
  if (is.null(las_crs_obj) || is.na(las_crs_obj)) {
    lidR::st_crs(las_full) <- crs_code   # set EPSG:5514 as CRS
    las_crs_obj <- lidR::st_crs(las_full)
    .mem_snap(paste("Polygon clip: LAS had no CRS, setting fallback:", crs_code))
  }
  
  poly <- sf::st_read(polygon_file, quiet = TRUE)
  if (is.null(sf::st_crs(poly))) stop("Polygon has no CRS; please define one.")
  
  # Clean geometry and reproject polygon to LAS CRS
  poly <- sf::st_make_valid(poly)
  suppressWarnings({
    poly <- sf::st_zm(poly, drop = TRUE, what = "ZM")  # drop Z/M if present
  })
  poly <- sf::st_transform(poly, crs = las_crs_obj)
  poly <- suppressWarnings(sf::st_union(poly))
  
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
CLASS_BUFFER_M <- as.numeric(Sys.getenv("CLASS_BUFFER_M", "30"))
NORM_CHUNK_M   <- as.numeric(Sys.getenv("NORM_CHUNK_M",   "80"))
NORM_BUFFER_M  <- as.numeric(Sys.getenv("NORM_BUFFER_M",  "2"))
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
    .mem_snap(paste("Step 1:", label, "small-hole fill start"))
    dtm_filled <- raster::focal(
      dtm_raster, w = matrix(1, 3, 3), pad = TRUE,
      fun = function(x, ...) {
        center <- x[ceiling(length(x) / 2)]
        if (is.na(center)) mean(x, na.rm = TRUE) else center
      }
    )
    dtm_raster <- raster::cover(dtm_raster, dtm_filled)
    na_cells <- raster::cellStats(is.na(dtm_raster), stat = "sum")
    .mem_snap(sprintf("%s NA cells after fill: %d", label, as.integer(na_cells)))
    .mem_snap(paste("Step 1:", label, "small-hole fill done"))
  } else {
    na_cells <- raster::cellStats(is.na(dtm_raster), stat = "sum")
    .mem_snap(sprintf("%s NA cells (no fill): %d", label, as.integer(na_cells)))
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

# Optional: small-hole fill for isolated NA patches (1 m DTM, used later)
dtm <- fill_dtm_if_needed(dtm, "DTM")

# Save DTM and reopen file-backed (this DTM is used for normalization)
dtm_tif <- file.path(results_folder_path, "DTM.tif")
if (file.exists(dtm_tif)) file.remove(dtm_tif)
writeRaster(dtm, dtm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
dtm <- raster(dtm_tif)
.mem_snap(paste("Step 1: DTM written", dtm_tif))

png(file.path(results_folder_path, "DTM.png"), width = 1000, height = 800)
plot(dtm, main = "Digital Terrain Model (DTM)")
dev.off()
.mem_snap("Step 1: DTM preview PNG")

# ---------- Additional 0.5 m DTM (output only) ----------
.mem_snap("Step 1: DTM_05 grid_terrain(knnidw) start")
if (requireNamespace("terra", quietly = TRUE)) {
  terra::tmpFiles(current = TRUE, old = TRUE, orphan = TRUE, remove = TRUE)
}

if (!is.null(ctg_class@output_options$drivers$SpatRaster)) {
  ctg_class@output_options$drivers$SpatRaster$param$overwrite <- TRUE
}
.mem_snap("Temporary terra files removed")
ctg_class <- set_raster_chunk_output(ctg_class, file.path(work_dir, "dtm05m_tiles"), "DTM05m")
dtm05 <- grid_terrain(
  ctg_class,
  res       = 0.5,
  algorithm = knnidw(k = 10, p = 2),
  overwrite = TRUE
)

gc()
.mem_snap("Step 1: DTM_05 grid_terrain(knnidw) done")

dtm05_tmp <- file.path(tmpdir, "DTM_05_tmp.tif")

# For terra backend:
if (inherits(dtm05, "SpatRaster")) {
  # remove any existing files with that base name
  if (file.exists(dtm05_tmp)) terra::rmf(dtm05_tmp)  # removes tif + aux files
  terra::writeRaster(dtm05, dtm05_tmp, overwrite = TRUE)
}

# For raster backend:
if (inherits(dtm05, "Raster")) {
  if (file.exists(dtm05_tmp)) unlink(paste0(tools::file_path_sans_ext(dtm05_tmp), ".*"),
                                     recursive = TRUE, force = TRUE)
  raster::writeRaster(dtm05, dtm05_tmp, overwrite = TRUE)
}

# Small-hole fill for isolated NA patches (0.5 m DTM, output only)
dtm05 <- fill_dtm_if_needed(dtm05, "DTM_05")

# Save 0.5 m DTM (no further use downstream)
dtm05_tif <- file.path(results_folder_path, "DTM_05.tif")
if (file.exists(dtm05_tif)) file.remove(dtm05_tif)
writeRaster(dtm05, dtm05_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 1: DTM_05 written", dtm05_tif))

# DTM_05 preview as PDF
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
    l <- normalize_height(l, dtm)
    if (!is.null(l$Z)) {
      bad <- l$Z < 0
      if (any(bad, na.rm=TRUE)) l$Z[bad] <- 0
    }
    l
  }, dtm = dtm)
  gc()
  # create sentinel AFTER success
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

# Compute the union extent of the normalized tiles (expected coverage)
norm_ext  <- union_extent_from_headers(norm_files)
width     <- norm_ext[2] - norm_ext[1]
height    <- norm_ext[4] - norm_ext[3]
site_area <- extent_area(norm_ext)

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
      attempt, tile_size, norm_ext[1], norm_ext[3]
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
# ---------- END HARDENED MERGE block ----------

# =========================================================
# STEP 3: CHM from normalized chunks — CHUNKED raster
# =========================================================
ctg_norm_for_chm <- readLAScatalog(norm_dir)
opt_chunk_size(ctg_norm_for_chm)   <- CHM_CHUNK_M
opt_chunk_buffer(ctg_norm_for_chm) <- CHM_BUFFER_M
opt_progress(ctg_norm_for_chm)     <- TRUE
opt_select(ctg_norm_for_chm)       <- "xyz"
.mem_snap("Step 3: grid_canopy(p2r) start")
chm <- grid_canopy(ctg_norm_for_chm, res = 1, algorithm = p2r(subcircle = 0.2))
gc()
.mem_snap("Step 3: grid_canopy(p2r) done")

chm_tif <- file.path(results_folder_path, "CHM.tif")
if (file.exists(chm_tif)) file.remove(chm_tif)
writeRaster(chm, chm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))
.mem_snap(paste("Step 3: CHM written", chm_tif))

png(file.path(results_folder_path, "CHM.png"), width = 1000, height = 800)
plot(chm, main = "Canopy Height Model (CHM)")
dev.off()
.mem_snap("Step 3: CHM preview PNG")

# ---------- CHM height histogram ----------
png(file.path(results_folder_path, "CHM_hist.png"), width = 1000, height = 800)
raster::hist(
  chm,
  main = "CHM Height Histogram",
  xlab = "Height (m)",
  ylab = "Frequency"
)
dev.off()
.mem_snap("Step 3: CHM histogram PNG")

# ---------- Cleanup work dir (always remove .work/<stem>) ----------
if (dir.exists(work_dir)) {
  try(unlink(work_dir, recursive = TRUE, force = TRUE), silent = TRUE)
  work_parent <- dirname(work_dir)  # results_folder_path/.work
  if (dir.exists(work_parent) && length(list.files(work_parent)) == 0L) {
    try(unlink(work_parent, recursive = TRUE, force = TRUE), silent = TRUE)
  }
}

.mem_snap("Processing complete")
