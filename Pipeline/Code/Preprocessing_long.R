suppressPackageStartupMessages({

  library(lidR)

  library(raster)

  library(ggplot2)

})



# ---------- Deterministic & single-thread ----------

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



# ---------- Force sequential if {future} exists ----------

.old_plan <- NULL

if (requireNamespace("future", quietly = TRUE)) {

  .old_plan <- future::plan()

  future::plan("sequential")

}

.onexit_future <- function() {

  if (!is.null(.old_plan)) try(future::plan(.old_plan), silent = TRUE)

}

on.exit(.onexit_future(), add = TRUE)



# ---------- Memory logger ----------

.memlog <- data.frame(step = character(), rss_MiB = numeric(), r_MiB = numeric(), stringsAsFactors = FALSE)

.as_mb <- function(x) if (is.na(x)) "NA" else sprintf("%.2f", x)

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



# Also dump last warnings

on.exit({

  ws <- warnings()

  if (!is.null(ws)) {

    cat("\n===== LAST WARNINGS (max 50) =====\n")

    print(utils::tail(ws, 50))

    cat("===================================\n")

  }

}, add = TRUE)



# ---------- Args ----------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {

  stop("Usage: Rscript preprocess_single_las.R <las_folder_path> <processed_folder_path> <results_folder_path>")

}

las_folder_path       <- args[1]

processed_folder_path <- args[2]

results_folder_path   <- args[3]



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



# ---------- Work dirs ----------

work_dir  <- file.path(results_folder_path, ".work", stem)  # ephemeral

RESUME    <- as.integer(Sys.getenv("RESUME", "1"))

if (RESUME == 0L && dir.exists(work_dir)) unlink(work_dir, recursive = TRUE, force = TRUE)

dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)



class_dir    <- file.path(work_dir, "class")     # classified chunks (ephemeral)

norm_tmp_dir <- file.path(work_dir, "norm_tmp")  # normalized tiles computed here, then moved

dir.create(class_dir,    showWarnings = FALSE, recursive = TRUE)

dir.create(norm_tmp_dir, showWarnings = FALSE, recursive = TRUE)



# ---------- Single-layer processed dataset dir ----------

# If processed_folder_path already ends with <stem>, use it directly; else use <processed>/<stem>.

# You can force flat behavior with USE_FLAT_DIR=1

processed_base <- tryCatch(basename(normalizePath(processed_folder_path, mustWork = FALSE)), error = function(e) NA_character_)

USE_FLAT_DIR   <- as.integer(Sys.getenv("USE_FLAT_DIR", ifelse(!is.na(processed_base) && processed_base == stem, "1", "0")))

dataset_dir    <- if (USE_FLAT_DIR == 1L) processed_folder_path else file.path(processed_folder_path, stem)

dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

.mem_snap(paste("Normalized tiles will be written to:", dataset_dir))



# ---------- Raster tmpdir ----------

run_tag <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", as.integer(Sys.getpid()))

tmpdir  <- file.path(work_dir, paste0("raster_tmp_", run_tag))

dir.create(tmpdir, recursive = TRUE, showWarnings = FALSE)

raster::rasterOptions(tmpdir = tmpdir, overwrite = TRUE, progress = "text")

.mem_snap(paste("Raster tmpdir:", tmpdir))



# ---------- Tunables (env overrides) ----------

CLASS_CHUNK_M  <- as.numeric(Sys.getenv("CLASS_CHUNK_M",  "300"))

CLASS_BUFFER_M <- as.numeric(Sys.getenv("CLASS_BUFFER_M", "30"))

NORM_CHUNK_M   <- as.numeric(Sys.getenv("NORM_CHUNK_M",   "80"))

NORM_BUFFER_M  <- as.numeric(Sys.getenv("NORM_BUFFER_M",  "2"))

CHM_CHUNK_M    <- as.numeric(Sys.getenv("CHM_CHUNK_M",    "300"))

CHM_BUFFER_M   <- as.numeric(Sys.getenv("CHM_BUFFER_M",   "30"))



NORM_LAZ       <- as.integer(Sys.getenv("NORM_LAZ", "0"))       # 0=LAS faster, 1=LAZ smaller

DTM_FILL_NA    <- as.integer(Sys.getenv("DTM_FILL_NA", "1"))    # fill isolated NAs in DTM

COVERAGE_OK    <- as.numeric(Sys.getenv("COVERAGE_OK", "0.98")) # reuse/accept threshold



.mem_snap(sprintf(

  "Chunk params: CLASS=%dm/%dm  NORM=%dm/%dm  CHM=%dm/%dm  NORM_LAZ=%d  USE_FLAT_DIR=%d",

  CLASS_CHUNK_M, CLASS_BUFFER_M, NORM_CHUNK_M, NORM_BUFFER_M, CHM_CHUNK_M, CHM_BUFFER_M, NORM_LAZ, USE_FLAT_DIR

))



# ---------- Helpers ----------

ensure_clean_or_resume <- function(dir, expected_row, recursive = TRUE) {

  meta <- file.path(dir, "_meta.tsv")

  have <- length(list.files(dir, pattern="\\.(las|laz)$", full.names=TRUE)) > 0L

  ok <- FALSE

  if (file.exists(meta)) {

    old <- tryCatch(read.delim(meta, header=TRUE, sep="\t", stringsAsFactors=FALSE), error=function(e) NULL)

    ok <- !is.null(old) && identical(colnames(old), colnames(expected_row)) &&

          all(as.character(old[1,]) == as.character(expected_row[1,]))

  }

  if (have && !ok) {

    if (recursive) {

      unlink(dir, recursive=TRUE, force=TRUE)

      dir.create(dir, recursive=TRUE, showWarnings=FALSE)

    } else {

      files <- list.files(dir, pattern="\\.(las|laz)$", full.names=TRUE)

      if (length(files)) try(file.remove(files), silent = TRUE)

      if (file.exists(meta)) try(file.remove(meta), silent = TRUE)

    }

  }

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)

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



# Meta guards: class tiles (work), normalized tiles (deliverable)

ensure_clean_or_resume(class_dir,  data.frame(kind="class", chunk=CLASS_CHUNK_M, buffer=CLASS_BUFFER_M), recursive = TRUE)

ensure_clean_or_resume(dataset_dir, data.frame(kind="norm",  chunk=NORM_CHUNK_M,  buffer=NORM_BUFFER_M), recursive = FALSE)



# Robust move (works across filesystems)

.move_files <- function(src_files, dst_dir) {

  dest <- file.path(dst_dir, basename(src_files))

  ok <- file.rename(src_files, dest)

  if (any(!ok)) {

    to_copy <- src_files[!ok]

    ok2 <- file.copy(to_copy, dst_dir, overwrite = TRUE)

    if (all(ok2)) unlink(to_copy)

    ok[!ok] <- ok2

  }

  all(ok)

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



# DTM from ground only

opt_filter(ctg_class)      <- "-keep_class 2"

opt_select(ctg_class)      <- "xyz"

.mem_snap("Step 1: DTM prep (filter=ground, select=xyz)")



.mem_snap("Step 1: grid_terrain(knnidw) start")

dtm <- grid_terrain(ctg_class, res = 1, algorithm = knnidw(k = 10, p = 2))

gc(); .mem_snap("Step 1: grid_terrain(knnidw) done")



# Optional small-hole fill

if (DTM_FILL_NA == 1L) {

  .mem_snap("Step 1: DTM small-hole fill start")

  dtm_filled <- raster::focal(

    dtm, w = matrix(1, 3, 3), pad = TRUE,

    fun = function(x, ...) {

      center <- x[ceiling(length(x) / 2)]

      if (is.na(center)) mean(x, na.rm = TRUE) else center

    }

  )

  dtm <- raster::cover(dtm, dtm_filled)

  na_cells <- raster::cellStats(is.na(dtm), stat = "sum")

  .mem_snap(sprintf("DTM NA cells after fill: %d", as.integer(na_cells)))

  .mem_snap("Step 1: DTM small-hole fill done")

} else {

  na_cells <- raster::cellStats(is.na(dtm), stat = "sum")

  .mem_snap(sprintf("DTM NA cells (no fill): %d", as.integer(na_cells)))

}



# Persist DTM

dtm_tif <- file.path(results_folder_path, "DTM.tif")

if (file.exists(dtm_tif)) file.remove(dtm_tif)

writeRaster(dtm, dtm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))

dtm <- raster(dtm_tif)

.mem_snap(paste("Step 1: DTM written", dtm_tif))



png(file.path(results_folder_path, "DTM.png"), width = 1000, height = 800)

plot(dtm, main = "Digital Terrain Model (DTM)")

dev.off()

.mem_snap("Step 1: DTM preview PNG")



# =========================================================

# STEP 2: Normalize — write tiles into single-layer processed dir (no merging)

# =========================================================

norm_done     <- file.path(dataset_dir, "_DONE")  # sentinel in deliverable dir

existing_norm <- list.files(dataset_dir, pattern="\\.(las|laz)$", full.names=TRUE)



need_normalize <- TRUE

if (length(existing_norm) > 0L && file.exists(norm_done)) {

  src_ext <- las_header_extent(las_file)

  dst_ext <- union_extent_from_headers(existing_norm)

  ar <- extent_area(dst_ext) / extent_area(src_ext)

  if (is.finite(ar) && ar >= COVERAGE_OK) {

    .mem_snap(sprintf("Step 2: reuse normalized tiles (coverage %.1f%%)", 100*ar))

    need_normalize <- FALSE

  } else {

    .mem_snap(sprintf("Step 2: normalized tiles incomplete (coverage %.1f%%) -> recomputing", 100*ar))

    files <- list.files(dataset_dir, pattern="\\.(las|laz)$", full.names=TRUE)

    if (length(files)) try(file.remove(files), silent = TRUE)

    if (file.exists(norm_done)) try(file.remove(norm_done), silent = TRUE)

  }

} else if (length(existing_norm) > 0L && !file.exists(norm_done)) {

  .mem_snap("Step 2: tiles present without sentinel -> recomputing")

  files <- list.files(dataset_dir, pattern="\\.(las|laz)$", full.names=TRUE)

  if (length(files)) try(file.remove(files), silent = TRUE)

}



if (need_normalize) {

  # Ensure temp dir is empty

  if (dir.exists(norm_tmp_dir)) unlink(norm_tmp_dir, recursive = TRUE, force = TRUE)

  dir.create(norm_tmp_dir, showWarnings = FALSE, recursive = TRUE)



  .mem_snap("Step 2: normalize_height (chunked) start")

  ctg_norm <- readLAScatalog(las_file)

  opt_chunk_size(ctg_norm)      <- NORM_CHUNK_M

  opt_chunk_buffer(ctg_norm)    <- NORM_BUFFER_M

  opt_progress(ctg_norm)        <- TRUE

  opt_output_files(ctg_norm)    <- file.path(norm_tmp_dir, "{XLEFT}_{YBOTTOM}_{ID}")  # write to temp

  opt_select(ctg_norm)          <- "*"   # keep attributes

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



  gc(); .mem_snap("Step 2: normalize_height (chunked) done (temp folder)")



  # Validate coverage before exposing tiles

  tmp_tiles <- list.files(norm_tmp_dir, pattern="\\.(las|laz)$", full.names=TRUE)

  if (!length(tmp_tiles)) stop("Normalization produced no tiles.")

  src_ext <- las_header_extent(las_file)

  tmp_ext <- union_extent_from_headers(tmp_tiles)

  ar_tmp  <- extent_area(tmp_ext) / extent_area(src_ext)

  .mem_snap(sprintf("Step 2: normalized temp coverage: %.1f%%", 100*ar_tmp))

  if (!is.finite(ar_tmp) || ar_tmp < COVERAGE_OK) {

    stop(sprintf("Normalized tiles cover only %.1f%% of source (< %.1f%%).", 100*ar_tmp, 100*COVERAGE_OK))

  }



  # Move tiles atomically into the single-layer processed dir

  .mem_snap("Step 2: moving normalized tiles to processed dir")

  ok <- .move_files(tmp_tiles, dataset_dir)

  if (!ok) stop("Failed to move all normalized tiles into processed dir.")



  # Drop temp folder and write sentinel

  unlink(norm_tmp_dir, recursive = TRUE, force = TRUE)

  file.create(norm_done)

  .mem_snap("Step 2: normalized tiles ready + sentinel written")

}



# Free DTM as we won't need it beyond here

rm(dtm); gc()

.mem_snap("Step 2: freed DTM")



# =========================================================

# STEP 3: CHM from normalized tiles (in single-layer processed dir)

# =========================================================

ctg_norm_for_chm <- readLAScatalog(dataset_dir)

opt_chunk_size(ctg_norm_for_chm)   <- CHM_CHUNK_M

opt_chunk_buffer(ctg_norm_for_chm) <- CHM_BUFFER_M

opt_progress(ctg_norm_for_chm)     <- TRUE

opt_select(ctg_norm_for_chm)       <- "xyz"



.mem_snap("Step 3: grid_canopy(p2r) start")

chm <- grid_canopy(ctg_norm_for_chm, res = 1, algorithm = p2r(subcircle = 0.2))

gc(); .mem_snap("Step 3: grid_canopy(p2r) done")



chm_tif <- file.path(results_folder_path, "CHM.tif")

if (file.exists(chm_tif)) file.remove(chm_tif)

writeRaster(chm, chm_tif, overwrite = TRUE, datatype = "FLT4S", options = c("COMPRESS=LZW"))

.mem_snap(paste("Step 3: CHM written", chm_tif))



png(file.path(results_folder_path, "CHM.png"), width = 1000, height = 800)

plot(chm, main = "Canopy Height Model (CHM)")

dev.off()

.mem_snap("Step 3: CHM preview PNG")



# ---------- Optional cleanup ----------

if (RESUME == 0L) try(unlink(work_dir, recursive = TRUE), silent = TRUE)



.mem_snap("Processing complete")

