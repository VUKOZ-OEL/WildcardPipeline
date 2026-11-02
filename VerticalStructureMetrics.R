suppressMessages(library(lidR))
suppressMessages(library(raster))   
suppressMessages(library(data.table))
suppressMessages(library(sf))

# ---------------------- CLI ----------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript VSC_UVC_Metrics.R <input_folder> <output_folder> <max_height_file> <output_prefix> <LOCAL|EPSG:code>")
}
las_folder <- args[1]
output_folder <- args[2]
max_height_file  <- args[3]   
base_name <- args[4]
coordinates <- trimws(args[5])

if (is.na(coordinates) || !nzchar(coordinates)) {
  stop("coordinates must be LOCAL or a valid CRS such as EPSG:3035.")
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

# Read global max height
max_height_global <- as.numeric(readLines(max_height_file, warn = FALSE))
if (is.na(max_height_global)) stop("Invalid max height value in the file.")

# Output raster folder
raster_output_path <- file.path(output_folder)
dir.create(raster_output_path, showWarnings = FALSE, recursive = TRUE)

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

las_crs_obj <- tryCatch(lidR::st_crs(las_combined), error = function(e) NULL)
if (local_mode) {
  lidR::st_crs(las_combined) <- sf::NA_crs_
  las_crs_obj <- sf::NA_crs_
  cat("Coordinate mode: LOCAL; output rasters will have no CRS.\n")
} else if (crs_missing(las_crs_obj)) {
  lidR::st_crs(las_combined) <- configured_crs
  las_crs_obj <- lidR::st_crs(las_combined)
  cat("LAS had no CRS; using configured fallback: ", coordinates, "\n", sep = "")
}

# Filter to non-ground points above 1 meter
las_nonground <- filter_poi(las_combined, Classification != 2 )

# Compute SHANNON indices rasters

compute_H <- function (z, by = 1, zmax = NULL) 
{
  if (is.null(zmax)) 
    zmax = max(z)
  if (zmax < 2 * by) 
    return(NA_real_)
  if (min(z) < 0) 
    return(NA_real_)
  bk <- seq(0, ceiling(zmax/by) * by, by)
  hist <- findInterval(z, bk)
  hist <- tabulate(hist, length(bk) - 1)
  hist <- hist/sum(hist)
  p <- hist[hist > 0]
  S <- -sum(p * log(p))
  return(S)
}

# compute_H10 <- function (z, by = 1, zmax = NULL, L = 10) 
# {
#   z <- z[is.finite(z)]
#   if (length(z) == 0) return(NA_real_)
#   
#   if (is.null(zmax)) 
#     zmax <- max(z)
#   if (zmax < 2 * by) 
#     return(NA_real_)
#   if (min(z) < 0) 
#     return(NA_real_)
#   
#   if (zmax == 0) {
#     return(0)
#   } else {
#     dz  <- zmax / L
#     bk  <- seq(0, zmax + dz, by = dz)
#     if (length(bk) < (L + 1L)) {
#       bk <- seq(0, zmax, length.out = L + 1L)
#     }
#     hist <- findInterval(z, bk)
#     hist <- tabulate(hist, length(bk) - 1)
#     hist <- hist/sum(hist)
#     p <- hist[hist > 0]
#     S <- -sum(p * log2(p))
#     return(S)
#   }
# }

r_H <- grid_metrics(las_nonground, ~ compute_H(Z, by = 1, zmax=max_height_global), res = 1)
# r_H10 <- grid_metrics(las_nonground, ~ compute_H10(Z, by = 1, zmax=max_height_global, L=10), res = 1) 

#Compute ENL rasters

compute_ENL1 <- function (z, by = 1, zmax = NULL) 
{
  if (is.null(zmax)) 
    zmax = max(z)
  if (zmax < 2 * by) 
    return(NA_real_)
  if (min(z) < 0) 
    return(NA_real_)
  bk <- seq(0, ceiling(zmax/by) * by, by)
  hist <- findInterval(z, bk)
  hist <- tabulate(hist, length(bk) - 1)
  hist <- hist/sum(hist)
  p <- hist[hist > 0]
  S <- -sum(p * log(p))
  ENL1<- exp(S)
  return(ENL1)
}

compute_ENL2 <- function (z, by = 1, zmax = NULL) 
{
  if (is.null(zmax)) 
    zmax = max(z)
  if (zmax < 2 * by) 
    return(NA_real_)
  if (min(z) < 0) 
    return(NA_real_)
  bk <- seq(0, ceiling(zmax/by) * by, by)
  hist <- findInterval(z, bk)
  hist <- tabulate(hist, length(bk) - 1)
  hist <- hist/sum(hist)
  p <- hist[hist > 0]
  ENL2<- 1 / sum((p)^2) 
  return(ENL2)
}

r_ENL1 <- grid_metrics(las_nonground, ~ compute_ENL1(Z, by = 1, zmax=max_height_global), res = 1)
r_ENL2 <- grid_metrics(las_nonground, ~ compute_ENL2(Z, by = 1, zmax=max_height_global), res = 1)

# Compute VCI 
r_VCI <- grid_metrics(las_nonground, ~ VCI(Z, max_height_global, by = 1), res = 1)

# Compute UCI
understory_top<- 3.0 #m
r_UCI <- grid_metrics(las_nonground, ~ VCI(Z, understory_top, by = 0.4), res = 1)

# Compute NRD
compute_NRD <-function(z,
         i = 0.5, j = 3.0) {
  z <- z[is.finite(z)]
  if (length(z) == 0) return(NA_real_)
  num <- sum(z >= i & z < j)
  den <- sum(z >= 0 & z < j)
  NRD <- if (den > 0) num / den else NA_real_
  return(NRD)
  }

r_NRD <- grid_metrics(las_nonground, ~ compute_NRD(Z,i = 0.5, j = 3.0), res = 1)

# Set CRS and save rasters. LOCAL outputs intentionally remain unreferenced;
# georeferenced outputs retain the actual LAS CRS, with the configured CRS used
# only as a fallback when the LAS metadata are absent.
r_H    <- set_output_crs(r_H,    las_crs_obj)
#r_H10 <- set_output_crs(r_H10,  las_crs_obj)
r_ENL1 <- set_output_crs(r_ENL1, las_crs_obj)
r_ENL2 <- set_output_crs(r_ENL2, las_crs_obj)
r_VCI  <- set_output_crs(r_VCI,  las_crs_obj)
r_UCI  <- set_output_crs(r_UCI,  las_crs_obj)
r_NRD  <- set_output_crs(r_NRD,  las_crs_obj)


writeRaster(r_H,      file.path(raster_output_path, paste0(base_name, "_H.tif")),     overwrite = TRUE)
#writeRaster(r_H10,   file.path(raster_output_path, paste0(base_name, "_H10.tif")),   overwrite = TRUE)
writeRaster(r_ENL1,        file.path(raster_output_path, paste0(base_name, "_ENL1.tif")),       overwrite = TRUE)
writeRaster(r_ENL2,      file.path(raster_output_path, paste0(base_name, "_ENL2.tif")),     overwrite = TRUE)
writeRaster(r_VCI,   file.path(raster_output_path, paste0(base_name, "_VCI.tif")),   overwrite = TRUE)
writeRaster(r_UCI,        file.path(raster_output_path, paste0(base_name, "_UCI.tif")),       overwrite = TRUE)
writeRaster(r_NRD,      file.path(raster_output_path, paste0(base_name, "_NRD.tif")),     overwrite = TRUE)

# Summary metrics 

metrics <- data.frame(
  H   = round(cloud_metrics(las_nonground, ~ compute_H(Z, by = 1, zmax = max_height_global))[[1]], 2),
  #H10 = round(cloud_metrics(las_nonground, ~ compute_H10(Z, by = 1, zmax = max_height_global, L = 10))[[1]], 2),
  ENL1 = round(cloud_metrics(las_nonground, ~ compute_ENL1(Z, by = 1, zmax = max_height_global))[[1]], 2),
  ENL2 = round(cloud_metrics(las_nonground, ~ compute_ENL2(Z, by = 1, zmax = max_height_global))[[1]], 2),
  VCI  = round(cloud_metrics(las_nonground, ~ VCI(Z, max_height_global, by = 1))[[1]], 2),
  UCI  = round(cloud_metrics(las_nonground, ~ VCI(Z, understory_top, by = 0.4))[[1]], 2),
  NRD  = round(cloud_metrics(las_nonground, ~ compute_NRD(Z, i = 0.5, j = 3.0))[[1]], 2)
)

write.table(metrics, file = file.path(output_folder, paste0(base_name, "_vertical_structure_metrics.txt")),
            row.names = FALSE, col.names = TRUE, sep = "        ", quote = FALSE)


# Save PNGs
png(file.path(output_folder, paste0(base_name, "_H.png")), width = 800, height = 600)
plot(r_H, main = "Shannon's entropy (1m-layers)", xlab = "X", ylab = "Y")
dev.off()

# png(file.path(output_folder, paste0(base_name, "_H10.png")), width = 800, height = 600)
# plot(r_H10, main = "Shannon's entropy (10 layers)", xlab = "X", ylab = "Y")
# dev.off()

png(file.path(output_folder, paste0(base_name, "_ENL1.png")), width = 800, height = 600)
plot(r_ENL1, main = "Effective Number of Layers (based on Shannon)", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_ENL2.png")), width = 800, height = 600)
plot(r_ENL2, main = "Effective Number of Layers (based on Simpson)", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_VCI.png")), width = 800, height = 600)
plot(r_VCI, main = "Vertical Complexity Index (VCI)", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_UCI.png")), width = 800, height = 600)
plot(r_UCI, main = "Understory Complexity Index (UCI)", xlab = "X", ylab = "Y")
dev.off()

png(file.path(output_folder, paste0(base_name, "_NRD.png")), width = 800, height = 600)
plot(r_NRD, main = "Normalized Relative Density (NRD)", xlab = "X", ylab = "Y")
dev.off()


cat("Vertical structure metrics processing complete. Results saved to:", output_folder, "\n")
