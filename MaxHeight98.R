library(lidR)
library(raster)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript calculate_max_height.R <input_folder> <output_folder> <output_prefix>")
}

input_folder <- args[1]
output_folder <- args[2]
output_prefix <- args[3]

cat("Step: Calculating Max Height (98th percentile of Z > 2m non-ground points per 1m column)...\n")

# List all non-empty LAS/LAZ files
las_files <- list.files(input_folder, pattern = "\\.(las|laz)$", full.names = TRUE)
las_files <- las_files[file.info(las_files)$size > 0]

max_heights <- c()

# Process each LAS/LAZ file
for (file in las_files) {
  las_data <- tryCatch(readLAS(file), error = function(e) NULL)
  if (!is.null(las_data) && !is.empty(las_data)) {
    
    # Filter nonground points above 2m
    las_filtered <- filter_poi(las_data, Classification != 2 & Z > 2)
    
    if (!is.empty(las_filtered)) {
      max_height_raster <- tryCatch({
        grid_metrics(las_filtered, ~quantile(Z, probs = 0.98, na.rm = TRUE), res = 1)
      }, error = function(e) {
        cat("Error in file:", basename(file), "- skipping.\n")
        NULL
      })
      
      if (!is.null(max_height_raster)) {
        max_heights <- c(max_heights, max(values(max_height_raster), na.rm = TRUE))
      }
    }
  }
}

# Compute global max height as max of all column-wise 98th percentiles
max_height_global <- max(max_heights, na.rm = TRUE)
cat("Global Max Height (98th percentile of columns):", max_height_global, "\n")

# Write result to output file
output_file <- file.path(output_folder, paste0(output_prefix, "_maxHeight.txt"))
write(max_height_global, file = output_file, ncolumns = 1)
