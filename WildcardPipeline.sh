#!/bin/bash
set -Eeuo pipefail
IFS=$'\n\t'

# ---------------------------
# Logging / timing
# ---------------------------
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
export LOG_FILE="./log_${TIMESTAMP}.log"

START_TS=$(date +%s)
LAST_TS=$START_TS
fmt_hms() { local s=$1; printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60)); }
timelog() {
  local now ts elapsed delta
  now=$(date +%s)
  ts="$(date +"%Y-%m-%d %H:%M:%S")"
  elapsed=$(fmt_hms $((now - START_TS)))
  delta=$(fmt_hms $((now - LAST_TS)))
  echo "[$ts] (+$elapsed, Δ$delta) $*" >> "$LOG_FILE"
  LAST_TS=$now
}

trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]} failed at line ${LINENO} (exit code: ${rc})." >> "$LOG_FILE"; exit "${rc}"' ERR

echo "$(date) node ready" >> "$LOG_FILE"
echo "SCRATCHDIR=${SCRATCHDIR-<unset>}" >> "$LOG_FILE"
echo "TMPDIR=${TMPDIR-<unset>}" >> "$LOG_FILE"
timelog "timer initialized"

# ---------------------------
# Params
# ---------------------------
if [ $# -lt 6 ]; then
  echo "Usage: $0 <full/path/to/file.laz> <full/path/to/polygon.shp> <nOfElements> <voxelSize> <site_id> <LOCAL|EPSG:code> [max_normalized_height]" >> "$LOG_FILE"
  timelog "exiting: missing parameters"
  exit 1
fi

inputPath="$1"
inputPolygonPath="$2"
nOfElements="$3"
voxelSize="$4"
coordinates="$6"
maxNormalizedHeight="${7:-}"

# Accept exactly the two supported coordinate modes and pass a normalized value
# to every R stage. LOCAL is case-insensitive; EPSG codes may use either case.
coordinates="$(printf '%s' "$coordinates" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
coordinatesUpper="${coordinates^^}"
if [ "$coordinatesUpper" = "LOCAL" ]; then
  coordinates="LOCAL"
  coordinateMode="LOCAL metric Cartesian coordinates (no CRS)"
elif [[ "$coordinatesUpper" =~ ^EPSG:([0-9]+)$ ]]; then
  coordinates="EPSG:${BASH_REMATCH[1]}"
  coordinateMode="georeferenced ($coordinates)"
else
  echo "ERROR: coordinates must be LOCAL or EPSG:<code>; got: '$coordinates'" >> "$LOG_FILE"
  timelog "exiting: invalid coordinate mode"
  exit 2
fi

# Empty means disabled. Otherwise require a finite, non-negative decimal value.
if [ -n "$maxNormalizedHeight" ]; then
  if ! awk -v value="$maxNormalizedHeight" 'BEGIN {
    if (value !~ /^[+]?[0-9]+([.][0-9]+)?([eE][+-]?[0-9]+)?$/) exit 1
    numeric = value + 0
    if (numeric < 0) exit 1
  }'; then
    echo "ERROR: max_normalized_height must be empty or a finite non-negative number; got: $maxNormalizedHeight" >> "$LOG_FILE"
    timelog "exiting: invalid max normalized height"
    exit 2
  fi
fi

if [ ! -f "$inputPath" ]; then
  echo "ERROR: Input LAZ not found: $inputPath" >> "$LOG_FILE"
  timelog "exiting: missing input LAZ"
  exit 2
fi

if [ ! -f "$inputPolygonPath" ]; then
  echo "ERROR: Polygon SHP not found: $inputPolygonPath" >> "$LOG_FILE"
  timelog "exiting: missing polygon SHP"
  exit 2
fi

originalFullFileName=$(basename "$inputPath")    # original file.laz
originalFileName="${originalFullFileName%.*}"    # original file stem
inputExt="${originalFullFileName##*.}"           # laz / las
inputDir=$(dirname "$inputPath")                 # original input dir

# Scratch-safe filename: the original storage file is unchanged,
# but the scratch copy and all downstream processing use this name.
fileName="$(printf '%s' "$originalFileName" | sed 's/[^A-Za-z0-9._-]/_/g')"
fullFileName="${fileName}.${inputExt}"

if [ "$fullFileName" != "$originalFullFileName" ]; then
  echo "Original LAZ filename: $originalFullFileName" >> "$LOG_FILE"
  echo "Scratch LAZ filename: $fullFileName" >> "$LOG_FILE"
fi

fullPolygonName=$(basename "$inputPolygonPath") # Velka_Ples.shp
polygonBase="${fullPolygonName%.*}"             # Velka_Ples
inputPolygonDir=$(dirname "$inputPolygonPath")  # original polygon dir

echo "Input LAZ: $inputPath" >> "$LOG_FILE"
echo "Polygon SHP: $inputPolygonPath" >> "$LOG_FILE"
echo "nOfElements: $nOfElements" >> "$LOG_FILE"
echo "voxelSize: $voxelSize" >> "$LOG_FILE"
echo "Coordinate system: $coordinates" >> "$LOG_FILE"
echo "Coordinate mode: $coordinateMode" >> "$LOG_FILE"
if [ -n "$maxNormalizedHeight" ]; then
  echo "Max normalized height: $maxNormalizedHeight m above ground" >> "$LOG_FILE"
else
  echo "Max normalized height: disabled" >> "$LOG_FILE"
fi
timelog "inputs logged"

# ---------------------------
# BASEDIR and scratch setup
# ---------------------------
BASEDIR="/storage/brno2/home/mihazivec/WildcardWP3"
cd "$BASEDIR"
timelog "changed directory to BASEDIR"

# Ensure SCRATCHDIR exists
if [ -z "${SCRATCHDIR-}" ] || [ ! -d "${SCRATCHDIR}" ]; then
  echo "ERROR: SCRATCHDIR is not set or does not exist. PBS scratch_local may not be active on this node." >> "$LOG_FILE"
  timelog "exiting: scratch not available"
  exit 2
fi

# Per-job workdir on scratch
# fileName is already sanitized above, so JOBSCRATCH contains no spaces.
JOBSCRATCH="$SCRATCHDIR/wildcard_${PBS_JOBID:-nojobid}_$fileName"
if [ -d "$JOBSCRATCH" ]; then
  echo "Wiping existing scratch workdir: $JOBSCRATCH" >> "$LOG_FILE"
  rm -rf "$JOBSCRATCH"
fi

SCR_IN="$JOBSCRATCH/input"
SCR_POLY="$JOBSCRATCH/polygon"
SCR_POLY_SEL="$JOBSCRATCH/polygon_selected"
SCR_HOM="$JOBSCRATCH/Homogenized/$fileName"
SCR_PROC="$JOBSCRATCH/Processed/$fileName"
SCR_RES="$JOBSCRATCH/Results/$fileName"

# Critical: R temp directories must not contain spaces.
SCR_TMP="$SCRATCHDIR/wildcard_tmp_${PBS_JOBID:-nojobid}_$fileName"

mkdir -p "$SCR_IN" "$SCR_POLY" "$SCR_POLY_SEL" "$SCR_HOM" "$SCR_PROC" "$SCR_RES" "$SCR_TMP"
echo "Scratch workdir: $JOBSCRATCH" >> "$LOG_FILE"
timelog "scratch directories created"

# Force temp dirs to scratch (helps R and some tools)
export TMPDIR="$SCR_TMP"
export TMP="$SCR_TMP"
export TEMP="$SCR_TMP"
export R_TEMP_DIR="$SCR_TMP"
export R_SESSION_TMPDIR="$SCR_TMP"

# Force temp dirs inside Singularity to /tmp. Each singularity exec below
# binds $SCR_TMP to /tmp, so container-side /tmp still writes to scratch.
export SINGULARITYENV_TMPDIR="/tmp"
export SINGULARITYENV_TMP="/tmp"
export SINGULARITYENV_TEMP="/tmp"
export SINGULARITYENV_R_TEMP_DIR="/tmp"
export SINGULARITYENV_R_SESSION_TMPDIR="/tmp"

echo "TMPDIR=$TMPDIR" >> "$LOG_FILE"
timelog "temp dirs set to scratch"

# ---------------------------
# Copy-back / cleanup (runs on success and on failure)
# ---------------------------
PIPESTATUS_OK=0

copy_back() {
  echo "Copy-back started" >> "$LOG_FILE"
  mkdir -p "$BASEDIR/Processed/$fileName" "$BASEDIR/Results/$fileName"
  cp -a "$SCR_PROC/." "$BASEDIR/Processed/$fileName/" >> "$LOG_FILE" 2>&1 || true

  # Never publish preprocessing intermediates. Keep them on scratch after a
  # failed run for debugging, but exclude them from the persistent Results tree.
  rm -rf "$BASEDIR/Results/$fileName/.work" >> "$LOG_FILE" 2>&1 || true
  (
    shopt -s nullglob dotglob
    for result_item in "$SCR_RES"/*; do
      [ "${result_item##*/}" = ".work" ] && continue
      cp -a "$result_item" "$BASEDIR/Results/$fileName/"
    done
  ) >> "$LOG_FILE" 2>&1 || true

  # Also remove a stale .work directory from an older pipeline version.
  rm -rf "$BASEDIR/Results/$fileName/.work" >> "$LOG_FILE" 2>&1 || true
  echo "Copy-back finished" >> "$LOG_FILE"
}

cleanup() {
  # Preserve scratch if something failed (useful for debugging)
  if [ "${PIPESTATUS_OK:-0}" -eq 1 ]; then
    echo "Cleaning scratch: $JOBSCRATCH" >> "$LOG_FILE"
    rm -rf "$JOBSCRATCH" >> "$LOG_FILE" 2>&1 || true
  else
    echo "Not cleaning scratch (failure/debug): $JOBSCRATCH" >> "$LOG_FILE"
  fi
}

on_exit() {
  rc=$?
  if [ $rc -eq 0 ]; then PIPESTATUS_OK=1; else PIPESTATUS_OK=0; fi
  copy_back || true
  cleanup || true
  exit $rc
}

trap on_exit EXIT

# ---------------------------
# Stage inputs on scratch
# ---------------------------
timelog "copying inputs to scratch"
# Copy original input to scratch under the sanitized filename.
cp -v "$inputPath" "$SCR_IN/$fullFileName" >> "$LOG_FILE" 2>&1

# Copy the shapefile *set* (shp, shx, dbf, prj, etc.)
cp -v "$inputPolygonDir/$polygonBase".* "$SCR_POLY/" >> "$LOG_FILE" 2>&1

# Basic sanity check that the .shp exists on scratch
if [ ! -f "$SCR_POLY/$fullPolygonName" ]; then
  echo "ERROR: Polygon SHP not present on scratch after copy: $SCR_POLY/$fullPolygonName" >> "$LOG_FILE"
  timelog "exiting: polygon staging failed"
  exit 2
fi

timelog "inputs copied to scratch"

# ---------------------------
# Singularity / tools
# ---------------------------
if ! command -v singularity >/dev/null 2>&1; then
  echo "ERROR: singularity not found in PATH." >> "$LOG_FILE"
  timelog "exiting: singularity missing"
  exit 127
fi
echo "Singularity path: $(command -v singularity)" >> "$LOG_FILE"
timelog "singularity path logged"

SIMG="/storage/brno2/home/krucek/bluecat/singularity_img/lidr.img"
if [ ! -f "$SIMG" ]; then
  echo "ERROR: Singularity image not found: $SIMG" >> "$LOG_FILE"
  timelog "exiting: missing singularity image"
  exit 2
fi
# ---------------------------
# Prepare polygon input for downstream scripts
# - If the input shapefile already contains a single feature, use it directly.
# - If it contains multiple features, try to select by site_id and union to one feature.
# ---------------------------
# SITE_ID can be provided as optional 5th argument (recommended) OR derived from the LAZ base name.
# Default derivation strips everything from SITE_ID_STRIP_MARKER onward (marker can be overridden via env var).
SITE_ID_STRIP_MARKER="${SITE_ID_STRIP_MARKER:-__ULS}"
SITE_ID="${5:-${originalFileName%%${SITE_ID_STRIP_MARKER}*}}"
echo "SITE_ID_STRIP_MARKER=$SITE_ID_STRIP_MARKER" >> "$LOG_FILE"
echo "Derived SITE_ID=$SITE_ID (from originalFileName=$originalFileName)" >> "$LOG_FILE"

# Keep SITE_ID unchanged for attribute matching, but sanitize the output shapefile name.
safeSiteID="$(printf '%s' "$SITE_ID" | sed 's/[^A-Za-z0-9._-]/_/g')"
SEL_SHP_NAME="selected_${safeSiteID}.shp"
SEL_SHP_HOST="$SCR_POLY_SEL/$SEL_SHP_NAME"

timelog "Preparing polygon input (single-feature direct use, otherwise select by site_id=$SITE_ID)"
cat > "$SCR_TMP/prepare_polygon.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
in_path  <- args[[1]]
out_path <- args[[2]]
site_id  <- args[[3]]
coordinate_mode <- toupper(trimws(args[[4]]))

suppressPackageStartupMessages({
  library(sf)
})

p <- st_read(in_path, quiet = TRUE)

# LOCAL polygons are deliberately unreferenced. Clear any attached label before
# geometry repair so sf uses planar Cartesian operations and never reprojects.
if (identical(coordinate_mode, "LOCAL")) {
  suppressWarnings(st_crs(p) <- NA_crs_)
  cat("Coordinate mode: LOCAL; polygon CRS left undefined.\n")
}

p <- st_make_valid(p)
suppressWarnings({
  p <- st_zm(p, drop = TRUE, what = "ZM")
})

if (nrow(p) == 0L || all(st_is_empty(p))) {
  stop(sprintf("Polygon layer contains no usable geometry: %s", in_path))
}

feature_count <- nrow(p)
cat(sprintf("Polygon feature count: %d\n", feature_count))
cat("Polygon geometry flattened to 2D before selection/write.\n")

site_col_idx <- which(tolower(names(p)) == "site_id")
site_col <- if (length(site_col_idx) > 0) names(p)[site_col_idx[[1]]] else NA_character_

if (feature_count == 1) {
  cat("Single-feature polygon detected; using it directly without site_id filtering.\n")
  out <- st_sf(geometry = st_geometry(p))
} else {
  if (is.na(site_col)) {
    stop(sprintf(
      "Polygon layer contains %d features but no site_id column. Cannot disambiguate composite polygons in %s",
      feature_count, in_path
    ))
  }

  ids <- sort(unique(as.character(p[[site_col]])))
  cat(sprintf("Using attribute column: %s\n", site_col))
  cat(sprintf("Unique site_id count: %d\n", length(ids)))
  for (id in ids) cat(id, "\n")

  sel <- p[as.character(p[[site_col]]) == site_id, ]
  if (nrow(sel) == 0) {
    stop(sprintf(
      "No polygon found with site_id='%s' in %s (column: %s)",
      site_id, in_path, site_col
    ))
  }

  cat(sprintf("Matched %d feature(s) for site_id=%s; unioning to single output feature.\n", nrow(sel), site_id))
  geom <- st_union(sel)
  out <- st_sf(site_id = site_id, geometry = geom, crs = st_crs(p))
}

out <- st_make_valid(out)
suppressWarnings({
  out <- st_zm(out, drop = TRUE, what = "ZM")
})

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
old_sidecars <- Sys.glob(paste0(tools::file_path_sans_ext(out_path), ".*"))
if (length(old_sidecars)) unlink(old_sidecars, recursive = TRUE, force = TRUE)

st_write(out, out_path, driver = "ESRI Shapefile", quiet = TRUE)
cat(sprintf("Prepared polygon written to %s\n", out_path))
RSCRIPT

singularity exec \
  -B "$SCR_POLY:/polygon" \
  -B "$SCR_POLY_SEL:/polygon_sel" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /tmp/prepare_polygon.R "/polygon/$fullPolygonName" "/polygon_sel/$SEL_SHP_NAME" "$SITE_ID" "$coordinates" >> "$LOG_FILE" 2>&1

if [ ! -f "$SEL_SHP_HOST" ]; then
  echo "ERROR: Prepared polygon shapefile not created: $SEL_SHP_HOST" >> "$LOG_FILE"
  timelog "exiting: polygon preparation failed"
  exit 3
fi

timelog "Polygon prepared: $SEL_SHP_HOST"

# Ensure homogenization script is executable
chmod +x "$BASEDIR/Code/CloudHomogenizationGeometricNN.sh" >> "$LOG_FILE" 2>&1 || true

# ---------------------------
# Homogenization (on scratch)
# ---------------------------
timelog "Starting homogenization"
# Homogenized filename tag: voxelSize * 100 (e.g., 0.2 -> h20, 0.5 -> h50)
VOX_TAG=$(awk -v v="$voxelSize" 'BEGIN{ printf "%d", (v*100)+0.5 }')
HOM_OUT="$SCR_HOM/${fileName}_h${VOX_TAG}.laz"

"$BASEDIR/Code/CloudHomogenizationGeometricNN.sh" \
  "$SCR_IN/$fullFileName" \
  "$HOM_OUT" \
  "$voxelSize" >> "$LOG_FILE" 2>&1

if [ ! -f "$HOM_OUT" ]; then
  echo "ERROR: Homogenized output not found: $HOM_OUT" >> "$LOG_FILE"
  timelog "exiting: homogenization failed"
  exit 3
fi

file_count_hom=$(find "$SCR_HOM" -type f | wc -l)
echo "Homogenized $fileName contains $file_count_hom files (scratch)" >> "$LOG_FILE"
timelog "finished homogenization"

# ---------------------------
# Pipeline steps (all on scratch)
# Preprocessing reads from homogenized directory
# ---------------------------
export RESUME=0

timelog "Starting preprocessing"
singularity exec \
  -B "$SCR_HOM:/input" \
  -B "$SCR_PROC:/output" \
  -B "$SCR_RES:/results" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$SCR_POLY:/polygon" \
  -B "$SCR_POLY_SEL:/polygon_sel" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/Preprocessing.R \
    /input \
    /output \
    /results \
    /polygon_sel/"$SEL_SHP_NAME" \
    "$coordinates" \
    "$maxNormalizedHeight" >> "$LOG_FILE" 2>&1
timelog "finished preprocessing"

# Preprocessing normally removes this through R on.exit(). Remove it again at
# the orchestration layer so successful runs cannot carry intermediates forward.
if [ -d "$SCR_RES/.work" ]; then
  rm -rf "$SCR_RES/.work"
  timelog "removed preprocessing .work directory"
fi

file_count_proc=$(find "$SCR_PROC" -type f | wc -l)
echo "Cloud $fileName contains $file_count_proc files (scratch processed)" >> "$LOG_FILE"
timelog "post-preprocessing file count logged"

timelog "Starting MaxHeight98"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_PROC:/output" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/MaxHeight98.R /input /output "$fileName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished MaxHeight98"

# Common maxHeight file bind target (created by MaxHeight98.R)
MAXHEIGHT_HOST_PATH="$SCR_PROC/${fileName}_maxHeight.txt"
if [ ! -f "$MAXHEIGHT_HOST_PATH" ]; then
  echo "ERROR: MaxHeight file not found at $MAXHEIGHT_HOST_PATH (MaxHeight98 step failed or wrote elsewhere)." >> "$LOG_FILE"
  timelog "exiting: maxHeight missing"
  exit 3
fi

timelog "Starting HeightVariabilityMetrics"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_RES:/output" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$MAXHEIGHT_HOST_PATH:/maxHeight" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/HeightVariabilityMetrics.R /input /output "$fileName" /maxHeight "$nOfElements" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished HeightVariabilityMetrics"

timelog "Starting GappinesMetrics"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_RES:/output" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$MAXHEIGHT_HOST_PATH:/maxHeight" \
  -B "$SCR_POLY:/polygon" \
  -B "$SCR_POLY_SEL:/polygon_sel" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/GappinesMetrics.R /input /output /maxHeight "$fileName" /polygon_sel/"$SEL_SHP_NAME" "$nOfElements" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished GappinesMetrics"

timelog "Starting CanopyDensity"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_RES:/output" \
  -B "$MAXHEIGHT_HOST_PATH:/maxHeight" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/CanopyDensity.R /input /output /maxHeight "$fileName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished CanopyDensity"

timelog "Starting CanopyVolumeMetrics"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_RES:/output" \
  -B "$MAXHEIGHT_HOST_PATH:/maxHeight" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/CanopyVolumeMetrics.R /input /output /maxHeight "$fileName" "$nOfElements" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished CanopyVolumeMetrics"

timelog "Starting VerticalStructureMetrics"
singularity exec \
  -B "$SCR_PROC:/input" \
  -B "$SCR_RES:/output" \
  -B "$BASEDIR/Code:/scripts" \
  -B "$MAXHEIGHT_HOST_PATH:/maxHeight" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /scripts/VerticalStructureMetrics.R /input /output /maxHeight "$fileName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished VerticalStructureMetrics"

timelog "pipeline completed (results will be copied back by trap)"
