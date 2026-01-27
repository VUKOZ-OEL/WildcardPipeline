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
if [ $# -lt 5 ]; then
  echo "Usage: $0 <full/path/to/file.laz> <full/path/to/polygon.shp> <nOfElements> <voxelSize> <site_id> <coordinates>" >> "$LOG_FILE"
  timelog "exiting: missing parameters"
  exit 1
fi

inputPath="$1"
inputPolygonPath="$2"
nOfElements="$3"
voxelSize="$4"
coordinates="$6"

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

fullFileName=$(basename "$inputPath")            # file.laz
fileName="${fullFileName%.*}"                   # file
inputDir=$(dirname "$inputPath")                # original input dir

fullPolygonName=$(basename "$inputPolygonPath") # Velka_Ples.shp
polygonBase="${fullPolygonName%.*}"             # Velka_Ples
inputPolygonDir=$(dirname "$inputPolygonPath")  # original polygon dir

echo "Input LAZ: $inputPath" >> "$LOG_FILE"
echo "Polygon SHP: $inputPolygonPath" >> "$LOG_FILE"
echo "nOfElements: $nOfElements" >> "$LOG_FILE"
echo "voxelSize: $voxelSize" >> "$LOG_FILE"
echo "Coordinate system: $coordinates" >> "$LOG_FILE"
timelog "inputs logged"

# ---------------------------
# BASEDIR and scratch setup
# ---------------------------
BASEDIR="/storage/brno2/home/mihazivec/Wildcard"
cd "$BASEDIR"
timelog "changed directory to BASEDIR"

# Ensure SCRATCHDIR exists
if [ -z "${SCRATCHDIR-}" ] || [ ! -d "${SCRATCHDIR}" ]; then
  echo "ERROR: SCRATCHDIR is not set or does not exist. PBS scratch_local may not be active on this node." >> "$LOG_FILE"
  timelog "exiting: scratch not available"
  exit 2
fi

# Per-job workdir on scratch
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
SCR_TMP="$JOBSCRATCH/tmp"

mkdir -p "$SCR_IN" "$SCR_POLY" "$SCR_POLY_SEL" "$SCR_HOM" "$SCR_PROC" "$SCR_RES" "$SCR_TMP"
echo "Scratch workdir: $JOBSCRATCH" >> "$LOG_FILE"
timelog "scratch directories created"

# Force temp dirs to scratch (helps R and some tools)
export TMPDIR="$SCR_TMP"
export TMP="$SCR_TMP"
export TEMP="$SCR_TMP"
export R_TEMP_DIR="$SCR_TMP"
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
  cp -a "$SCR_RES/." "$BASEDIR/Results/$fileName/" >> "$LOG_FILE" 2>&1 || true
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
cp -v "$inputPath" "$SCR_IN/" >> "$LOG_FILE" 2>&1

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
# Select the correct polygon by site_id (global step)
# - The input shapefile may contain multiple polygons.
# - We select rows where site_id == $fileName, union them, and write a single-feature shapefile.
# ---------------------------
# SITE_ID can be provided as optional 5th argument (recommended) OR derived from the LAZ base name.
# Default derivation strips everything from SITE_ID_STRIP_MARKER onward (marker can be overridden via env var).
SITE_ID_STRIP_MARKER="${SITE_ID_STRIP_MARKER:-__ULS}"
SITE_ID="${5:-${fileName%%${SITE_ID_STRIP_MARKER}*}}"
echo "SITE_ID_STRIP_MARKER=$SITE_ID_STRIP_MARKER" >> "$LOG_FILE"
echo "Derived SITE_ID=$SITE_ID (from fileName=$fileName)" >> "$LOG_FILE"
SEL_SHP_NAME="selected_${SITE_ID}.shp"
SEL_SHP_HOST="$SCR_POLY_SEL/$SEL_SHP_NAME"

timelog "Selecting polygon for site_id=$SITE_ID"
# Optional debug: print all unique site_id values found in the polygon layer
timelog "Listing available site_id values from polygon file"
cat > "$SCR_TMP/list_site_ids.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
in_path <- args[[1]]

suppressPackageStartupMessages({
  library(sf)
})

p <- st_read(in_path, quiet = TRUE)

if (!("site_id" %in% names(p))) {
  stop(sprintf("Column 'site_id' not found in polygon file: %s", in_path))
}

ids <- sort(unique(as.character(p$site_id)))
cat(sprintf("Unique site_id count: %d\n", length(ids)))
for (id in ids) cat(id, "\n")
RSCRIPT

singularity exec \
  -B "$SCR_POLY:/polygon" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /tmp/list_site_ids.R "/polygon/$fullPolygonName" >> "$LOG_FILE" 2>&1

cat > "$SCR_TMP/select_polygon.R" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
in_path  <- args[[1]]
out_path <- args[[2]]
site_id  <- args[[3]]

suppressPackageStartupMessages({
  library(sf)
})

p <- st_read(in_path, quiet = TRUE)

if (!("site_id" %in% names(p))) {
  stop(sprintf("Column 'site_id' not found in polygon file: %s", in_path))
}

sel <- p[p$site_id == site_id, ]
if (nrow(sel) == 0) {
  stop(sprintf("No polygon found with site_id='%s' in %s", site_id, in_path))
}

# Union to ensure single output even if multiple features match site_id
geom <- st_union(sel)
out <- st_sf(site_id = site_id, geometry = geom, crs = st_crs(p))

# Overwrite if exists
st_write(out, out_path, delete_dsn = TRUE, quiet = TRUE)
RSCRIPT

singularity exec \
  -B "$SCR_POLY:/polygon" \
  -B "$SCR_POLY_SEL:/polygon_sel" \
  -B "$SCR_TMP:/tmp" \
  "$SIMG" \
  Rscript /tmp/select_polygon.R "/polygon/$fullPolygonName" "/polygon_sel/$SEL_SHP_NAME" "$SITE_ID" >> "$LOG_FILE" 2>&1

if [ ! -f "$SEL_SHP_HOST" ]; then
  echo "ERROR: Selected polygon shapefile not created: $SEL_SHP_HOST" >> "$LOG_FILE"
  timelog "exiting: polygon selection failed"
  exit 3
fi

timelog "Polygon selected: $SEL_SHP_HOST"


# Ensure homogenization script is executable
chmod +x "$BASEDIR/Code/CloudHomogenization.sh" >> "$LOG_FILE" 2>&1 || true

# ---------------------------
# Homogenization (on scratch)
# ---------------------------
timelog "Starting homogenization"
# Homogenized filename tag: voxelSize * 100 (e.g., 0.2 -> h20, 0.5 -> h50)
VOX_TAG=$(awk -v v="$voxelSize" 'BEGIN{ printf "%d", (v*100)+0.5 }')
HOM_OUT="$SCR_HOM/${fileName}_h${VOX_TAG}.laz"

"$BASEDIR/Code/CloudHomogenization.sh" \
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
  Rscript /scripts/Preprocessing.R /input /output /results /polygon_sel/"$SEL_SHP_NAME" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished preprocessing"

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
