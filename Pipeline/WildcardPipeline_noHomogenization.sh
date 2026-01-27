#!/bin/bash
set -euo pipefail

# ---------------------------
# Logging / timing
# ---------------------------
TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
export LOG_FILE="./log_$TIMESTAMP.log"

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

timelog "timer initialized"
echo "$(date) node ready" >> "$LOG_FILE"
timelog "node ready"

# ---------------------------
# Params
# ---------------------------
if [ $# -lt 5 ]; then
  echo "Usage: $0 <full/path/to/file.laz> <full/path/to/polygon.shp> <nOfElements> <voxelSize> <coordinates>" >> "$LOG_FILE"
  timelog "exiting: missing parameters"
  exit 1
fi

inputPath="$1"
inputPolygonPath="$2"
nOfElements="$3"
voxelSize="$4"
coordinates="$5"

fullFileName=$(basename "$inputPath")          # file.laz
fileName="${fullFileName%.*}"                 # file
inputDir=$(dirname "$inputPath")              # original input dir

fullPolygonName=$(basename "$inputPolygonPath")  # Velka_Ples.shp
polygonBase="${fullPolygonName%.*}"              # Velka_Ples
inputPolygonDir=$(dirname "$inputPolygonPath")   # original polygon dir

echo "Input LAZ: $inputPath" >> "$LOG_FILE"
echo "Polygon SHP: $inputPolygonPath" >> "$LOG_FILE"
echo "nOfElements: $nOfElements" >> "$LOG_FILE"
echo "voxelSize: $voxelSize" >> "$LOG_FILE"
echo "coordinates: $coordinates" >> "$LOG_FILE"
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
SCR_PROC="$JOBSCRATCH/Processed/$fileName"
SCR_RES="$JOBSCRATCH/Results/$fileName"
SCR_TMP="$JOBSCRATCH/tmp"

mkdir -p "$SCR_IN" "$SCR_POLY" "$SCR_PROC" "$SCR_RES" "$SCR_TMP"
echo "Scratch workdir: $JOBSCRATCH" >> "$LOG_FILE"
timelog "scratch directories created"

# Force temp dirs to scratch (helps R and some tools)
export TMPDIR="$SCR_TMP"
export TMP="$SCR_TMP"
export TEMP="$SCR_TMP"
export R_TEMP_DIR="$SCR_TMP"

echo "TMPDIR=$TMPDIR" >> "$LOG_FILE"
timelog "temp dirs set to scratch"

# Copy inputs to scratch
cp -v "$inputPath" "$SCR_IN/" >> "$LOG_FILE" 2>&1
# Copy the shapefile *set* (shp, shx, dbf, prj, etc.)
cp -v "$inputPolygonDir/$polygonBase".* "$SCR_POLY/" >> "$LOG_FILE" 2>&1
timelog "inputs copied to scratch"

# For container paths
SCR_LAZ="/input/$fullFileName"
SCR_SHP="/polygon/$fullPolygonName"

# ---------------------------
# Singularity
# ---------------------------
echo "Singularity path: $(which singularity)" >> "$LOG_FILE"
timelog "singularity path logged"

SIMG="/storage/brno2/home/krucek/bluecat/singularity_img/lidr.img"

# Copy-back function (runs on success and on failure)
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
# Pipeline steps (all on scratch)
# ---------------------------
export RESUME=0
echo "Listing scratch outputs before preprocessing:" >> "$LOG_FILE"
find "$SCR_PROC" "$SCR_RES" -type f -maxdepth 3 -print >> "$LOG_FILE" 2>&1 || true
timelog "Starting preprocessing"
singularity exec \
   -B "$SCR_IN:/input" \
   -B "$SCR_PROC:/output" \
   -B "$SCR_RES:/results" \
   -B "$BASEDIR/Code:/scripts" \
   -B "$SCR_POLY:/polygon" \
   -B "$SCR_TMP:/tmp" \
   "$SIMG" \
   Rscript /scripts/Preprocessing.R /input /output /results /polygon/"$fullPolygonName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished preprocessing"

file_count2=$(find "$SCR_PROC" -type f | wc -l)
echo "Cloud $fileName contains $file_count2 files (scratch)" >> "$LOG_FILE"
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

# Require the file to exist and be non-empty (-s). Use -f if you only care about existence.
if [ ! -s "$MAXHEIGHT_HOST_PATH" ]; then
  echo "ERROR: MaxHeight file missing or empty at $MAXHEIGHT_HOST_PATH (MaxHeight98 step failed or wrote elsewhere)." >> "$LOG_FILE"
  echo "Debug: listing scratch processed outputs (first 100 files):" >> "$LOG_FILE"
  find "$SCR_PROC" -maxdepth 3 -type f | head -n 100 >> "$LOG_FILE" 2>&1 || true
  timelog "exiting: maxHeight missing"
  exit 3
fi

timelog "Starting HeightVariabilityMetrics"
singularity exec \
    -B "$SCR_PROC:/input" \
    -B "$SCR_RES:/output" \
    -B "$BASEDIR/Code:/scripts" \
    -B "$SCR_PROC/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$SCR_TMP:/tmp" \
    "$SIMG" \
    Rscript /scripts/HeightVariabilityMetrics.R /input /output "$fileName" /maxHeight "$nOfElements" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished HeightVariabilityMetrics"

timelog "Starting GappinesMetrics"
singularity exec \
    -B "$SCR_PROC:/input" \
    -B "$SCR_RES:/output" \
    -B "$SCR_PROC/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$BASEDIR/Code:/scripts" \
    -B "$SCR_POLY:/polygon" \
    -B "$SCR_TMP:/tmp" \
    "$SIMG" \
    Rscript /scripts/GappinesMetrics.R /input /output /maxHeight "$fileName" /polygon/"$fullPolygonName" "$nOfElements" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished GappinesMetrics"

timelog "Starting CanopyDensity"
singularity exec \
    -B "$SCR_PROC:/input" \
    -B "$SCR_RES:/output" \
    -B "$SCR_PROC/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$BASEDIR/Code:/scripts" \
    -B "$SCR_TMP:/tmp" \
    "$SIMG" \
    Rscript /scripts/CanopyDensity.R /input /output /maxHeight "$fileName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished CanopyDensity"

timelog "Starting CanopyVolumeMetrics"
singularity exec \
    -B "$SCR_PROC:/input" \
    -B "$SCR_RES:/output" \
    -B "$SCR_PROC/${fileName}_maxHeight.txt:/maxHeight" \
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
    -B "$SCR_PROC/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$SCR_TMP:/tmp" \
    "$SIMG" \
    Rscript /scripts/VerticalStructureMetrics.R /input /output /maxHeight "$fileName" "$coordinates" >> "$LOG_FILE" 2>&1
timelog "finished VerticalStructureMetrics"

timelog "pipeline completed (results will be copied back by trap)"
