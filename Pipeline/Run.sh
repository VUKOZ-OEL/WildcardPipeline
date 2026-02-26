#!/bin/bash
#PBS -N Wildcard_job
#PBS -l select=1:ncpus=5:mem=128gb:scratch_local=300gb
#PBS -l walltime=3:30:00
#PBS -v AWS_SHARED_CREDENTIALS_FILE=/storage/brno2/home/mihazivec/.aws/credentials

set -Eeuo pipefail
IFS=$'\n\t'

trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]} failed at line ${LINENO} (exit code: ${rc})." >&2; exit "${rc}"' ERR

# ------------------------------------------
# MODIFY THESE VARIABLES BEFORE EACH RUN
S3_SUBFOLDER="$1"                          # Folder/prefix in wp2 bucket ie "FVA-BW" | "VUK"
S3_FILENAME="$2"      # LAZ file to process ie: "FVA-BW__14__Wildseemoor__1123__ULS.laz"
POLYGON_SHP_FILENAME="$3"


COORDINATES=${3:-"EPSG:5514"}   
VOXELSIZE=${4:-0.2}   # is number was text
NOFELEMENTS=${5:-10}                   
                            

# If the polygons prefix contains MULTIPLE .shp files, set the exact .shp filename.


# ---- SITE_ID DERIVATION (EDIT IF NEEDED) ----
# Polygon layers contain multiple sites identified by the "site_id" column.
# We derive SITE_ID from the LAZ base name by stripping the suffix starting at this marker.
# Examples:
#   FVA-BW__10__Zweribach__NA__ULS_2   -> FVA-BW__10__Zweribach__NA
#   FVA-BW__14__Wildseemoor__1123__ULS -> FVA-BW__14__Wildseemoor__1123
SITE_ID_STRIP_MARKER="__ULS"
# --------------------------------------------
# ------------------------------------------

# Base directory
BASEDIR="/storage/brno2/home/mihazivec/Wildcard"

# s5cmd configuration
S5_PROFILE="s3cmd_access"
S5_ENDPOINT="https://s3.cl4.du.cesnet.cz"
S3_BUCKET="wp2"

s5() { s5cmd --profile "$S5_PROFILE" --endpoint-url="$S5_ENDPOINT" "$@"; }

# ---- Preconditions ----
command -v s5cmd >/dev/null 2>&1 || { echo "ERROR: s5cmd not found in PATH." >&2; exit 127; }
[ -d "$BASEDIR" ] || { echo "ERROR: BASEDIR does not exist: $BASEDIR" >&2; exit 2; }

echo "SCRATCHDIR=${SCRATCHDIR-<unset>}"
echo "TMPDIR=${TMPDIR-<unset>}"

cd "$BASEDIR"

S3_PATH_REL="${S3_SUBFOLDER%/}/$S3_FILENAME"
FILEBASE="${S3_FILENAME%.*}"

# Derive SITE_ID for polygon selection (passed to the pipeline)
SITE_ID="${FILEBASE%%${SITE_ID_STRIP_MARKER}*}"

LOCAL_TARGET_DIR="$BASEDIR/Data/$FILEBASE"
LOCAL_TARGET_PATH="$LOCAL_TARGET_DIR/$S3_FILENAME"

LOCAL_POLYGON_DIR="$LOCAL_TARGET_DIR/polygons"
mkdir -p "$LOCAL_TARGET_DIR" "$LOCAL_POLYGON_DIR"

echo "Selected input:"
echo "  LAZ      : s3://$S3_BUCKET/$S3_PATH_REL"
echo "  polygons : s3://$S3_BUCKET/${S3_SUBFOLDER%/}/polygons/"
echo "  Params   : NOFELEMENTS=$NOFELEMENTS, VOXELSIZE=$VOXELSIZE"
echo

# ---- Correct bucket/prefix listing (S3 needs wildcard) ----
echo "Listing: s3://$S3_BUCKET/${S3_SUBFOLDER%/}/*"
s5 ls "s3://$S3_BUCKET/${S3_SUBFOLDER%/}/*" || true
echo

# ---- Verify LAZ exists ----
if ! s5 ls "s3://$S3_BUCKET/$S3_PATH_REL" >/dev/null 2>&1; then
  echo "ERROR: LAZ not found in S3: s3://$S3_BUCKET/$S3_PATH_REL" >&2
  echo "Please verify S3_SUBFOLDER and S3_FILENAME (case-sensitive)." >&2
  exit 2
fi

# ---- Verify polygons exist (must use wildcard) ----
POLY_PREFIX="s3://$S3_BUCKET/${S3_SUBFOLDER%/}/polygons/"
if ! s5 ls "${POLY_PREFIX}*" >/dev/null 2>&1; then
  echo "ERROR: polygons prefix not found or empty in S3: ${POLY_PREFIX}" >&2
  echo "Tried listing: ${POLY_PREFIX}*" >&2
  exit 2
fi

echo "Remote polygon listing candidates:"
s5 ls "${POLY_PREFIX}*"     || true
s5 ls "${POLY_PREFIX}*/*"   || true
s5 ls "${POLY_PREFIX}*/*/*" || true
echo

# ---- Download LAZ ----
echo "Downloading LAZ -> $LOCAL_TARGET_PATH"
s5 cp "s3://$S3_BUCKET/$S3_PATH_REL" "$LOCAL_TARGET_PATH"

# ---- Download polygons ----
# NOTE: your s5cmd build and/or object layout may require deeper patterns than polygons/*.
# We try several depths and then hard-fail if no .shp appears locally.
echo "Downloading polygons -> $LOCAL_POLYGON_DIR"
for pat in "*" "*/*" "*/*/*" "*/*/*/*"; do
  s5 cp "${POLY_PREFIX}${pat}" "$LOCAL_POLYGON_DIR/" >/dev/null 2>&1 || true
done

if ! find "$LOCAL_POLYGON_DIR" -type f -name '*.shp' -print -quit | grep -q .; then
  echo "ERROR: polygons download produced no .shp files under: $LOCAL_POLYGON_DIR" >&2
  echo "This usually means polygons are nested deeper than 4 levels, or the files are not .shp." >&2
  echo "Check the 'Remote polygon listing candidates' above and extend the patterns if needed." >&2
  exit 2
fi

echo "Local polygon files (first 50):"
find "$LOCAL_POLYGON_DIR" -type f | head -n 50
echo

# ---- Resolve polygon .shp to use ----
POLYGON_PATH=""
if [ -n "$POLYGON_SHP_FILENAME" ]; then
  POLYGON_PATH_FOUND="$(find "$LOCAL_POLYGON_DIR" -maxdepth 10 -type f -name "$POLYGON_SHP_FILENAME" | head -n 1 || true)"
  if [ -z "$POLYGON_PATH_FOUND" ]; then
    echo "ERROR: Requested polygon shapefile not found: $POLYGON_SHP_FILENAME under $LOCAL_POLYGON_DIR" >&2
    exit 2
  fi
  POLYGON_PATH="$POLYGON_PATH_FOUND"
else
  mapfile -t SHP_FILES < <(find "$LOCAL_POLYGON_DIR" -maxdepth 10 -type f -name '*.shp' | sort)
  if [ "${#SHP_FILES[@]}" -eq 0 ]; then
    echo "ERROR: No .shp found under $LOCAL_POLYGON_DIR" >&2
    exit 2
  fi
  if [ "${#SHP_FILES[@]}" -gt 1 ]; then
    echo "ERROR: Multiple .shp files found under $LOCAL_POLYGON_DIR. Set POLYGON_SHP_FILENAME explicitly." >&2
    printf '  - %s\n' "${SHP_FILES[@]}" >&2
    exit 2
  fi
  POLYGON_PATH="${SHP_FILES[0]}"
fi

echo "Using polygon shapefile: $POLYGON_PATH"
echo

# ---- Run pipeline ----
PIPELINE="$BASEDIR/WildcardPipeline.sh"
[ -f "$PIPELINE" ] || { echo "ERROR: Pipeline script not found: $PIPELINE" >&2; exit 2; }

chmod +x "$PIPELINE"
"$PIPELINE" "$LOCAL_TARGET_PATH" "$POLYGON_PATH" "$NOFELEMENTS" "$VOXELSIZE" "$SITE_ID" "$COORDINATES"
