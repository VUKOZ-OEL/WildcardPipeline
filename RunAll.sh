#!/bin/bash
#PBS -N Wildcard_job
#PBS -l select=1:ncpus=5:mem=256gb:scratch_local=450gb
#PBS -l walltime=48:00:00
#PBS -v AWS_SHARED_CREDENTIALS_FILE=/storage/brno2/home/mihazivec/.aws/credentials

set -Eeuo pipefail
IFS=$'\n\t'

trap 'rc=$?; echo "ERROR: ${BASH_SOURCE[0]} failed at line ${LINENO} (exit code: ${rc})." >&2; exit "${rc}"' ERR

# ------------------------------------------
# MODIFY THESE VARIABLES BEFORE EACH RUN
S3_SUBFOLDER="UNIPA"                 # Folder/prefix in wp3 bucket containing the LAZ files

NOFELEMENTS="10"
VOXELSIZE="0.2"                      # meters

# Coordinate mode for every site in this batch:
#   LOCAL       = aligned local Cartesian X/Y coordinates in metres; no CRS
#   EPSG:<code> = georeferenced data (for example EPSG:3035)
# Change this setting between batches when your source collection contains
# both local-coordinate and georeferenced groups.
COORDINATES="EPSG:3035"

# Optional hard ceiling in metres above ground after height normalization.
# Points with normalized Z greater than this value are removed from the
# scratch working cloud before CHM generation and downstream metrics.
# Leave empty to disable the hard height cutoff.
MAX_NORMALIZED_HEIGHT=""

# ---- SITE_ID DERIVATION (EDIT IF NEEDED) ----
# This preserves the same SITE_ID derivation used by the original script.
# With ".laz" as the marker, SITE_ID is simply the LAZ filename without .laz.
SITE_ID_STRIP_MARKER=".laz"
echo "${SITE_ID_STRIP_MARKER}"
# --------------------------------------------

# Normalize and validate the coordinate setting before downloading any data.
COORDINATES="$(printf '%s' "$COORDINATES" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
COORDINATES_UPPER="${COORDINATES^^}"
if [ "$COORDINATES_UPPER" = "LOCAL" ]; then
  COORDINATES="LOCAL"
  COORDINATE_MODE="LOCAL metric Cartesian coordinates (no CRS)"
elif [[ "$COORDINATES_UPPER" =~ ^EPSG:([0-9]+)$ ]]; then
  COORDINATES="EPSG:${BASH_REMATCH[1]}"
  COORDINATE_MODE="georeferenced ($COORDINATES)"
else
  echo "ERROR: COORDINATES must be LOCAL or EPSG:<code>; got: '$COORDINATES'" >&2
  exit 2
fi

# Base directory
BASEDIR="/storage/brno2/home/mihazivec/WildcardWP3"

# s5cmd configuration
S5_PROFILE="s3cmd_access"
S5_ENDPOINT="https://s3.cl4.du.cesnet.cz"
S3_BUCKET="wp3"

s5() { s5cmd --profile "$S5_PROFILE" --endpoint-url="$S5_ENDPOINT" "$@"; }

# ---- Preconditions ----
command -v s5cmd >/dev/null 2>&1 || { echo "ERROR: s5cmd not found in PATH." >&2; exit 127; }
[ -d "$BASEDIR" ] || { echo "ERROR: BASEDIR does not exist: $BASEDIR" >&2; exit 2; }

PIPELINE="$BASEDIR/WildcardPipeline.sh"
[ -f "$PIPELINE" ] || { echo "ERROR: Pipeline script not found: $PIPELINE" >&2; exit 2; }
chmod +x "$PIPELINE"

echo "SCRATCHDIR=${SCRATCHDIR-<unset>}"
echo "TMPDIR=${TMPDIR-<unset>}"

cd "$BASEDIR"

FILEBASE="${S3_SUBFOLDER%/}"
LOCAL_TARGET_DIR="$BASEDIR/Data/$FILEBASE"
LOCAL_POLYGON_DIR="$LOCAL_TARGET_DIR/polygons"
POLY_PREFIX="s3://$S3_BUCKET/${S3_SUBFOLDER%/}/polygons"

mkdir -p "$LOCAL_TARGET_DIR" "$LOCAL_POLYGON_DIR"

if [ -n "$MAX_NORMALIZED_HEIGHT" ]; then
  HEIGHT_CUTOFF_LABEL="$MAX_NORMALIZED_HEIGHT m above ground"
else
  HEIGHT_CUTOFF_LABEL="disabled"
fi

echo "Batch input:"
echo "  LAZ folder : s3://$S3_BUCKET/${S3_SUBFOLDER%/}/"
echo "  polygons   : $POLY_PREFIX/"
echo "  Params     : NOFELEMENTS=$NOFELEMENTS, VOXELSIZE=$VOXELSIZE"
echo "  Coordinates: $COORDINATE_MODE"
echo "  Max normalized height: $HEIGHT_CUTOFF_LABEL"
echo

# ---- Find every LAZ file in the selected S3 folder ----
# s5cmd ls output ends with the S3 object path, so awk extracts the last field.
echo "Finding LAZ files: s3://$S3_BUCKET/${S3_SUBFOLDER%/}/*.laz"
LAZ_LISTING="$(s5 ls "s3://$S3_BUCKET/${S3_SUBFOLDER%/}/*.laz" || true)"

mapfile -t S3_LAZ_PATHS < <(
  printf '%s\n' "$LAZ_LISTING" | awk 'NF {print $NF}' | sort
)

if [ "${#S3_LAZ_PATHS[@]}" -eq 0 ]; then
  echo "ERROR: No .laz files found under s3://$S3_BUCKET/${S3_SUBFOLDER%/}/" >&2
  exit 2
fi

echo "Found ${#S3_LAZ_PATHS[@]} LAZ file(s)."
echo

PROCESSED=0
FAILED=0
SKIPPED=0
FAILED_FILES=()
SKIPPED_FILES=()

# ---- Process all LAZ files sequentially ----
for S3_LAZ_PATH in "${S3_LAZ_PATHS[@]}"; do
  S3_FILENAME="$(basename "$S3_LAZ_PATH")"
  FILE_STEM="${S3_FILENAME%.laz}"
  S3_PATH_REL="${S3_SUBFOLDER%/}/$S3_FILENAME"

  # Matching polygon has exactly the same base name as the LAZ file.
  POLYGON_SHP_FILENAME="${FILE_STEM}.shp"
  REMOTE_POLYGON_BASE="$POLY_PREFIX/$FILE_STEM"

  # Preserve the original SITE_ID logic.
  SITE_ID="${S3_FILENAME%%${SITE_ID_STRIP_MARKER}*}"

  LOCAL_TARGET_PATH="$LOCAL_TARGET_DIR/$S3_FILENAME"
  POLYGON_PATH="$LOCAL_POLYGON_DIR/$POLYGON_SHP_FILENAME"

  echo "============================================================"
  echo "Processing: $S3_FILENAME"
  echo "  LAZ      : s3://$S3_BUCKET/$S3_PATH_REL"
  echo "  polygon  : ${REMOTE_POLYGON_BASE}.shp"
  echo "  SITE_ID  : $SITE_ID"
  echo "  CRS mode : $COORDINATES"
  echo "============================================================"

  # ---- Verify matching polygon exists ----
  if ! s5 ls "${REMOTE_POLYGON_BASE}.shp" >/dev/null 2>&1; then
    echo "WARNING: Matching polygon not found; skipping this LAZ:" >&2
    echo "  ${REMOTE_POLYGON_BASE}.shp" >&2
    SKIPPED=$((SKIPPED + 1))
    SKIPPED_FILES+=("$S3_FILENAME")
    echo
    continue
  fi

  # ---- Download LAZ ----
  echo "Downloading LAZ -> $LOCAL_TARGET_PATH"
  if ! s5 cp "s3://$S3_BUCKET/$S3_PATH_REL" "$LOCAL_TARGET_PATH"; then
    echo "ERROR: Failed to download LAZ: $S3_FILENAME" >&2
    FAILED=$((FAILED + 1))
    FAILED_FILES+=("$S3_FILENAME (LAZ download failed)")
    echo
    continue
  fi

  # ---- Download matching polygon and all available shapefile sidecars ----
  # This fetches e.g. .shp, .shx, .dbf, .prj, .cpg with the same base name.
  # Remove cached sidecars for this exact polygon first. This is important when
  # a local-coordinate shapefile has no .prj: a stale .prj from an older copy
  # must not be mistaken for current CRS metadata.
  shopt -s nullglob
  OLD_POLYGON_SIDECARS=("$LOCAL_POLYGON_DIR/$FILE_STEM".*)
  if [ "${#OLD_POLYGON_SIDECARS[@]}" -gt 0 ]; then
    rm -f -- "${OLD_POLYGON_SIDECARS[@]}"
  fi
  shopt -u nullglob
  echo "Downloading polygon files -> $LOCAL_POLYGON_DIR"
  if ! s5 cp "${REMOTE_POLYGON_BASE}.*" "$LOCAL_POLYGON_DIR/"; then
    echo "ERROR: Failed to download polygon files for: $S3_FILENAME" >&2
    FAILED=$((FAILED + 1))
    FAILED_FILES+=("$S3_FILENAME (polygon download failed)")
    echo
    continue
  fi

  if [ ! -f "$POLYGON_PATH" ]; then
    echo "ERROR: Expected local polygon was not produced: $POLYGON_PATH" >&2
    FAILED=$((FAILED + 1))
    FAILED_FILES+=("$S3_FILENAME (local .shp missing)")
    echo
    continue
  fi

  echo "Using polygon shapefile: $POLYGON_PATH"
  echo

  # ---- Run pipeline ----
  if "$PIPELINE" \
    "$LOCAL_TARGET_PATH" \
    "$POLYGON_PATH" \
    "$NOFELEMENTS" \
    "$VOXELSIZE" \
    "$SITE_ID" \
    "$COORDINATES" \
    "$MAX_NORMALIZED_HEIGHT"; then

    PROCESSED=$((PROCESSED + 1))
    echo
    echo "SUCCESS: $S3_FILENAME"
  else
    rc=$?
    FAILED=$((FAILED + 1))
    FAILED_FILES+=("$S3_FILENAME (pipeline exit code $rc)")
    echo
    echo "ERROR: Pipeline failed for $S3_FILENAME with exit code $rc; continuing with next file." >&2
  fi

  echo
done

# ---- Batch summary ----
echo "============================================================"
echo "Batch complete"
echo "  Successful : $PROCESSED"
echo "  Failed     : $FAILED"
echo "  Skipped    : $SKIPPED"
echo "  Total      : ${#S3_LAZ_PATHS[@]}"
echo "============================================================"

if [ "${#FAILED_FILES[@]}" -gt 0 ]; then
  echo "Failed files:"
  printf '  - %s\n' "${FAILED_FILES[@]}"
fi

if [ "${#SKIPPED_FILES[@]}" -gt 0 ]; then
  echo "Skipped files (matching .shp not found):"
  printf '  - %s\n' "${SKIPPED_FILES[@]}"
fi

# Return non-zero if anything actually failed. Missing polygons are reported as skipped.
if [ "$FAILED" -gt 0 ]; then
  exit 1
fi
