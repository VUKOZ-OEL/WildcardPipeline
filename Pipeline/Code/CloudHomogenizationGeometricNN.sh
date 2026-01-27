#!/bin/bash

# Strict mode
set -euo pipefail

#==============================
# 1) Args
#==============================
export SOURCE_FILE=${1:-}
export RESULT_FILE=${2:-}
export VOXEL_SPEC=${3:-0.005}  # "xy" or "xy,z" (e.g., "0.5" or "0.5,0.2")

if [[ -z "${SOURCE_FILE}" || -z "${RESULT_FILE}" ]]; then
  echo "Usage: $0 <input_file> <result_file> [voxel_size or xy,z]" >&2
  exit 1
fi

# Parse VOXEL_SPEC into RES_XY and RES_Z
if [[ "$VOXEL_SPEC" == *","* ]]; then
  export RES_XY="${VOXEL_SPEC%%,*}"
  export RES_Z="${VOXEL_SPEC#*,}"
else
  export RES_XY="$VOXEL_SPEC"
  export RES_Z="$VOXEL_SPEC"
fi

# Validate > 0
awk -v a="$RES_XY" -v b="$RES_Z" 'BEGIN{ if(a<=0 || b<=0) exit 1 }' || {
  echo "[FATAL] voxel sizes must be positive: xy=$RES_XY z=$RES_Z" >&2; exit 2; }

#==============================
# 2) Logging
#==============================
export LOG_FILE="$RESULT_FILE.log"
echo "$(date) node ready" >> "$LOG_FILE"
echo "Scratch directory: $SCRATCHDIR" >> "$LOG_FILE"
echo "Input: $SOURCE_FILE" >> "$LOG_FILE"
echo "Output: $RESULT_FILE" >> "$LOG_FILE"
echo "Voxel: $VOXEL_SPEC  (xy=$RES_XY, z=$RES_Z)" >> "$LOG_FILE"

#==============================
# 3) Prep scratch
#==============================
cd "$SCRATCHDIR"
cp "$SOURCE_FILE" "$SCRATCHDIR/in.laz" &>> "$LOG_FILE"

#==============================
# 4) Singularity env
#==============================
echo "$(date) [INFO] Attempting to load Singularity module..." >> "$LOG_FILE"
if module add singul/ &>> "$LOG_FILE"; then
  echo "$(date) [INFO] Module 'singul/' loaded." >> "$LOG_FILE"
else
  echo "$(date) [WARN] Could not load 'singul/' module." >> "$LOG_FILE"
fi

PDAL_IMG_SOURCE="/storage/brno2/home/krucek/bluecat/singularity_img/pdal.img"
PDAL_IMG_DEST="$SCRATCHDIR/pdal.img"
echo "$(date) [INFO] Copying PDAL image..." >> "$LOG_FILE"
cp "$PDAL_IMG_SOURCE" "$PDAL_IMG_DEST" &>> "$LOG_FILE"

# Workaround for LD_LIBRARY_PATH
if [[ -z "${LD_LIBRARY_PATH:-}" ]]; then
  export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu"
else
  export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:$(echo "$LD_LIBRARY_PATH" | sed 's/::/:/g' | sed 's/:$//')"
fi
echo "$(date) [INFO] LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> "$LOG_FILE"

#==============================
# 5) Build PDAL pipeline JSON
#==============================
INPUT_FILE="in.laz"
OUTPUT_FILE="cloud.laz"

# anisotropic support: scale Z so we can use a single cubic cell (RES_XY)
# zscale = RES_XY / RES_Z
zscale=$(awk -v xy="$RES_XY" -v zz="$RES_Z" 'BEGIN{ printf "%.17g", xy/zz }')
invz=$(awk -v s="$zscale" 'BEGIN{ printf "%.17g", 1.0/s }')

# pre/post Z translation by ±0.5*RES_Z
zshift=$(awk -v zz="$RES_Z" 'BEGIN{ printf "%.17g", 0.5*zz }')

cat > pdal_pipeline.json <<EOF
{
  "pipeline": [
    { "type": "readers.las", "filename": "$INPUT_FILE" },

    // Pre Z shift: -0.5 * resZ
    { "type": "filters.transformation",
      "matrix": "1 0 0 0  0 1 0 0  0 0 1 -$zshift  0 0 0 1" },

    // Pre Z scale (anisotropy): scale Z by (resXY/resZ)
    { "type": "filters.transformation",
      "matrix": "1 0 0 0  0 1 0 0  0 0 $zscale 0  0 0 0 1" },

    // Centroid-of-points per voxel; output = nearest real point to centroid
    { "type": "filters.voxelcentroidnearestneighbor",
      "cell": $RES_XY },

    // Undo Z scale
    { "type": "filters.transformation",
      "matrix": "1 0 0 0  0 1 0 0  0 0 $invz 0  0 0 0 1" },

    // Post Z shift: +0.5 * resZ
    { "type": "filters.transformation",
      "matrix": "1 0 0 0  0 1 0 0  0 0 1 $zshift  0 0 0 1" },

    // Writer (LAS/LAZ inferred by extension)
    { "type": "writers.las",
      "filename": "$OUTPUT_FILE",
      "minor_version": 2,
      "dataformat_id": 1
    }
  ]
}
EOF

#==============================
# 6) Run PDAL
#==============================
singularity exec -B "$SCRATCHDIR":/data "$PDAL_IMG_DEST" \
  pdal pipeline /data/pdal_pipeline.json

#==============================
# 7) Stage back the result
#==============================
cp "$OUTPUT_FILE" "$RESULT_FILE" &>> "$LOG_FILE"
echo "$(date) [OK] wrote $RESULT_FILE" >> "$LOG_FILE"
