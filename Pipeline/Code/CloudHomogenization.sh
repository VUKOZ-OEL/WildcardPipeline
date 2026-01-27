#!/bin/bash

# Enable strict error handling:
# -e: Exit on error
# -u: Treat unset variables as an error
# -o pipefail: Return non-zero if any part of a pipeline fails
set -euo pipefail



#==============================
# 1. Parse and validate inputs
#==============================

# First argument: input point cloud file (e.g., data.laz)
export SOURCE_FILE=${1:-}
# Second argument: output result file (e.g., output.laz)
export RESULT_FILE=${2:-}
# Third argument: optional voxel size, default is 0.005
export VOXEL_SIZE=${3:-0.005}

# Check that mandatory arguments are provided
if [[ -z "$SOURCE_FILE" || -z "$RESULT_FILE" ]]; then
  echo "Usage: $0 <input_file> <result_file> [voxel_size]" >&2
  exit 1
fi

#==============================
# 2. Setup logging
#==============================

# Define log file based on result file name
export LOG_FILE="$RESULT_FILE.log"
echo "$(date) node ready" >> $LOG_FILE
echo "Scratch directory: $SCRATCHDIR" >> $LOG_FILE

#==============================
# 3. Prepare scratch environment
#==============================

# Move to temporary scratch directory (high-speed disk)
cd "$SCRATCHDIR"

# Log the current directory (for debugging)
cd "$SCRATCHDIR" &>> "$LOG_FILE"

# Copy input LAZ file to scratch and rename it to in.laz
cp "$SOURCE_FILE" "$SCRATCHDIR/in.laz" &>> "$LOG_FILE"

#==============================
# 4. Load Singularity module
#==============================
echo "$(date) [INFO] Attempting to load Singularity module..." >> "$LOG_FILE"

# Load Singularity module — adjust this if the module name is wrong
if module add singul/ &>> "$LOG_FILE"; then
  echo "$(date) [INFO] Module 'singul/' loaded successfully." >> "$LOG_FILE"
else
  echo "$(date) [WARN] Failed to load module 'singul/'. You may need to check the name." >> "$LOG_FILE"
fi

# Copy PDAL Singularity image to scratch directory
# NOTE: update this path if needed (ask Martin for the current one)
PDAL_IMG_SOURCE="/storage/brno2/home/krucek/bluecat/singularity_img/pdal.img"
PDAL_IMG_DEST="$SCRATCHDIR/pdal.img"

echo "$(date) [INFO] Copying PDAL Singularity image to scratch..." >> "$LOG_FILE"
cp "$PDAL_IMG_SOURCE" "$PDAL_IMG_DEST" &>> "$LOG_FILE"
echo "$(date) [INFO] Image copied to: $PDAL_IMG_DEST" >> "$LOG_FILE"

# Workaround for broken LD_LIBRARY_PATH
echo "$(date) [INFO] Sanitizing LD_LIBRARY_PATH..." >> "$LOG_FILE"
if [[ -z "${LD_LIBRARY_PATH:-}" ]]; then
  export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu"
else
  export LD_LIBRARY_PATH="/lib/x86_64-linux-gnu:$(echo "$LD_LIBRARY_PATH" | sed 's/::/:/g' | sed 's/:$//')"
fi
echo "$(date) [INFO] LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH" >> "$LOG_FILE"

#==============================
# 5. Build PDAL pipeline JSON
#==============================

# Set input and output filenames used inside the container
INPUT_FILE="in.laz"
OUTPUT_FILE="cloud.laz"

# Generate a PDAL pipeline JSON file for voxel downsampling
cat <<EOF > pdal_pipeline.json
{
  "pipeline": [
    {
      "type": "readers.las",
      "filename": "$INPUT_FILE"
    },
    {
      "type": "filters.voxeldownsize",
      "cell": $VOXEL_SIZE,
      "mode": "center"
    },
    {
      "type": "writers.las",
      "dataformat_id": 1,
      "minor_version": 2,
      "filename": "$OUTPUT_FILE"
    }
  ]
}
EOF

#==============================
# 6. Run PDAL pipeline
#==============================
# Execute PDAL inside the container with the mounted scratch directory
# -B mounts $SCRATCHDIR to /data inside the container
singularity exec -B "$SCRATCHDIR":/data ./pdal.img pdal pipeline /data/pdal_pipeline.json

#==============================
# 7. Copy result to destination
#==============================

# Copy processed LAZ file back to the desired result path
cp "$OUTPUT_FILE" "$RESULT_FILE" &>> "$LOG_FILE"
