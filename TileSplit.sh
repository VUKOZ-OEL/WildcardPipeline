#!/bin/bash
set -e
set -u
set -o pipefail

# Input arguments
SOURCE_FILE=$1           # Path to homogenized .laz file
OUTPUT_DIR=$2            # Directory to store tile outputs
TILE_SIZE=$3

# Create working directory
SCRATCHDIR=$(mktemp -d)
cd "$SCRATCHDIR"

# Copy input and prepare folders
cp "$SOURCE_FILE" in.laz
mkdir tiles

# Download and extract LAStools
wget -q https://downloads.rapidlasso.de/LAStools.tar.gz
tar -xvzf LAStools.tar.gz > /dev/null
rm LAStools.tar.gz

# Index input file
./bin/lasindex64 -i in.laz

# Perform tiling
echo "Tiling:" >> $LOG_FILE
./bin/lastile64 -i in.laz -odir tiles -tile_size $TILE_SIZE -buffer 5 -cores 5 -demo

# Validate output
shopt -s nullglob
las_files=(tiles/*.las)
if [ ${#las_files[@]} -eq 0 ]; then
    echo "No LAS tiles were created — check input file and LAStools output" >&2
    exit 1
fi

# Copy to destination
mkdir -p "$OUTPUT_DIR"
cp "${las_files[@]}" "$OUTPUT_DIR"

# Cleanup
rm -rf "$SCRATCHDIR"
