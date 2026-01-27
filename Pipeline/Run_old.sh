#!/bin/bash
#PBS -N Wildcard_job_homogenization
#PBS -l select=1:ncpus=5:mem=256gb:scratch_local=300gb
#PBS -l walltime=36:00:00

set -euo pipefail

S3_SUBFOLDER="VUK__6__Eustaska__NA__ULS_h10"
S3_FILENAME="VUK__6__Eustaska__NA__ULS_h10.laz"

POLYGON_SUBFOLDER="Polygon/Eustaska"
POLYGON_FILE="EU_boundaries_cropped_3035.shp"

BASEDIR="/storage/brno2/home/mihazivec/Wildcard"
NOFELEMENTS="10"
VOXELSIZE="0.1"  #Must be written in meters
COORDINATES="EPSG:3035"                                    # EPSG system in use for the data; Number only

cd "$BASEDIR"

LOCAL_TARGET_PATH="$BASEDIR/Data/$S3_SUBFOLDER/$S3_FILENAME"
POLYGON_PATH="$BASEDIR/$POLYGON_SUBFOLDER/$POLYGON_FILE"

echo "SCRATCHDIR=${SCRATCHDIR-<unset>}"
echo "TMPDIR=${TMPDIR-<unset>}"

chmod +x "$BASEDIR/WildcardPipeline_noHomogenization.sh"
"$BASEDIR/WildcardPipeline_noHomogenization.sh" "$LOCAL_TARGET_PATH" "$POLYGON_PATH" "$NOFELEMENTS" "$VOXELSIZE" "$COORDINATES"

