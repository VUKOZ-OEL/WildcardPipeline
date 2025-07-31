#!/bin/bash
set -euo pipefail

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S") 
export LOG_FILE="./log_$TIMESTAMP.log"

echo "$(date) node ready" >> $LOG_FILE
echo "Scratch directory: $SCRATCHDIR" >> $LOG_FILE

# Check param
if [ $# -lt 1 ]; then
  echo "Usage: $0 <full/path/to/file.laz>"  >> $LOG_FILE
  exit 1
fi

# Input full path: e.g., /storage/brno2/home/mihazivec/Wildcard/Data/<folder>/<file.laz>
inputPath="$1"
fullFileName=$(basename "$inputPath")      # file.laz
fileName="${fullFileName%.*}"             # file
inputDir=$(dirname "$inputPath")          # full dir path

echo "Full file: $fullFileName" >> $LOG_FILE
echo "Stripped name: $fileName" >> $LOG_FILE

BASEDIR=/storage/brno2/home/mihazivec/Wildcard
voxelSize=0.5
cd $BASEDIR

echo "Singularity path: $(which singularity)" >> $LOG_FILE
#chmod +x Code/CloudHomogenization.sh
#chmod +x Code/TileSplit.sh

#mkdir -p ./Homogenized/$fileName
#./Code/CloudHomogenization.sh \
#    "$inputPath" \
#    "$BASEDIR/Homogenized/$fileName/${fileName}_h20.laz" \
#    $voxelSize

#file_count=$(find "./Homogenized/$fileName" -type f | wc -l)
#echo "Cloud $fileName contains $file_count files" >> $LOG_FILE

#INPUT_FILE="$BASEDIR/Homogenized/$fileName/${fileName}_h20.laz"
#echo "Checking input file for TileSplit: $INPUT_FILE" >> $LOG_FILE
#if [ ! -f "$INPUT_FILE" ]; then
#    echo "ERROR: Homogenized file not found at $INPUT_FILE" >> $LOG_FILE
#else
#    echo "Homogenized file found, proceeding with TileSplit..." >> $LOG_FILE
#fi

#mkdir -p ./Cut/$fileName
#./Code/TileSplit.sh \
#    "$BASEDIR/Homogenized/$fileName/${fileName}_h20.laz" \
#    "$BASEDIR/Cut/$fileName/" \
#    10 #>> $LOG_FILE 2>&1

#file_count1=$(find "./Cut/$fileName" -type f | wc -l)
#echo "Cloud $fileName contains $file_count1 files" >> $LOG_FILE

SIMG="/storage/brno2/home/krucek/bluecat/singularity_img/lidr.img"
#Make directory for processed files
mkdir -p ./Processed/$fileName
#Make directory for results
mkdir -p ./Results/$fileName

echo "Starting preprocessing" >> $LOG_FILE
singularity exec \
   -B "$inputPath:/input" \
   -B "$BASEDIR/Processed/$fileName/:/output" \
   -B "$BASEDIR/Results/$fileName/:/results" \
   -B "$(pwd)/Code:/scripts" \
   $SIMG \
   Rscript /scripts/Preprocessing.R /input /output /results >> $LOG_FILE

file_count2=$(find "./Processed/$fileName" -type f | wc -l)
echo "Cloud $fileName contains $file_count2 files" >> $LOG_FILE

singularity exec \
    -B "$BASEDIR/Processed/$fileName/:/input" \
    -B "$BASEDIR/Processed/$fileName/:/output" \
    -B "$(pwd)/Code:/scripts" \
    $SIMG \
    Rscript /scripts/MaxHeight98.R /input /output $fileName

#mkdir -p ./Results/$fileName
singularity exec \
    -B "$BASEDIR/Processed/$fileName/:/input" \
    -B "$BASEDIR/Results/$fileName/:/output" \
    -B "$BASEDIR/Code:/scripts" \
    $SIMG \
    Rscript /scripts/HeightVariabilityMetrics.R /input /output $fileName >> $LOG_FILE

singularity exec \
    -B "$BASEDIR/Processed/$fileName/:/input" \
    -B "$BASEDIR/Results/$fileName/:/output" \
    -B "$BASEDIR/Processed/$fileName/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$BASEDIR/Code:/scripts" \
    $SIMG \
    Rscript /scripts/GappinesMetrics.R /input /output /maxHeight $fileName >> $LOG_FILE

singularity exec \
    -B "$BASEDIR/Processed/$fileName/:/input" \
    -B "$BASEDIR/Results/$fileName/:/output" \
    -B "$BASEDIR/Processed/$fileName/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$BASEDIR/Code:/scripts" \
    $SIMG \
    Rscript /scripts/CanopyDensity.R /input /output /maxHeight $fileName >> $LOG_FILE

singularity exec \
    -B "$BASEDIR/Processed/$fileName/:/input" \
    -B "$BASEDIR/Results/$fileName/:/output" \
    -B "$BASEDIR/Processed/$fileName/${fileName}_maxHeight.txt:/maxHeight" \
    -B "$BASEDIR/Code:/scripts" \
    $SIMG \
    Rscript /scripts/CanopyVolumeMetrics.R /input /output /maxHeight $fileName $voxelSize >> $LOG_FILE
