#!/bin/bash
#PBS -N Wildcard_job
#PBS -l select=1:ncpus=5:mem=32gb:scratch_local=100gb
#PBS -l walltime=2:30:00
#PBS -N Wildcard_job
#PBS -v AWS_SHARED_CREDENTIALS_FILE=/storage/brno2/home/mihazivec/.aws/credentials

# ------------------------------------------
# MODIFY THESE TWO VARIABLES BEFORE EACH RUN
S3_SUBFOLDER="UNIUD"                                         # Folder in wp2 bucket
S3_FILENAME="UNIUD_1_Cansiglio_AF-NM_ULS.laz"                # File to process
# ------------------------------------------

# Set base directory
BASEDIR="/storage/brno2/home/mihazivec/Wildcard"
cd $BASEDIR

# Create full S3 path and local target
S3_PATH_REL="$S3_SUBFOLDER/$S3_FILENAME"
FILEBASE="${S3_FILENAME%.*}"                                  # Remove .laz extension
LOCAL_TARGET_DIR="$BASEDIR/Data/$FILEBASE"
LOCAL_TARGET_PATH="$LOCAL_TARGET_DIR/$S3_FILENAME"


# Create data directory
mkdir -p "$LOCAL_TARGET_DIR"

# Download file from S3
echo "Downloading from S3: s3://wp2/$S3_PATH_REL"
#list all files
s5cmd --profile s3cmd_access \
    --endpoint-url=https://s3.cl4.du.cesnet.cz \
    ls s3://wp2/VUK/

# Check if file exists
s5cmd \
  --profile s3cmd_access \
  --endpoint-url=https://s3.cl4.du.cesnet.cz \
  ls "s3://wp2/$S3_PATH_REL" > /dev/null

#exit 0

if [ $? -ne 0 ]; then
  echo "ERROR: File not found in S3: s3://wp2/$S3_PATH_REL"
  echo "Please verify the folder and filename (case-sensitive)."
  exit 2
fi
# ------------------------------------------------------------

# Proceed to download
echo "File found. Downloading..."
s5cmd \
  --profile s3cmd_access \
  --endpoint-url=https://s3.cl4.du.cesnet.cz \
  cp "s3://wp2/$S3_PATH_REL" "$LOCAL_TARGET_PATH"

# Run the actual pipeline
chmod +x $BASEDIR/WildcardPipeline.sh
./WildcardPipeline.sh "$LOCAL_TARGET_PATH"
