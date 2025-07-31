#!/bin/bash
#PBS -N Wildcard_job                                                                                            
#PBS -l select=1:ncpus=4:mem=60gb:scratch_local=100gb                                                            
#PBS -l walltime=1:00:00

S3_FILENAME="Tls_zofin_34_leafOff_h50.laz"  
# Set base directory (change to your actual full path!)                                                          
BASEDIR="/storage/brno2/home/mihazivec/Wildcard"

# Move to the project directory                                                                                  
cd $BASEDIR

LOCAL_TARGET_DIR="$BASEDIR/Data"
LOCAL_TARGET_PATH="$LOCAL_TARGET_DIR/$S3_FILENAME"

#Run the actual job
chmod +x $BASEDIR/WildcardPipeline_noHomogenization.sh

./WildcardPipeline.sh "$LOCAL_TARGET_PATH"
