#! /bin/bash

#SBATCH --job-name="hfMC_sys"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-1162
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/EEC/slurm-%A_%a.out

FILE_PATHS='/software/users/blianggi/pyjetty/pyjetty/alihfjets/dev/hfjet/FileList/files_D0count_pp5TeV_enhanceMC_758.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Define the configuration type
CONFIG_TYPE="config4"
echo "Config type: $CONFIG_TYPE"


# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 1162 + 1 ))
echo "Files per job: $FILES_PER_JOB"

STOP=$(( SLURM_ARRAY_TASK_ID * FILES_PER_JOB ))
START=$(( $STOP - $(( $FILES_PER_JOB - 1 )) ))

if (( $STOP > $NFILES ))
then
  STOP=$NFILES
fi

echo "START=$START"
echo "STOP=$STOP"

for (( JOB_N = $START; JOB_N <= $STOP; JOB_N++ ))
do
  echo "Starting script..."
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  srun EEC_hf_LHC18a4a2.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID $CONFIG_TYPE
done
