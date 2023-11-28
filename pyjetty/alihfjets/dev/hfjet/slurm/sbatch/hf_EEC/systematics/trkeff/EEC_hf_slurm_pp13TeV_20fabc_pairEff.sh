#! /bin/bash

#SBATCH --job-name="hfMC"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-4000
#SBATCH --output=/rstorage/alice/AnalysisResults/preeti/EEC/slurm-%A_%a.out

#FILE_PATHS='/home/preeti/analysis/pyjetty/pyjetty/alihfjets/dev/hfjet/FileList/files_D0count_pp5TeV_Data.txt'
FILE_PATHS='/rstorage/alice/data/LHC20fabc/file_list.txt'

NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 4000 + 1 ))


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
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  srun EEC_hf_pp13TeV_20fabc_pairEff.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
