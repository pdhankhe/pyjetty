#! /bin/bash

#SBATCH --job-name="processdata_dEECs"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-640
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/dEEC/slurm-%A_%a.out


FILE_PATHS='/software/users/blianggi/mypyjetty/dEEC/filelist_LHC17pq_779.txt'
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=$(( $NFILES / 640 + 1 ))
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
  pwd
  cd /software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/slurm/sbatch/dEEC
  pwd

  # FILE=/global/cfs/projectdirs/alice/alicepro/hiccup$(sed -n "$JOB_N"p $FILE_PATHS)
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  srun dEEC_LHC17pq_hiccup.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
