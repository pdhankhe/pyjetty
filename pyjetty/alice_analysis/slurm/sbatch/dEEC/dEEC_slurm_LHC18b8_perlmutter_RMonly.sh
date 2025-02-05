#! /bin/bash

#SBATCH --job-name="processMC_dEECs"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=12:00:00
#SBATCH --array=1-437
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/slurm-%A_%a.out
#SBATCH --mem=16GB

# FILE_PATHS='/global/cfs/projectdirs/alice/blianggi/mypyjetty/dEEC/filelist_LHC18b8_charge_804.txt' #pass 1 mc production
FILE_PATHS='/global/cfs/projectdirs/alice/blianggi/mypyjetty/dEEC/filelist_LHC23a3_806.txt' #using pass 2 version here
# FILE_PATHS='/global/cfs/projectdirs/alice/blianggi/mypyjetty/dEEC/replacement_filelist_LHC23a3_806.txt' #using pass 2 version here
NFILES=$(wc -l < $FILE_PATHS) #4362 files in pass2
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
# FILES_PER_JOB=1 #1-679
FILES_PER_JOB=10 #$(( $NFILES / 640 + 1 )) #array 1-450
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
  cd /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/slurm/sbatch/dEEC
  pwd
  
  # FILE=/global/cfs/projectdirs/alice/alicepro/hiccup$(sed -n "$JOB_N"p $FILE_PATHS)
  FILE=$(sed -n "$JOB_N"p $FILE_PATHS)
  srun dEEC_LHC18b8_perlmutter_RMonly.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done
