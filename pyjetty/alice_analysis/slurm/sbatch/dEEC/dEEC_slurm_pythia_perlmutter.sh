#! /bin/bash

#SBATCH --job-name="processpythia_dEECs"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=2:00:00
#SBATCH --array=1-1000
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/slurm-%A_%a.out
#SBATCH --mem=20GB

FILE_PATHS='/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/generators/pythia_alice/tree_fastsim/1143757/files.txt' #5 TeV
# FILE_PATHS='/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_fastsim/47946306/47835863/files.txt' #5 TeV, CTEQ 5L PDF
# FILE_PATHS='/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/pythiagen/tree_gen/26652369/files.txt' #13 TeV
NFILES=$(wc -l < $FILE_PATHS)
echo "N files to process: ${NFILES}"

# Currently we have 8 nodes * 20 cores active
FILES_PER_JOB=5 #$(( $NFILES / 640 + 1 ))
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
  srun dEEC_pythia_perlmutter.sh $FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
done


# Move stdout to appropriate folder
OUTPUT_BASEPATH="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice"
OUTPUT_PREFIX="AnalysisResults/blianggi/dEEC/$SLURM_ARRAY_JOB_ID"

mkdir -p $OUTPUT_BASEPATH/$OUTPUT_PREFIX/slurm-output
mv $OUTPUT_BASEPATH/AnalysisResults/blianggi/dEEC/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $OUTPUT_BASEPATH/$OUTPUT_PREFIX/slurm-output
