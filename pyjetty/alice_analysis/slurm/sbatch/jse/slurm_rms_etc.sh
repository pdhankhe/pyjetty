#!/bin/bash

#SBATCH --job-name="RMs_AxB"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=6:00:00
#SBATCH --array=1-715
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out
#SBATCH --mem=20GB


# THIS IS FOR PERLMUTTER!!!

FILELIST=/global/cfs/cdirs/alice/blianggi/mypyjetty/jse/LHC24ppRef_filelist.txt #17874 files, on perly... 8826 files?

NUM_FILES_PER_JOB=13 #25
TOTAL_FILES=$(wc -l < "$FILELIST")


# Load environment
source /global/homes/b/blianggi/pyjetty_env.sh

TASK_ID=$(( SLURM_ARRAY_TASK_ID - 1 ))
START=$(( TASK_ID * NUM_FILES_PER_JOB + 1 ))
STOP=$(( START + NUM_FILES_PER_JOB - 1 ))
if [ "$STOP" -gt "$TOTAL_FILES" ]; then
    STOP=$TOTAL_FILES
fi
echo "START:" $START "STOP:" $STOP

# from my testing, takes <5 minutes each to run
for i in $(seq "$START" "$STOP")
do

    # Define input/output directories
    INFILE=$(sed -n "${i}p" "$FILELIST")
    OUTPUT_SUBDIR=$(echo "$INFILE" | awk -F/ '{print $(NF-3)"/"$(NF-2)"/"$(NF-1)}')
    OUTDIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/$SLURM_ARRAY_JOB_ID/${OUTPUT_SUBDIR}"
    echo "INFILE:" $INFILE "OUTDIR:" $OUTDIR
    mkdir -p $OUTDIR

    cd $OUTDIR
    echo $PWD

    # Find jets
    SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/anchmc_find_jets.py
    echo "Processing line $i: $INFILE"
    python $SCRIPT $INFILE

    # Process jets, make response matrices
    PARQUET_FILE=$OUTDIR/jets_out.parquet
    RM_SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/make_rms_etc.py
    python $RM_SCRIPT $PARQUET_FILE
    
done

# Move stdout to appropriate folder
cd /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/${SLURM_ARRAY_JOB_ID}/
mkdir -p slurm-output
mv /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/rms/${SLURM_ARRAY_JOB_ID}/slurm-output/
