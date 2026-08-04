#!/bin/bash

#SBATCH --job-name="RMs_AxB"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=6:00:00
#SBATCH --array=1-715
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out

# THIS IS FOR HICCUP!!!

FILELIST=/software/users/blianggi/mypyjetty/jse/LHC24ppRef_filelist.txt #17874 files

NUM_FILES_PER_JOB=25
TOTAL_FILES=$(wc -l < "$FILELIST")


# Load environment
source /home/blianggi/activate_pyjetty.sh

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
    OUTDIR="/rstorage/alice/AnalysisResults/blianggi/jse/rms/$SLURM_ARRAY_JOB_ID/${OUTPUT_SUBDIR}"
    echo "INFILE:" $INFILE "OUTDIR:" $OUTDIR
    mkdir -p $OUTDIR

    # modify as needed but keep /scratch/u/$USER in front, operate on the node's local /scratch ...
    TEMP_OUTPUT_DIRECTORY="/scratch/u/${USER}/${SLURM_ARRAY_JOB_ID}/${OUTPUT_SUBDIR}"
    mkdir -p ${TEMP_OUTPUT_DIRECTORY}

    cd $TEMP_OUTPUT_DIRECTORY #$OUTDIR
    echo $PWD

    # Find jets
    SCRIPT=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/anchmc_find_jets.py
    echo "Processing line $i: $INFILE"
    python $SCRIPT $INFILE

    # Process jets, make response matrices
    PARQUET_FILE=$TEMP_OUTPUT_DIRECTORY/jets_out.parquet
    RM_SCRIPT=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/make_rms.py
    python $RM_SCRIPT $PARQUET_FILE


    # now copy files/results to /rstorage
    cp -r ${TEMP_OUTPUT_DIRECTORY}/* ${OUTDIR}/
    # in general we'd want to delete the temp dir after job is done
    # make sure the copy was succesful...
    if [ $? -eq 0 ]; then
        echo "copy done - removing the temp output dir."
        rm -rf ${TEMP_OUTPUT_DIRECTORY}
    else
        echo "The copy command failed - not deleting the temp output dir."
    fi
    cd $OUTDIR
    
done

# Move stdout to appropriate folder
cd /rstorage/alice/AnalysisResults/blianggi/jse/rms/${SLURM_ARRAY_JOB_ID}/
mkdir -p slurm-output
mv /rstorage/alice/AnalysisResults/blianggi/jse/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /rstorage/alice/AnalysisResults/blianggi/jse/rms/${SLURM_ARRAY_JOB_ID}/slurm-output/
