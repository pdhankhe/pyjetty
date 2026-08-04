#!/bin/bash


#SBATCH --job-name="dataprocess_jse"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=9:00:00
#SBATCH --array=1-170%50        # %50 run 170 tasks, up to 50 concurrent
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out


# Load your environment
cd ~/
source /home/blianggi/activate_pyjetty.sh
cd analysis/

FILELIST=/rstorage/alice/run3/data/LHC24_ppref/BerkeleyTrees/tree_list.txt
INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILELIST)
TAG=$(basename $(dirname "$INFILE"))   # or any unique tag
OUTDIR=/rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/${TAG}

mkdir -p $OUTDIR

TEMP_OUTPUT_DIRECTORY="/scratch/u/${USER}/${SLURM_ARRAY_JOB_ID}/${TAG}"
mkdir -p ${TEMP_OUTPUT_DIRECTORY}

# Find jets
SCRIPT_FIND_JETS=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/data_find_jets.py
PARQUET_OUTFILE=${TEMP_OUTPUT_DIRECTORY}/DataJetsForAnalysis.parquet
python ${SCRIPT_FIND_JETS} "$INFILE" "$PARQUET_OUTFILE" #--ptmin 50.0

# Process jets
SCRIPT_PROCESS_JETS=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/data_process_jets.py
ROOT_OUTFILE=${TEMP_OUTPUT_DIRECTORY}/AnalysisResults.root
python ${SCRIPT_PROCESS_JETS} ${PARQUET_OUTFILE} ${ROOT_OUTFILE} --zcuts 0.1 0.2 --no-maxkt --add-noweight --binning groomed


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


# Move stdout to appropriate folder
cd /rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/ 
mkdir -p slurm-output/ #${SLURM_ARRAY_JOB_ID}
mv /rstorage/alice/AnalysisResults/blianggi/jse/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/slurm-output/ #${SLURM_ARRAY_JOB_ID}

