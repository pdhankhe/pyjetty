#!/bin/bash

#SBATCH --job-name="dataprocess_jse"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=6:00:00
#SBATCH --array=1-170         # %50 run 170 tasks, up to 50 concurrent
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out
#SBATCH --mem=20GB

# Load your environment
cd ~/
source pyjetty_env.sh
cd analysis/

FILELIST=/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/run3/data/LHC24_ppref/BerkeleyTrees/tree_list.txt
INFILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FILELIST)
TAG=$(basename $(dirname "$INFILE"))   # or any unique tag
OUTDIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/${TAG}

mkdir -p $OUTDIR

# Find jets
SCRIPT_FIND_JETS=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/data_find_jets.py
PARQUET_OUTFILE=${OUTDIR}/DataJetsForAnalysis.parquet
python ${SCRIPT_FIND_JETS} "$INFILE" "$PARQUET_OUTFILE" --ptmin 50.0

# Process jets
SCRIPT_PROCESS_JETS=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/data_process_jets.py
ROOT_OUTFILE=${OUTDIR}/AnalysisResults.root
python ${SCRIPT_PROCESS_JETS} ${PARQUET_OUTFILE} ${ROOT_OUTFILE} --zcuts 0.1 0.2 --no-maxkt

# Move stdout to appropriate folder
cd /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/ #before was hepmc
mkdir -p slurm-output
mv /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data/${SLURM_ARRAY_JOB_ID}/slurm-output/

