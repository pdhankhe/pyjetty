#!/bin/bash

#SBATCH --job-name="datachecks_jse"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=2:00:00
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
OUTDIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/data_checks/${SLURM_ARRAY_JOB_ID}/${TAG}

mkdir -p $OUTDIR
OUTFILE=${OUTDIR}/HistsDataCheck.root

SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/data_checks.py
python $SCRIPT "$INFILE" "$OUTFILE"