#! /bin/bash

#SBATCH --job-name="generateherwig_jse_justjets"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=1:00:00
#SBATCH --array=1-2000
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-%A_%a.out
#SBATCH --mem=20GB


# I want to generate files files with 150K events. I need about 100x the events.
# Before, I generated 8000 events per file. --> Going to try 5k-->10k-->30k events per file. --> 500 files per pt hat bin --> could shorten time to 3:00:00? (taking ~30 min on perlmutter)

JET_PTS=(50 100 200 500)
echo "Number of pT-hat bins: ${#JET_PTS[@]}"

NEVENTS_PER_FILE=30000

FILL_IN_JOBID=55293842

PT_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) / 500 + 1 )) # 500 files per pt hat bin
CORE_IN_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) % 500 + 1 )) # 500 files per pt hat bin

CURRENT_PT=${JET_PTS[$((PT_BIN - 1))]}
# PT_HAT_MIN=$(( CURRENT_PT * 80 / 100 ))
SEED=$(( ($CORE_IN_BIN - 1) * NEVENTS_PER_FILE + 1111 ))
EVNUM_START=$(( NEVENTS_PER_FILE * (CORE_IN_BIN - 1) ))


OUTDIR="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/hepmc/$FILL_IN_JOBID/${CURRENT_PT}gev/$CORE_IN_BIN"
# mkdir -p $OUTDIR


source /global/homes/b/blianggi/pyjetty_env.sh
which root

cd $OUTDIR
# pushd .

cd /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/

# Now extract just the jets
NEW_OUTPUT_DIR="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/$FILL_IN_JOBID/${CURRENT_PT}gev/$CORE_IN_BIN"
cd $NEW_OUTPUT_DIR
JETS_SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/herwig_make_jets.py
INPUT_FILE=${NEW_OUTPUT_DIR}/AnalysisResultsGen.root
CONFIG=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/config/jse/pp/configcuts_ptbin.yaml
EVNUM_START=$(( NEVENTS_PER_FILE * (CORE_IN_BIN - 1) ))
echo "running jets script: $JETS_SCRIPT -i $INPUT_FILE -c $CONFIG --output-dir $NEW_OUTPUT_DIR --ev-num-base $EVNUM_START"
python $JETS_SCRIPT -i $INPUT_FILE -c $CONFIG --output-dir $NEW_OUTPUT_DIR --ev-num-base $EVNUM_START

# Move stdout to appropriate folder
cd /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${FILL_IN_JOBID}/ #before was hepmc
mkdir -p slurm-output
mv /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${FILL_IN_JOBID}/slurm-output/
