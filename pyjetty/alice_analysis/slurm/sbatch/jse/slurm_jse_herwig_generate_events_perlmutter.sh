#! /bin/bash

#SBATCH --job-name="generateherwig_jse"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=5:00:00
#SBATCH --array=1-2000
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-%A_%a.out
#SBATCH --mem=20GB


# I want to generate files files with 150K events. I need about 100x the events.
# Before, I generated 8000 events per file. --> Going to try 5k-->10k-->30k events per file. --> 500 files per pt hat bin --> could shorten time to 3:00:00? (taking ~30 min on perlmutter)

JET_PTS=(50 100 200 500)
echo "Number of pT-hat bins: ${#JET_PTS[@]}"

NEVENTS_PER_FILE=30000

PT_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) / 500 + 1 )) # 500 files per pt hat bin
CORE_IN_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) % 500 + 1 )) # 500 files per pt hat bin

CURRENT_PT=${JET_PTS[$((PT_BIN - 1))]}
# PT_HAT_MIN=$(( CURRENT_PT * 80 / 100 ))
SEED=$(( ($CORE_IN_BIN - 1) * NEVENTS_PER_FILE + 1111 ))
EVNUM_START=$(( NEVENTS_PER_FILE * (CORE_IN_BIN - 1) ))

# Load Herwig environment
source /global/homes/b/blianggi/herwig_pyjetty_env.sh

# HERWIG_SCRIPT_MPI="/home/james/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$BIN/LHC_5020_MPI.run"
# OUTDIR="/rstorage/generators/herwig_alice/hepmc/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
HERWIG_SCRIPT_MPI="/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$PT_BIN/LHC_5360_MPI_jse.run"
OUTDIR="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/hepmc/$SLURM_ARRAY_JOB_ID/${CURRENT_PT}gev/$CORE_IN_BIN"
mkdir -p $OUTDIR

# Generate events
cd $OUTDIR
echo $PWD
echo "Running Herwig7 with MPI switched on..."
Herwig run $HERWIG_SCRIPT_MPI -d2 -N $NEVENTS_PER_FILE -s $SEED

# Clean up
rm *-EvtGen.log #don't remove other log yet

# Now convert hepmc to root files and delete hepmc for space reasons
# module purge
module unload herwig_with_deps
source /global/homes/b/blianggi/pyjetty_env.sh
which root

cd $OUTDIR
pushd .
HEPMC_FILE=($(ls -1 *.hepmc))
HERWIG_LOG_FILE=($(ls -1 *S${SEED}.log)) #i.e. LHC_5020_MPI_jse-S1.log
cd /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/
srun process_convert_herwig_hepmc_perlmutter.sh ${OUTDIR}/$HEPMC_FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID false true ${OUTDIR}/$HERWIG_LOG_FILE # last three args are D0_tree, JSE_tree, herwig log file

popd
rm $HEPMC_FILE
rm $HERWIG_LOG_FILE
echo "hepmc file and log files are removed"

# Now extract just the jets
NEW_OUTPUT_DIR="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/$SLURM_ARRAY_JOB_ID/${CURRENT_PT}gev/$CORE_IN_BIN"
cd $NEW_OUTPUT_DIR
JETS_SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/herwig_make_jets.py
INPUT_FILE=${NEW_OUTPUT_DIR}/AnalysisResultsGen.root
CONFIG=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/config/jse/pp/configcuts_ptbin.yaml
EVNUM_START=$(( NEVENTS_PER_FILE * (CORE_IN_BIN - 1) ))
echo "running jets script: $JETS_SCRIPT -i $INPUT_FILE -c $CONFIG --output-dir $NEW_OUTPUT_DIR --ev-num-base $EVNUM_START"
python $JETS_SCRIPT -i $INPUT_FILE -c $CONFIG --output-dir $NEW_OUTPUT_DIR --ev-num-base $EVNUM_START

# Move stdout to appropriate folder
cd /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${SLURM_ARRAY_JOB_ID}/ #before was hepmc
mkdir -p slurm-output
mv /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${SLURM_ARRAY_JOB_ID}/slurm-output/
