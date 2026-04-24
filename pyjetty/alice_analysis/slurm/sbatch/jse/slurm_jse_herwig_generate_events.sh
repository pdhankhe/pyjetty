#! /bin/bash

#SBATCH --job-name="generateherwig_jse"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=3:00:00
#SBATCH --array=1-600
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-%A_%a.out

# I want to generate files files with 150K events. I need about 100x the events but I will just do 5x.
# Before, I generated 8000 events per file. --> Going to try 5k events per file. --> 150 files per pt hat bin

JET_PTS=(50 100 200 500)
echo "Number of pT-hat bins: ${#JET_PTS[@]}"

NEVENTS_PER_FILE=5000

PT_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) / 150 + 1 )) # 150 files per pt hat bin
CORE_IN_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) % 150 + 1 )) # 150 files per pt hat bin

CURRENT_PT=${JET_PTS[$((PT_BIN - 1))]}
PT_HAT_MIN=$(( CURRENT_PT * 80 / 100 ))
SEED=$(( ($CORE_IN_BIN - 1) * NEVENTS_PER_FILE + 1111 ))
EVNUM_START=$(( NEVENTS_PER_FILE * CORE_IN_BIN ))

# Load Herwig environment
source /home/blianggi/activate_pyjetty.sh
module load herwig_with_deps

# HERWIG_SCRIPT_MPI="/home/james/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$BIN/LHC_5020_MPI.run"
# OUTDIR="/rstorage/generators/herwig_alice/hepmc/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
HERWIG_SCRIPT_MPI="/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$PT_BIN/LHC_5020_MPI_jse.run"
OUTDIR="/rstorage/generators/herwig_alice/hepmc/$SLURM_ARRAY_JOB_ID/${CURRENT_PT}gev/$CORE_IN_BIN"
mkdir -p $OUTDIR

# Generate events
cd $OUTDIR
echo $PWD
echo "Running Herwig7 with MPI switched on..."
Herwig run $HERWIG_SCRIPT_MPI -d2 -N $NEVENTS_PER_FILE -s $SEED

# Clean up
rm *.log

# Now convert hepmc to root files and delete hepmc for space reasons
module purge
source /home/blianggi/activate_pyjetty.sh
which root

cd $OUTDIR
pushd .
HEPMC_FILE=($(ls -1 *.hepmc))
cd /software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/
srun process_convert_herwig_hepmc.sh ${OUTDIR}/$HEPMC_FILE $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID false true # last two args are D0_tree, JSE_tree

popd
rm $HEPMC_FILE
echo "hepmc file removed"

# Now extract just the jets
NEW_OUTPUT_DIR="/rstorage/generators/herwig_alice/tree_gen/$SLURM_ARRAY_JOB_ID/${CURRENT_PT}gev/$CORE_IN_BIN"
cd $NEW_OUTPUT_DIR
JETS_SCRIPT=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/herwig_make_jets.py
INPUT_FILE=${NEW_OUTPUT_DIR}/AnalysisResultsGen.root
CONFIG=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/config/jse/pp/configcuts_ptbin.yaml
EVNUM_START=$(( NEVENTS_PER_FILE * CORE_IN_BIN ))
python $JETS_SCRIPT -i $INPUT_FILE -c $CONFIG --output-dir $NEW_OUTPUT_DIR --ev-num-base $EVNUM_START

# Move stdout to appropriate folder
cd /rstorage/generators/herwig_alice/tree_gen/${SLURM_ARRAY_JOB_ID}/ #before was hepmc
mkdir -p slurm-output
mv /rstorage/alice/AnalysisResults/blianggi/herwig/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out /rstorage/generators/herwig_alice/tree_gen/${SLURM_ARRAY_JOB_ID}/slurm-output/
