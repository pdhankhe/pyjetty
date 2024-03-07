#! /bin/bash

#SBATCH --job-name="HerwigHepMC"
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=1-2500
#SBATCH --output=/rstorage/alice/AnalysisResults/blianggi/herwig/slurm-%A_%a.out

# Number of events per pT-hat bin (for statistics)
NEV_DESIRED=2000000

# Lower edges of the pT-hat bins
# PTHAT_BINS=(7 9 12 16 21 28 36 45 57 70 85 99 115 132 150 169 190 212 235 260)
PTHAT_BINS=(5 7 9 12 16 21 28 36 45 57)
echo "Number of pT-hat bins: ${#PTHAT_BINS[@]}"

# Currently we have 8 nodes * 20 cores active
NCORES=5000
NEV_PER_JOB=$(( $NEV_DESIRED * ${#PTHAT_BINS[@]} / $NCORES ))
echo "Number of events per job: $NEV_PER_JOB"
NCORES_PER_BIN=$(( $NCORES / ${#PTHAT_BINS[@]} ))
echo "Number of cores per pT-hat bin: $NCORES_PER_BIN"

BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) / $NCORES_PER_BIN + 1))
CORE_IN_BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) % $NCORES_PER_BIN + 1))
PTHAT_MIN=${PTHAT_BINS[$(( $BIN - 1 ))]}
if [ $BIN -lt ${#PTHAT_BINS[@]} ]; then
	PTHAT_MAX=${PTHAT_BINS[$BIN]}
	echo "Calculating bin $BIN (pThat=[$PTHAT_MIN,$PTHAT_MAX]) with core number $CORE_IN_BIN"
else
	echo "Calculating bin $BIN (pThat_min=$PTHAT_MIN) with core number $CORE_IN_BIN"
fi

SEED=$(( ($CORE_IN_BIN - 1) * NEV_PER_JOB + 1111 ))

# HERWIG_SCRIPT_MPI="/home/james/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$BIN/LHC_5020_MPI.run"
# OUTDIR="/rstorage/generators/herwig_alice/hepmc/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
HERWIG_SCRIPT_MPI="/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/herwig/run/$BIN/LHC_13000_HF_MPI.run"
OUTDIR="/rstorage/generators/herwig_alice/hepmc/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
mkdir -p $OUTDIR

# Load Herwig environment and generate events
# source /software/users/james/herwig/bin/activate
source /home/blianggi/activate_pyjetty.sh
module load herwig_with_deps
cd $OUTDIR
echo $PWD
echo "Running Herwig7 with MPI switched on..."
Herwig run $HERWIG_SCRIPT_MPI -d2 -N $NEV_PER_JOB -s $SEED

# Clean up
rm *.log

# Move stdout to appropriate folder
cd /rstorage/generators/herwig_alice/hepmc/${SLURM_ARRAY_JOB_ID}/
mkdir -p slurm-output
mv /rstorage/alice/AnalysisResults/blianggi/herwig/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out slurm-output/
