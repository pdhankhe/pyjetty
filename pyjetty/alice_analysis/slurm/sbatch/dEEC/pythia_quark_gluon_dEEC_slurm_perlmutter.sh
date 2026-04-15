#!/bin/bash
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --account=alice
#SBATCH --job-name=pythiagen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --array=1-300
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/slurm-%A_%a.out
#SBATCH --exclude=nid004104,nid004160,nid004149


# Center of mass energy in GeV
ECM=5020  #13000

# Number of events per pT-hat bin (for statistics)
NEV_DESIRED=10000000 # 10 M per pt-hat bin

# Lower edges of the pT-hat bins
#PTHAT_BINS=(5 7 9 12 16 21 28 36 45 57 70 85 99 115 132 150 169 190 212 235)
# PTHAT_BINS=(5 7 9 12 16 21 28 36 45 57)
TARGET_PTHAT_BINS=(100 200 500)
PTHAT_MIN_BINS=(80 160 400)
# PTHAT_MAX_BINS=(120 240 600)
echo "Number of pT-hat bins: ${#TARGET_PTHAT_BINS[@]}"

# 100 cores per pt-hat bin, w/ 3 pt-hat bins
NCORES=300 
NEV_PER_JOB=$(( $NEV_DESIRED * ${#TARGET_PTHAT_BINS[@]} / $NCORES ))
echo "Number of events per job: $NEV_PER_JOB"
NCORES_PER_BIN=$(( $NCORES / ${#TARGET_PTHAT_BINS[@]} ))
echo "Number of cores per pT-hat bin: $NCORES_PER_BIN"

BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) / $NCORES_PER_BIN + 1))
CORE_IN_BIN=$(( ($SLURM_ARRAY_TASK_ID - 1) % $NCORES_PER_BIN + 1))
# PTHAT_MIN=${PTHAT_BINS[$(( $BIN - 1 ))]}
PTHAT_MIN=${PTHAT_MIN_BINS[$BIN - 1]}

# PTHAT_MAX=${PTHAT_MAX_BINS[$BIN - 1]}
echo "Calculating bin $BIN (pThat=[$PTHAT_MIN,infinity]) with core number $CORE_IN_BIN"


SEED=$(( ($CORE_IN_BIN - 1) * NEV_PER_JOB + 1111 ))

# Do the PYTHIA simulation & matching
OUTDIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/$SLURM_ARRAY_JOB_ID/$BIN/$CORE_IN_BIN"
mkdir -p $OUTDIR

source /global/homes/b/blianggi/pyjetty_env.sh
echo "python is" $(which python)

BASEDIR="/global/cfs/cdirs/alice/blianggi/mypyjetty"
cd ${BASEDIR}/analysis/
SCRIPT="${BASEDIR}/pyjetty/pyjetty/alice_analysis/analysis/user/dEEC/new_code/pythia_quark_gluon_EEC_inclusive_foremily.py"
CONFIG="${BASEDIR}/pyjetty/pyjetty/alice_analysis/config/dEEC/pp/configcuts_ptbin.yaml"


echo "pipenv run python $SCRIPT -c $CONFIG --output-dir $OUTDIR --user-seed $SEED --py-pthatmin $PTHAT_MIN --py-ecm $ECM --nev $NEV_PER_JOB --pythiaopts HardQCD:all=on,TimeShower:pTmin=0.2"
pipenv run python $SCRIPT -c $CONFIG --output-dir $OUTDIR --user-seed $SEED \
    --py-pthatmin $PTHAT_MIN --py-ecm $ECM --nev $NEV_PER_JOB \
    --pythiaopts HardQCD:all=on,TimeShower:pTmin=0.2 --weightON 1

SLURM_OUTDIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/$SLURM_ARRAY_JOB_ID/slurm-output
mkdir -p $SLURM_OUTDIR
mv /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/dEEC/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $SLURM_OUTDIR/