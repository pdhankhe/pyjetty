#! /bin/bash

#SBATCH --job-name="generatepythiaotf_jse"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=6:00:00
#SBATCH --array=1-400
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out
#SBATCH --mem=20GB



# I want to generate files files with 150K events. I need about 100x the events. So I will run 100x this amount for each pt hat bin

JET_PTS=(50 100 200 500)
# NUM_PT_BINS=${#JET_PTS[@]}

NEVENTS_PER_FILE=150000 #100k

PT_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) / 100 + 1 ))
CORE_IN_BIN=$(( (SLURM_ARRAY_TASK_ID - 1) % 100 + 1 ))

CURRENT_PT=${JET_PTS[$((PT_BIN - 1))]}
PT_HAT_MIN=$(( CURRENT_PT * 80 / 100 ))
SEED=$(( ($CORE_IN_BIN - 1) * NEVENTS_PER_FILE + 1111 ))
EVNUM_START=$(( NEVENTS_PER_FILE * CORE_IN_BIN ))

OUTPUT_DIR=/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/${SLURM_ARRAY_JOB_ID}/${CURRENT_PT}gev/${CORE_IN_BIN}
SCRIPT=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/process/user/jse/pythia_quark_gluon_jse.py
CONFIG=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/config/jse/pp/configcuts_ptbin.yaml

cd ~/
source pyjetty_env.sh
cd analysis/

mkdir -p $OUTPUT_DIR
echo "running python $SCRIPT -c $CONFIG --output-dir $OUTPUT_DIR --user-seed $SEED --py-pthatmin $PT_HAT_MIN --py-ecm 5020 --nev $NEVENTS_PER_FILE --ev-num-base $EVNUM_START --pythiaopts HardQCD:all=on,TimeShower:pTmin=0.2"
python $SCRIPT -c $CONFIG --output-dir $OUTPUT_DIR --user-seed $SEED --py-pthatmin $PT_HAT_MIN --py-ecm 5020 --nev $NEVENTS_PER_FILE --ev-num-base $EVNUM_START --pythiaopts HardQCD:all=on,TimeShower:pTmin=0.2


# Move stdout to appropriate folder
SLURM_OUTPUT_DIR=/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/pythia_otf/${SLURM_ARRAY_JOB_ID}/slurm-output

mkdir -p $SLURM_OUTPUT_DIR
mv /global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out $SLURM_OUTPUT_DIR
