#!/bin/bash

#SBATCH --job-name="RMs_AxB"
#SBATCH --nodes=1 --ntasks=1
#SBATCH --account=alice
#SBATCH --qos=shared
#SBATCH --constraint=cpu
#SBATCH --time=24:00:00
#SBATCH --array=0-3
#SBATCH --exclude=nid004104,nid004160,nid004149
#SBATCH --output=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jse/slurm-%A_%a.out


# THIS IS FOR PERLMUTTER!!!

OBJS=(AA AB BB rad)
OBJ=${OBJS[$SLURM_ARRAY_TASK_ID]}

cd ~/
source pyjetty_env.sh

cd /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/analysis/user/jse/3Dunfolding
OUTPUT_TXT_FILE=/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/jse/rootfiles/unfolding/txt_output
root -l -b -q "unfolding_tests.C(\"${OBJ}\")" > ${OUTPUT_TXT_FILE}/unfolding_test_output_${OBJ}.txt 2>&1

