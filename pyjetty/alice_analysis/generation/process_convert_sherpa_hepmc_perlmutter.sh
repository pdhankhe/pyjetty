#! /bin/bash

# This script takes an input HepMC file path as an argument, and runs a python script to
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

if [ "$1" != "" ]; then
  INPUT_FILE=$1
  #echo "Input file: $INPUT_FILE"
else
  echo "Wrong command line arguments"
fi

if [ "$2" != "" ]; then
  JOB_ID=$2
  echo "Job ID: $JOB_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$3" != "" ]; then
  TASK_ID=$3
  echo "Task ID: $TASK_ID"
else
  echo "Wrong command line arguments"
fi

if [ "$4" != "" ]; then
  GENERATION_TYPE=$4
  echo "Generation type: $GENERATION_TYPE"
else
  echo "Wrong command line arguments"
fi


# Load modules
source /global/homes/b/blianggi/pyjetty_env.sh

# Run main script
ALICEANALYSIS_DIR=/global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis
OUTPUT_DIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/sherpagen/$GENERATION_TYPE/tree_gen/$JOB_ID/$TASK_ID"
python $ALICEANALYSIS_DIR/generation/hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g sherpa --no-progress-bar -d --hepmc 3
