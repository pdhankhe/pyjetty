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

# if [ "$2" != "" ]; then
#   JOB_ID=$2
#   echo "Job ID: $JOB_ID"
# else
#   echo "Wrong command line arguments"
# fi

# if [ "$3" != "" ]; then
#   TASK_ID=$3
#   echo "Task ID: $TASK_ID"
# else
#   echo "Wrong command line arguments"
# fi


# Load modules
source /home/blianggi/activate_pyjetty.sh

# Run main script
ALICEANALYSIS_DIR=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis
python $ALICEANALYSIS_DIR/generation/hepmc2antuple_tn.py -i $INPUT_FILE -o ./AnalysisResultsGen.root -g sherpa --no-progress-bar -d --hepmc 3
