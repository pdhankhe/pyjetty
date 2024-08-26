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
  SAVE_D0=$4
  echo "Save the D0: $SAVE_D0"
else
  echo "Wrong command line arguments"
fi


# Define output path from relevant sub-path of input file
# Note: suffix depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f6-8)
echo "OUTPUT_SUFFIX SUPPOSED TO BE:"
echo $OUTPUT_SUFFIX
# OUTPUT_SUFFIX=${TASK_ID}
OUTPUT_DIR="/rstorage/generators/herwig_alice/tree_gen/$JOB_ID/$OUTPUT_SUFFIX"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
# module use /software/users/james/heppy/modules
# module load heppy/1.0
# module use /software/users/james/pyjetty/modules
# module load pyjetty/1.0
# module list
source /home/blianggi/activate_pyjetty.sh
module load herwig_with_deps

# Run main script
cd /software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation
if [ "$SAVE_D0" = true ] ; then
    echo 'Running with saving D0!'
    echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d"
    python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d
else
    echo 'Running inclusive'
    echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar"
    python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar
fi


# Move stdout to appropriate folder
# mkdir -p /rstorage/generators/herwig_alice/tree_gen/${JOB_ID}/slurm-output
# mv /rstorage/generators/herwig_alice/tree_gen/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/generators/herwig_alice/tree_gen/${JOB_ID}/slurm-output/


