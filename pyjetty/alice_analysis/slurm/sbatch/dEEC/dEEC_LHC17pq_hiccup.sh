#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

# Take three command line arguments -- (1) input file path, (2) job ID, (3) task ID
if [ "$1" != "" ]; then
  INPUT_FILE=$1
  echo "Input file: $INPUT_FILE"
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

# Define output path from relevant sub-path of input file
OUTPUT_BASEPATH="/rstorage/alice"
OUTPUT_PREFIX="AnalysisResults/blianggi/dEEC/$JOB_ID"
# Note: depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f7-9) # results in i.e. "LHC17q_FAST/000282367/0085" #results in i.e. "child_1/2029"
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="$OUTPUT_BASEPATH/$OUTPUT_PREFIX/$OUTPUT_SUFFIX/"
mkdir -p $OUTPUT_DIR
echo "Output dir: $OUTPUT_DIR"

# Load modules
source /home/blianggi/activate_pyjetty.sh
module list

# Run python script via pipenv
cd /software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/
python process/user/dEEC/process_data_dEEC.py -c config/dEEC/pp/process_pp.yaml -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mkdir -p $OUTPUT_BASEPATH/$OUTPUT_PREFIX/slurm-output
mv $OUTPUT_BASEPATH/AnalysisResults/blianggi/dEEC/slurm-${JOB_ID}_${TASK_ID}.out $OUTPUT_BASEPATH/$OUTPUT_PREFIX/slurm-output
