#! /bin/bash

# This script takes an input file path as an argument, and runs a python script to 
# process the input file and write an output ROOT file.
# The main use is to give this script to a slurm script.

# Take three command line arguments -- (1) input file path, (2) job ID, (3) task ID, (4) config type
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

if [ "$4" != "" ]; then
  CONFIG_TYPE=$4
  echo "Config type: $CONFIG_TYPE"
else
  echo "Wrong command line arguments"
fi


# Define output path from relevant sub-path of input file
OUTPUT_PREFIX="AnalysisResults/blianggi/EEC/$JOB_ID"
# Note: depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f6-10)
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="/rstorage/alice/$OUTPUT_PREFIX/$OUTPUT_SUFFIX/"
mkdir -p $OUTPUT_DIR
echo "Output dir: $OUTPUT_DIR"


# Load modules
module use /software/users/blianggi/heppy/modules
module load heppy/1.0
module use /software/users/blianggi/pyjetty/modules
module load pyjetty/1.0
module list

# Run python script via pipenv
cd /software/users/blianggi/pyjetty/pyjetty/alihfjets/dev/hfjet/
python process/user/hf_EEC/process_mc_hfjet_EEC.py -c config/hf_EEC/configcuts_ptbin.yaml -f $INPUT_FILE -o $OUTPUT_DIR -ct $CONFIG_TYPE

# Move stdout to appropriate folder
mv /rstorage/alice/AnalysisResults/blianggi/EEC/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/alice/AnalysisResults/blianggi/EEC/${JOB_ID}
