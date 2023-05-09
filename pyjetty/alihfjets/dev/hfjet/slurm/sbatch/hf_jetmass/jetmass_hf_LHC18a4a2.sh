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
OUTPUT_PREFIX="AnalysisResults/preeti/hf_jetmass/$JOB_ID"
# Note: depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f6-10)
#echo $OUTPUT_SUFFIX
OUTPUT_DIR="/rstorage/alice/$OUTPUT_PREFIX/$OUTPUT_SUFFIX/"
mkdir -p $OUTPUT_DIR
echo "Output dir: $OUTPUT_DIR"

# Load modules
module use /home/preeti/analysis/heppy/modules
module load heppy/1.0
module use /home/preeti/analysis/pyjetty/modules
module load pyjetty/1.0
module list

# Run python script via pipenv
cd /home/preeti/analysis/pyjetty/pyjetty/alihfjets/dev/hfjet/
python process/user/hf_jetmass/process_mc_hfjet_jetmass.py -c config/hf_ang/configcuts_ptbin.yaml -f $INPUT_FILE -o $OUTPUT_DIR

# Move stdout to appropriate folder
mv /rstorage/alice/AnalysisResults/preeti/jetmass/slurm-${JOB_ID}_${TASK_ID}.out /rstorage/alice/AnalysisResults/preeti/jetmass/${JOB_ID}
