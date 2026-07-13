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

if [ "$5" != "" ]; then
  JSE_TREE=$5
  echo "JSE_TREE: $JSE_TREE"
else
  echo "Wrong command line arguments"
fi

if [ "$6" != "" ]; then
  HERWIG_LOG_FILE=$6
  echo "Herwig log file: $HERWIG_LOG_FILE"
else
  echo "Wrong command line arguments"
fi

# Define output path from relevant sub-path of input file
# Note: suffix depends on file structure of input file -- need to edit appropriately for each dataset
OUTPUT_SUFFIX=$(echo $INPUT_FILE | cut -d/ -f15-16)
echo "OUTPUT_SUFFIX SUPPOSED TO BE:"
echo $OUTPUT_SUFFIX
# OUTPUT_DIR="/global/cfs/cdirs/alice/blianggi/rstorage/alice/generation/blianggi/herwiggen/tree_gen/$JOB_ID/$OUTPUT_SUFFIX/"
OUTPUT_DIR="/global/cfs/cdirs/alice/alicepro/hiccup/rstorage/alice/generation/blianggi/herwiggen/tree_gen/$JOB_ID/$OUTPUT_SUFFIX/"
echo "Output dir: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR

# Load modules
# module use /software/users/james/heppy/modules
# module load heppy/1.0
# module use /software/users/james/pyjetty/modules
# module load pyjetty/1.0
# module list
# source /global/homes/b/blianggi/pyjetty_env.sh
# module load herwig_with_deps

# Run main script
cd /global/cfs/cdirs/alice/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation

if [ "$SAVE_D0" = true ] ; then
    echo 'Running with saving D0!'
    echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d"
    python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d
    # echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d"
    # python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar -d
else
    echo 'Running inclusive'
    if [ "$JSE_TREE" = true ] ; then
        echo "Running with JSE tree structure!"
        echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --jse --add-herwig-parton -l $HERWIG_LOG_FILE --no-progress-bar"
        python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --jse --add-herwig-parton -l $HERWIG_LOG_FILE --no-progress-bar 
    else
        echo "python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar"
        python hepmc2antuple_tn.py -i $INPUT_FILE -o $OUTPUT_DIR/AnalysisResultsGen.root -g herwig --no-progress-bar
    fi
fi


# # Move stdout to appropriate folder
# mkdir -p /global/cfs/cdirs/alice/blianggi/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${JOB_ID}/slurm-output
# mv /global/cfs/cdirs/alice/blianggi/rstorage/alice/generation/blianggi/herwiggen/tree_gen/slurm-${JOB_ID}_${TASK_ID}.out /global/cfs/cdirs/alice/blianggi/rstorage/alice/generation/blianggi/herwiggen/tree_gen/${JOB_ID}/slurm-output/


