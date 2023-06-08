#! /bin/bash
# Script to merge output ROOT files

JOB_ID=1244121
#1161137 MC run 
#1162319 Data run
FILE_DIR="/rstorage/alice/AnalysisResults/preeti/EEC/$JOB_ID/LHC_16/"
FILES=$( find "$FILE_DIR" -name "*.root" )
echo "Number of files: $(wc -l $FILES)"

OUTPUT_DIR=/rstorage/alice/AnalysisResults/preeti/EEC/$JOB_ID/LHC_16
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
