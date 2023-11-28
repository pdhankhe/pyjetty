#! /bin/bash
# Script to merge output ROOT files

JOB_ID=1272019

#1161137 MC run 
#1162319 Data run
FILE_DIR="/rstorage/alice/AnalysisResults/preeti/EEC/$JOB_ID/"
FILES=$( find "$FILE_DIR" -name "*.root" )
echo "Number of files: $(wc -l $FILES)"

OUTPUT_DIR=/rstorage/alice/AnalysisResults/preeti/EEC/$JOB_ID/
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
