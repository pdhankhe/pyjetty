#! /bin/bash
# Script to merge output ROOT files

JOB_ID=13028172
#1161137 MC run 
#1162319 Data run
FILE_DIR="/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/$JOB_ID"
FILES=$( find "$FILE_DIR" -name "*.root" )
echo "Number of files: $(wc -l $FILES)"

OUTPUT_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/blianggi/jetaxis/$JOB_ID
hadd -f -j 20 $OUTPUT_DIR/AnalysisResultsFinal.root $FILES
