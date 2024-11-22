#! /bin/bash
# # Script to merge output ROOT files from all pt-hat bins together, in stages


# JOB_ID=$1
# FILE_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/wenqing/$JOB_ID
# OUTPUT_DIR=/global/cfs/projectdirs/alice/alicepro/hiccup/rstorage/alice/AnalysisResults/wenqing/$JOB_ID
FILE_DIR="histograms_from_tuples"


# # Merge all output files from each pt-hat bin
hadd -f $FILE_DIR/RawHistsAfterScaling.root $FILE_DIR/*/RawHists_*.root 



echo "Done"

