#!/bin/bash
#this file is currently only compatible for hiccup.

pushd .

cd histograms_from_tuples/

cd ../../../../
ALICEANALYSIS_DIR=$PWD
# ALICEANALYSIS_DIR="/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis"
cd -

# python ${ALICEANALYSIS_DIR}/slurm/utils/james/scaleHistograms.py -c /software/users/blianggi/mypyjetty/analysis/scalefactors/herwig_scaleFactors.yaml
python ../scaleHistograms_beatrice.py -c /rstorage/generators/pythia_alice/tree_fastsim/scaleFactors.yaml

popd