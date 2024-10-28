#!/bin/bash

# should be sitting in directory hf_corr.
# Input needed: ./make_plot_dirs.sh directory_name
# this is now set for DATA

if [ "$1" != "" ]; then
  ATTEMPT_DIR=$1
  echo "Directory name to be saved in: $ATTEMPT_DIR"
else
  echo "Wrong command line arguments"
  ATTEMPT_DIR="attempt_sub"
  # or exit??
fi

pushd .

mkdir -p plots
cd plots

mkdir -p $1
cd $1

PTBINS=("20-40" "40-60" "60-80")
NORMS=("unnormalized" "self_normalized" "norm_by_jets")
OBSERVABLES=("deltap" "deltapt" "deltapl" "charge" "chargeratio" "rc" "weights" "weights_vs_deltapt")
WEIGHTED=("unweighted" "weighted") #don't make yet...

for ptbin in "${PTBINS[@]}"; do
    echo "PT BIN: $ptbin"
    mkdir -p $ptbin
    for norm in "${NORMS[@]}"; do
        echo "NORM: $norm"
        mkdir -p $ptbin/$norm

        # mkdir -p $ptbin/individuals
        for obs in "${OBSERVABLES[@]}"; do
        echo "OBS: $obs"
            if [[ ($obs == "charge" || $obs == "chargeratio" || $obs == "rc") && $norm == "self_normalized" ]]; then
                echo "in if statement"
                continue
            fi
            mkdir -p $ptbin/$norm/$obs
            mkdir -p $ptbin/$norm/$obs/individuals
        done
        
    done

done

popd
