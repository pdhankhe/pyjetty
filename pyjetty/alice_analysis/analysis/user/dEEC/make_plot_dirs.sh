#!/bin/bash

# should be sitting in directory hf_corr.
# Input needed: ./make_plot_dirs.sh directory_name

if [ "$1" != "" ]; then
  ATTEMPT_DIR=$1
  echo "Directory name to be saved in: $ATTEMPT_DIR"
else
  echo "Wrong command line arguments"
  ATTEMPT_DIR="thirdattempt"
  # or exit??
fi

pushd .

mkdir -p plots
cd plots

mkdir -p $1
cd $1

PTBINS=("20-40" "40-60" "60-80")
NORMS=("unnormalized" "self_normalized" "norm_by_jets")
OBSERVABLES=("deltap" "deltapt" "deltapl" "charge" "chargeratio" "energyweights" "deltaptvsEW")
WEIGHTED=("unweighted" "weighted") #don't make yet...

for ptbin in "${PTBINS[@]}"; do
    mkdir -p $ptbin
    for norm in "${NORMS[@]}"; do
        mkdir -p $ptbin/$norm
    done

    mkdir -p $ptbin/individuals
    for obs in "${OBSERVABLES[@]}"; do
        mkdir -p $ptbin/individuals/$obs
        for norm in "${NORMS[@]}"; do
            if [[ ($obs == "charge" || $obs == "chargeratio") && $norm == "self_normalized" ]]; then
                continue
            fi
            mkdir -p $ptbin/individuals/$obs/$norm
        done
    done
done

popd
