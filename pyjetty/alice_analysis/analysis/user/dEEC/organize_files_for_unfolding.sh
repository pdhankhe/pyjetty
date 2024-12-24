#!/bin/bash

cd /global/cfs/cdirs/alice/blianggi/mypyjetty #if perlmutter

cd storage/dEEC/unfolding

PTBINS=("PTBIN0" "PTBIN1" "PTBIN2")
RLBINS=("RLBIN0" "RLBIN1" "RLBIN2" "RLBIN3" "RLBIN4")
OBSERVABLES=("deltap" "deltapt" "deltapl" "charge" "weights")


# make file directories and put 2 root files in each: 
for obs in "${OBSERVABLES[@]}"; do
    echo "OBS: $obs"
    mkdir -p $obs

    for ptbin in "${PTBINS[@]}"; do
        echo "PT BIN: $ptbin"
        mkdir -p $ptbin

        for rlbin in "${RLBINS[@]}"; do
            echo "RL BIN: $rlbin"
            mkdir -p $rlbin

            root -l organize_data_files_for_unfolding.cpp # this is on hiccup...
            root -l organize_mc_files_for_unfolding.cpp # this is on perlmutter...

        done

    done

done

