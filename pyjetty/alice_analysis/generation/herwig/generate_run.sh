#!/bin/bash

BASE_DIR=/software/users/blianggi/mypyjetty/pyjetty/pyjetty/alice_analysis/generation/herwig

for BIN in $(seq 1 20);
do
    echo "Generating bin: $BIN"
    cd $BASE_DIR/run/$BIN
    Herwig read $BASE_DIR/config/$BIN/LHC_13000_HF_MPI.in
    # Herwig read $BASE_DIR/config/$BIN/LHC_13000_MPI.in
done
