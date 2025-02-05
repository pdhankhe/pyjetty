#!/bin/bash

# should be sitting in directory hf_corr.
# Input needed: ./make_plot_dirs.sh directory_name -p
# -p is if perly is used, if hiccup, use -h
# this is now set for DATA

if [ "$1" != "" ]; then
  ATTEMPT_DIR=$1
  echo "Directory name to be saved in: $ATTEMPT_DIR"
else
  echo "Wrong command line arguments"
  ATTEMPT_DIR="attempt_sub"
  # or exit??
fi

if [ "$2" == "-p" ]; then
  BASEDIR="/global/cfs/cdirs/alice/blianggi/mypyjetty/storage/dEEC"
  echo "Using Perlmutter base directory"
elif [ "$2" == "-h" ]; then
  BASEDIR="/software/users/blianggi/mypyjetty/storage/dEEC"
  echo "Using hiccup base directory"
elif [ "$2" == "-l" ]; then
  BASEDIR="/Volumes/WORK USB/dEEC/storage"
  echo "Using local base directory"
else
  echo "Exiting - h/p/l not specified"
  exit
fi

pushd .

cd "${BASEDIR}/plots"

mkdir -p "${BASEDIR}/rootfiles/$1"

# mkdir -p plots
# cd plots

mkdir -p $1
cd $1

# PTBINS=("20-40" "40-60" "60-80")
# NORMS=("unnormalized" "self_normalized" "norm_by_jets")
OBSERVABLES=("deltap" "deltapt" "deltapl" "deltajt" "deltajl" "charge" "rc")
# WEIGHTED=("unweighted" "weighted") #don't make yet...

for obs in "${OBSERVABLES[@]}"; do
  echo "OBS: $obs"
  mkdir -p $obs

done

popd
