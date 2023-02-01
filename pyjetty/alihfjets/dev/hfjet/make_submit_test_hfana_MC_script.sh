#!/bin/bash

cpwd=${PWD}

pp5TeVMCfiles=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/files_D0count_pp5TeV_MC569.txt
enrichpp5TeVMCfiles=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/files_D0count_pp5TeV_enhanceMC.txt
flist=${enrichpp5TeVMCfiles}
#flist=${pp5TeVMCfiles}
dname=$(date +"%Y-%m-%d-%H-%M")

flistdname=$(dirname ${flist})
#outputdir=$(basename ${flistdname})
#outputdir=LHC18b8
outputdir=LHC18a4a2
outputdir=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/MC_hfana/${outputdir}/${dname}
mkdir -pv ${outputdir}
cd ${outputdir}
pwd

cp -v ${flist} .
#split --additional-suffix=.flist -d -l 10 -a 5 ${flist}
split --additional-suffix=.flist -d -l 1 -a 5 ${flist}

job_lists=$(find $PWD -name "*.flist" | sort)

cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/test_hfana.sh .
cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/process_mc_hfjet.py .
cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/configcuts_ptbin.yaml .
cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/bitwise.py .

executable=${PWD}/process_mc_hfjet.py
configFile=${PWD}/configcuts_ptbin.yaml
submit_script=${PWD}/submit_all.sh
rm -f ${submit_script}

echo "sbatch --chdir=${PWD}  ${PWD}/test_hfana.sh ${executable} ${configFile} ${PWD}/ ${outputdir}" | tee -a ${submit_script}
chmod +x ${submit_script}

echo "[i] created: ${submit_script}"

cd ${cpwd}
