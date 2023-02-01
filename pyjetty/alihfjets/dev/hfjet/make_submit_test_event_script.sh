#!/bin/bash

cpwd=${PWD}

pp5TeVfiles=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/files_D0count_pp5TeV_Data.txt
#pp5TeVfiles=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/files_pp5TeV_Data_full.txt
enrichpp5TeVMCfiles=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/files_D0count_pp5TeV_enhanceMC.txt
#flist=${pp5TeVfiles}
flist=${enrichpp5TeVMCfiles}

dname=$(date +"%Y-%m-%d-%H-%M")

flistdname=$(dirname ${flist})
#outputdir=$(basename ${flistdname})
outputdir=LHC17pq_MC_event
outputdir=/home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/${outputdir}/${dname}
mkdir -pv ${outputdir}
cd ${outputdir}
pwd

cp -v ${flist} .
split --additional-suffix=.flist -d -l 10 -a 5 ${flist}

job_lists=$(find $PWD -name "*.flist" | sort)

cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/test_hfana.sh .
#cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/process_data_hfjet_sparse.py .
cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/process_events.py .
cp -v /home/preeti/analysisenvHF/pyjetty/pyjetty/alihfjets/dev/hfjet/DmesonJet/configcuts_ptbin.yaml .

#executable=${PWD}/process_data_hfjet_sparse.py
executable=${PWD}/process_events.py
configFile=${PWD}/configcuts_ptbin.yaml
submit_script=${PWD}/submit_all.sh
rm -f ${submit_script}

echo "sbatch --chdir=${PWD}  ${PWD}/test_hfana.sh ${executable} ${configFile} ${PWD}/ ${outputdir}" | tee -a ${submit_script}
chmod +x ${submit_script}

echo "[i] created: ${submit_script}"

cd ${cpwd}
