#!/bin/bash -l

#SBATCH --job-name=hfana_MC
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --partition=std
#SBATCH --time=24:00:00
#SBATCH --array=0-1000


id
hostname
date

source scl_source enable devtoolset-7
module use /home/preeti/analysisenvHF/pyjetty/modules
module avail
module load pyjetty/1.0
module list


case_num=$(printf %05d $SLURM_ARRAY_TASK_ID)
case_num=x${case_num}.flist
echo $case_num


set -x
executable=${1}
configFile=${2}
flist=${3}$case_num
outputdir=${4}
output=${outputdir}/$(basename ${flist})

if [ -e ${flist} ]; then
	mkdir -p ${outputdir}
	if [ -d ${outputdir} ]; then
		cd ${outputdir}
		if [ -x ${executable} ]; then
			${HEPPY_DIR}/scripts/pipenv_heppy.sh run ${executable} -c ${configFile} -f ${flist} -o ${output}
		else
			echo "[e] no executable: ${executable}"
		fi
	else
		echo "[e] no outputdir: ${outputdir}"
	fi
else
	echo "[e] no file list: ${flist}"
fi

date
