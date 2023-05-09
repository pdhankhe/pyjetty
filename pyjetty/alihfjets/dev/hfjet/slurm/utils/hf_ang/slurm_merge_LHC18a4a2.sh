#! /bin/bash

#BATCH --job-name=mergehfMC
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=10
#SBATCH --partition=quick
#SBATCH --time=2:00:00
#SBATCH --array=1-20
#SBATCH --output=/rstorage/alice/AnalysisResults/preeti/ang/slurm-%A_%a.out

srun merge_LHC18a4a2.sh $SLURM_ARRAY_JOB_ID $SLURM_ARRAY_TASK_ID
