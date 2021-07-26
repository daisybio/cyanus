#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=pbmc_bencmarking_cytof
#SBATCH --output=%x.%j.txt
#SBATCH --mem=20G

Rscript ../benchmarking_script.R DataGeneration/cytof_workflow_SCE.rds DEComparison/pbmc_benchmarking condition patient_id TRUE FALSE merging1
