#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=addEff
#SBATCH --output=%x.%j.txt
#SBATCH --mem=20G

Rscript ../benchmarking_script.R /nfs/home/students/jbernett/cytof/cytof/data/cytof_workflow_SCE.rds /nfs/home/students/ga89koc/hiwi/cytof/DEComparison/pbmc_benchmarking condition patient_id TRUE FALSE merging1
