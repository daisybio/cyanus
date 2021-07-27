#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_dual_platelets
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=70G

Rscript ../benchmarking_script.R DataGeneration/platelets_dual DEComparison/dual_platelets platelets patient_id TRUE FALSE

