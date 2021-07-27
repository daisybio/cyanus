#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_covid_spike
#SBATCH --output=%x.%j.txt
#SBATCH --mem=50G

Rscript ../benchmarking_script.R DataGeneration/covid DEComparison/downsampled_covid_spike base_spike patient_id TRUE FALSE
