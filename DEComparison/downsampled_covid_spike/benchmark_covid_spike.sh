#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_covid_spike
#SBATCH --output=%x.%j.txt
#SBATCH --mem=50G

Rscript ../benchmarking_script.R /nfs/home/students/jbernett/cytof/covid /nfs/home/students/ga89koc/hiwi/cytof/DEComparison/downsampled_covid_spike base_spike patient_id TRUE FALSE
