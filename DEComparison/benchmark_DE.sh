#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=cytofDE_first_benchmark.txt
#SBATCH --job-name=cytofDE_first_benchmark

Rscript benchmarking_script.R /nfs/home/students/l.arend/data/covid_spiked/downsampled_files /nfs/home/students/ga89koc/hiwi/cytof/DEComparison/first_benchmark TRUE FALSE
