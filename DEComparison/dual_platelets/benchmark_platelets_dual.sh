#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_dual_platelets
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=70G

Rscript ../benchmarking_script.R /nfs/home/students/l.arend/data/platelets_dual /nfs/home/students/ga89koc/hiwi/cytof/DEComparison/dual_platelets platelets patient_id TRUE FALSE

