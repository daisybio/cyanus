#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_triple_platelets
#SBATCH --output=%x.%j.txt
#SBATCH --mem-per-cpu=60G

Rscript Untitled.R
