#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytoGLMM
#SBATCH --output=%x.%j.txt
#SBATCH --mem=50G

Rscript ../benchmarking_script.R DataGeneration/cytoGLMM_simulated DEComparison/simulatedCytoGLMM condition patient_id TRUE FALSE
