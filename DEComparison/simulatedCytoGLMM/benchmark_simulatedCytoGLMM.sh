#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=cytofDE_cytoGLMM
#SBATCH --output=%x.%j.txt
#SBATCH --mem=50G

Rscript ../benchmarking_script.R /nfs/home/students/l.arend/data/cytoGLMM_simulated /nfs/home/students/ga89koc/hiwi/cytof/DEComparison/simulatedCytoGLMM condition patient_id TRUE FALSE
