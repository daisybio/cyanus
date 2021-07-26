#!/bin/bash
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name=covid_spike_downsample
#SBATCH --output=%x.%j.txt
#SBATCH --mem=50G

Rscript ../benchmarking_script.R DataGeneration/covid_spiked DEComparison/downsampled_covid_spike base_spike patient_id TRUE FALSE
