library(CATALYST)
library(diffcyt)
library(uwot)
library(ggplot2)
library(dimRed)
library(RANN)
library(ggvenn)

makeSampleSelection <- function(sce=sce, deselected_samples){
  # All Samples
  all_samples <- as.character(unique(colData(sce)$sample_id))
  
  # Samples to keep
  samples <- all_samples[!all_samples %in% deselected_samples]
  
  # Filter Single Cell Experiment
  sce <- filterSCE(sce, sample_id %in% samples)
  
  return (sce)
}

makePatientSelection <- function(sce, deselected_patients){
  # All patients
  all_patients <- as.character(unique(colData(sce)$patient_id))
  
  # Patients to keep
  patients <- all_patients[!all_patients %in% deselected_patients]
  
  # Filter Single Cell Experiment
  sce <- filterSCE(sce, patient_id %in% patients)
  
  return (sce)
}