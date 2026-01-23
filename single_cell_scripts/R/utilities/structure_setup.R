################################################################################
## Title        : Structure setup
##
## Description  : Recreates the directory structure required to reproduce the
##                analysis seamlessly in Docker within the mounted directory
##
## Author       : Giacomo S Frattari
################################################################################

library(tidyverse)

# --- 1. Create raw_data directory ---

print("Organizing matrix data")

# Create raw data directories
dir.create("data/raw_data/single_cell_data/sc_data_by_donor",
           recursive = T)

# If data is available
if(file.exists("data/sc_data_by_donor")){
  
  # Move single cell data by donor under raw data
  file.rename("data/sc_data_by_donor", 
              "data/raw_data/single_cell_data/sc_data_by_donor") 
  
  # Split into by-sample directories
  # Identify sample names
  list.files("data/raw_data/single_cell_data/sc_data_by_donor") %>%
    sub("(aim_ID....*_v...{0,1})_.+", "\\1", .) %>%
    unique() %>%
    # For each sample name
    for(sample in .){
      
      # Create a directory
      dir.create(paste0(raw_data, sample))
      
      # Move filtered matrix
      file.rename(paste0(raw_data, sample, "_sample_filtered_feature_bc_matrix.h5"), 
                  paste0(raw_data, "/", sample, "/sample_filtered_feature_bc_matrix.h5"))
      
      # Move raw matrix
      file.rename(paste0(raw_data, sample, "_raw_feature_bc_matrix.h5"), 
                  paste0(raw_data, "/", sample, "/raw_feature_bc_matrix.h5"))
      
      # Move tcr annotations
      file.rename(paste0(raw_data, sample, "_filtered_contig_annotations.csv"), 
                  paste0(raw_data, "/", sample, "/filtered_contig_annotations.csv"))
      
    }
  
}

# --- 1. Create processed_data directory ---

print("Creating processed data structure")

dir.create("data/processed_data/")

# --- 2. Create steps directories ---

# Main directory
dir.create("data/processed_data/steps/")

# Directories with shared names between steps and objects
sapply(c("adt", "dab", "deg", "objects", "qc"), 
       function(x) dir.create(paste0(processed_data, "steps/", x)))

# Steps-specific directories
dir.create(paste0(processed_data, "steps/adt/dsb"))
dir.create(paste0(processed_data, "steps/tcr"))

# --- 3. Create results directories ---

# Main directory
dir.create("data/processed_data/results/")

# Directories with shared names between steps and results
sapply(c("adt", "dab", "deg", "objects", "qc"), 
       function(x) dir.create(paste0(processed_data, "results/", x)))

# Result-specific directories
dir.create(paste0(processed_data, "results/figures"))

# --- 4. Copy the provided final object to its end location - allows seamless
# execution of figure scripts without running the entire preprocessing

print("Copying Seurat object to results/objects (if Seurat object available)")

if(file.exists("data/aim_final.rds")){
  
  file.rename(from = "data/aim_final.rds",
              to = paste0(processed_data, "results/objects/aim_final.rds"))
  
}

print("Setup done!")
