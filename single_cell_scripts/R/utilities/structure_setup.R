################################################################################
## Title        : Structure setup
##
## Description  : Recreates the directory structure required to reproduce the
##                analysis seamlessly in Docker within the mounted directory
##
## Author       : Giacomo S Frattari
################################################################################

# --- 0. If not already present, create raw_data directory ---

if(dir.exists("data/raw_data") == FALSE){
  
  dir.create("data/raw_data/single_cell_data/",
             recursive = T)
  
}

# --- 1. Create processed_data directory ---
dir.create("data/processed_data/")

# --- 2. Create steps directories ---

# Main directory
dir.create("data/processed_data/steps/")

# Directories with shared names between steps and objects
common.dir <- c("adt", "dab", "deg", "objects", "qc")

sapply(common.dir, function(x) dir.create(paste0(processed_data, "steps/", x)))

# Steps-specific directories
dir.create(paste0(processed_data, "steps/adt/dsb"))
dir.create(paste0(processed_data, "steps/tcr"))

# --- 3. Create results directories ---

# Main directory
dir.create("data/processed_data/results/")

# Directories with shared names between steps and results
sapply(common.dir, function(x) dir.create(paste0(processed_data, "results/", x)))

# Result-specific directories
dir.create(paste0(processed_data, "results/figures"))

# --- 4. Copy the provided final object to its end location - allows seamless
# execution of figure scripts without running the entire preprocessing

if(file.exists("data/reproduce_results/aim_final.rds")){
  
  file.copy(from = "data/reproduce_results/aim_final.rds",
            to = paste0(processed_data, "results/objects/"))
  
}
