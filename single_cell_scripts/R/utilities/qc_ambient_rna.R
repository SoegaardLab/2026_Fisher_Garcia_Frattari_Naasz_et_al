################################################################################
## Title        : Ambient RNA quantification
##
## Input        : Raw Seurat objects just after low quality cell removal
##
## Output       : Quantification of ambient RNA contamination per sample
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading required libraries:##################################################

library(Seurat)
library(SoupX)
library(tidyverse)

## Environment setup:###########################################################

# Get data location from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

## AIM #########################################################################

# Read objects
aim <- readRDS(paste0(processed_data, "steps/objects/aim_single_objects.rds"))

# Create empty ambient RNA data frame
rho.tab <- data.frame(orig.ident = c(), rho = c())

for (seur.obj in aim){
  
  # Removes ambient RNA, returns object with corrected raw counts
  
  sample.id <- unique(seur.obj$orig.ident)
  data.slot <- paste0(raw_data, sample.id)
  
  print(sample.id)
  
  # Get clusters from Seurat object without modifying Seurat object
  get_clusters <- function(obj){
    
    return(obj@meta.data[["seurat_clusters"]])
    
  }
  
  clust <- get_clusters(seur.obj)
  
  # Add cluster info to Seurat object
  seur.obj$soup_group <- clust
  
  # Read raw matrix in
  raw <- Read10X(paste0(data.slot, "/raw_feature_bc_matrix"))
  
  # Keep only gene data if ADT assay also present
  if(is.list(raw)){
    
    raw <- raw$`Gene Expression`
    
  }
  
  # Keep only cells that are already in the object (i.e. high quality)
  raw <- raw[rownames(seur.obj),]
  
  # Create Soup Channel, estimate and remove ambient RNA
  sc <- SoupChannel(raw, seur.obj[["RNA"]]$counts)
  
  sc <- setClusters(sc, seur.obj$soup_group)
  
  set.seed(123)
  sc <- autoEstCont(sc, doPlot = T, forceAccept = T)
  
  # Add new row in ambient RNA estimation data frame
  new.rho <- data.frame(orig.id = sample.id, rho = sc$fit$rhoEst)
  
  rho.tab <- rbind(rho.tab, new.rho)
  
}

aim_rho <- rho.tab

# Check the result for ID112 v01 - unstable rho, optimize parameters -----------

for (seur.obj in aim["aim_ID112_v01"]){
  
  # Removes ambient RNA, returns object with corrected raw counts
  
  sample.id <- unique(seur.obj$orig.ident)
  data.slot <- paste0(raw_data, sample.id)
  
  print(sample.id)
  
  # Get clusters from Seurat object without modifying Seurat object
  get_clusters <- function(obj){
    
    return(obj@meta.data[["seurat_clusters"]])
    
  }
  
  clust <- get_clusters(seur.obj)
  
  # Add cluster info to Seurat object
  seur.obj$soup_group <- clust
  
  # Read raw matrix in
  raw <- Read10X(paste0(data.slot, "/raw_feature_bc_matrix"))
  
  # Keep only gene data if ADT assay also present
  if(is.list(raw)){
    
    raw <- raw$`Gene Expression`
    
  }
  
  # Keep only cells that are already in the object (i.e. high quality)
  raw <- raw[rownames(seur.obj),]
  
  # Create Soup Channel, estimate and remove ambient RNA
  sc <- SoupChannel(raw, seur.obj[["RNA"]]$counts)
  
  sc <- setClusters(sc, seur.obj$soup_group)
  
  set.seed(123)
  sc <- autoEstCont(sc, doPlot = T, forceAccept = T, tfidfMin = 0.5)
  
  # Add new row in ambient RNA estimation data frame
  new.rho <- data.frame(orig.id = sample.id, rho = sc$fit$rhoEst)
  
}

# Update estimate
aim_rho[aim_rho$orig.id == "aim_ID112_v01",]$rho <- new.rho$rho

# Collect the results, export for Supplementary File ---------------------------
aim_rho %>% 
  rename(Sample = orig.id,
         `Estimated pct. ambient RNA` = rho) %>% 
  write.table(paste0(processed_data, "results/qc/ambient_rna.txt"),
              dec = ",",
              row.names = F)
