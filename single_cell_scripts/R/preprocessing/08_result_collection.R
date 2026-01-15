################################################################################
## Title        : Result collection
##
## Input        : Results from clustering (AIM object), TCR, cluster, and gene
##                expression analysis
##
## Output       : Seurat object with all results collected for figure generation
##                and submission
##
## Author       : Giacomo S Frattari
################################################################################
## Loading required libraries:##################################################

library(Seurat)
library(tidyverse)

## Read and collect various results ############################################

AIM <- readRDS(paste0(processed_data, "steps/objects/aim_annotated_adt.rds"))

# ADT results
adt <- readRDS(paste0(processed_data, "steps/adt/adt_results.rds"))

# TCR results
tcr <- readRDS(paste0(processed_data, "steps/tcr/tcr_data.rds"))

# Cluster proportions
dab <- readRDS(paste0(processed_data, "steps/dab/dab_results.rds"))

# Gene expression
rna <- readRDS(paste0(processed_data, "steps/deg/gene_expression.rds"))

# Collect into final object
AIM@misc <- list(rna = rna,
                 adt = adt,
                 tcr = tcr,
                 dab = dab)

# Add module score to metadata
AIM <- AddMetaData(AIM,
                   AIM@misc$rna$antiviral_cd8_module$module_score,
                   col.name = "antigen_specific_cd8")

# Add barcode information as independent metadata
AIM$barcode <- rownames(AIM@meta.data)

# Ensure that w_clusters (working clusters = annotation indices) are standard
# idents for the object

Idents(AIM) <- "w_clusters"

# Save results
saveRDS(AIM, paste0(processed_data, "results/objects/aim_final.rds"))

