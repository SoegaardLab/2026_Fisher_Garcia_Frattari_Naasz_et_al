################################################################################
## Title        : Protein coding genes extraction
##
## Input        : GTF file containing human gene annotations used for Cell
##                Ranger (Gencode v32) from the 10x Genomics reference 
##                `refdata-gex-GRCh38-2020-A`. The reference can be downloaded 
##                from https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads
##                and should be placed in `additional_data/`
##
## Output       : RDS with protein coding genes for PCA
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading required libraries:##################################################

library(tidyverse)

## Environment setup:###########################################################

# Get data location from environment
processed_data <- Sys.getenv("PROCESSED_DATA")

## Extract protein coding genes

v32 <- rtracklayer::readGFF("additional_data/genes.gtf") 

v32 %>% 
  # Keep only protein coding (removes lncRNA and TCR genes)
  filter(gene_type == "protein_coding") %>% 
  # Pull genes and save
  pull(gene_name) %>% 
  unique() %>% 
  sort() %>%
  saveRDS("additional_data/protein_coding_genes.rds")