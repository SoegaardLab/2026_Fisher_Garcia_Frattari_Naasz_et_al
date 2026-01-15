################################################################################
## Title        : Object processing functions
##
## Description  : Collection of functions shared across object processing 
##                scripts
##
## Author       : Giacomo S Frattari
################################################################################
##
## Log normalize and scale #####################################################

# Normalizes, finds 2000 Variable Features and scales RNA assay

lognorm.and.scale <- function(seur.obj){
  
  print(unique(seur.obj$orig.ident))
  
  DefaultAssay(seur.obj) <- "RNA"
  
  seur.obj <- NormalizeData(seur.obj)
  seur.obj <- FindVariableFeatures(seur.obj)
  seur.obj <- ScaleData(seur.obj)
  
  return(seur.obj)
  
}

## PCA calculator ##############################################################

# PCA calculation for RNA/SCT assay
PCA.calculator <- function(seur.obj, 
                           pcs = 50,
                           n.var.features = 2000){
  
  # Calculates PCA based on protein coding genes in the Gencode v32 annotation
  # (used for alignment) without mitochondrial genes
  
  print(unique(seur.obj$orig.ident))
  
  gene_whitelist <- readRDS("additional_data/protein_coding_genes.rds")
  
  # Remove mitochondrial genes
  mt.genes <- grep("MT-", rownames(seur.obj), value = T)
  
  pca.features <- setdiff(gene_whitelist, mt.genes)
  
  # Select number of Variable Features
  var.features <- VariableFeatures(seur.obj)[1:n.var.features]
  
  # Run PCA
  seur.obj <- RunPCA(seur.obj, 
                     features = var.features[var.features %in% pca.features],
                     npcs = pcs)
  
  return(seur.obj)
  
}
