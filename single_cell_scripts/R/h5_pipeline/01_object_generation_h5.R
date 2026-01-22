################################################################################
## Title        : Object generation - RNA only, AIM+ cells - H5 version
##
## Input        : 10x raw and filtered bc_feature matrices from AIM+ sorted
##                single-cell RNAseqas .h5 files
##
## Output       : Integrated object of high quality singlets
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading libraries and data:##################################################

library(Seurat)
library(scDblFinder)
library(SeuratWrappers)
library(tidyverse)

## Environment setup:###########################################################

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

# Common object processing functions
source("R/utilities/object_processing_functions.R")

# Activate Python environment
reticulate::use_virtualenv(Sys.getenv("PYTHON_ENV"))

# Define samples
samples <- c(
  paste("aim", "ID104", c("v01", "v13", "v24"), sep = "_"),
  paste("aim", "ID107", c("v01", "v13", "v24", "v28", "v34", "v37", "v40"), sep = "_"),
  paste("aim", "ID112", c("v01", "v13", "v24"), sep = "_"),
  paste("aim", "ID142", c("v02", "v11b", "v11c"), sep = "_")
)

## Pre-flight check #############################################################

# Check that all files are there
file.check <- function(sample.id){
  
  # Produce file vector for a sample
  files <- paste0(raw_data,
                  sample.id,
                  c("/filtered_contig_annotations.csv",
                    "/raw_feature_bc_matrix.h5",
                    "/sample_filtered_feature_bc_matrix.h5"))
  
  # Check that files are there
  all <- lapply(files, function(x) file.exists(x))
  
  # Flatten the list
  out <- unlist(all) %>% 
    unique()
  
  # If all are TRUE, mark sample as OK, otherwise flag it
  ifelse(length(out) == 1 & out == TRUE,
         paste0(sample.id, ": OK!"),
         paste("Check sample", sample.id))
  
}

# Check samples
sapply(c(samples), file.check)

## Read and filter matrices, create seurat object ##############################

matrix.cleaner <- function(sample.id){
  
  # Takes feature-bc 10x matrices from a sample, returns list of 2 matrices:
  # 1. rna counts
  # 2. metadata
  # where low-quality droplets are removed .
  # These matrices are then used to generate a Seurat object.
  
  print(sample.id)
  
  data.slot <- paste0(raw_data, sample.id, "/")
  
  ### Matrix generation --------------------------------------------------------
  
  # Read data
  raw <- Read10X_h5(paste0(data.slot, "raw_feature_bc_matrix.h5"))
  
  # If the matrix comes from an ADT experiment, keep only RNA data
  if(is.list(raw)){
    
    rna <- raw$`Gene Expression`
    
  } else {
    
    rna <- raw
    
  }
  
  cells <- Read10X_h5(paste0(data.slot, "sample_filtered_feature_bc_matrix.h5"))
  
  # If the matrix comes from an ADT experiment, keep only RNA data
  if(is.list(cells)){
    
    cells <- cells$`Gene Expression`
    
  }
  
  # Modify cell names to later match add.tcr
  name.changer <- function(mtx){
    
    colnames(mtx) <- paste(sample.id, colnames(mtx), sep= "_")
    
    return(mtx)
    
  }
  
  rna <- name.changer(rna)
  
  cells <- name.changer(cells)
  
  # Divide barcodes in those that have been called as cells from 10x 
  # and background
  called_cells <-  colnames(cells)
  
  background <-  setdiff(colnames(rna), called_cells)
  
  ### Metadata generation ------------------------------------------------------
  
  # Define mitochondrial genes
  mtgene <-  grep(pattern = "^MT-", rownames(rna), value = TRUE)
  
  # Create metadata
  md <-  data.frame(log10_nCount_RNA = log10(Matrix::colSums(rna)), 
                    nFeature_RNA = Matrix::colSums(rna > 0), 
                    percent.mt = (Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna))*100
  )
  
  md$drop.class <-  ifelse(rownames(md) %in% called_cells, 'cell', 'background')
  
  # Filter barcodes with no counts
  md <- md %>% filter (log10_nCount_RNA > 0)
  
  # Filter metadata
  cellmd <- md %>% filter(drop.class == "cell")
  
  # Export csv's for quality check
  write.table(cellmd,
              paste0(processed_data, "steps/qc/", sample.id, "_cells_features.csv"),
              quote = F,
              sep = ",")
  
  ### Filter matrix ------------------------------------------------------------
  
  # Define filtering limits
  rna.lower <- log10(1000)
  gen.lower <- 500
  mt.upper <- 7.5
  
  # Find high quality cells
  qc_cells <- cellmd %>% 
    filter(log10_nCount_RNA > rna.lower &
             cellmd$nFeature_RNA > gen.lower &
             cellmd$percent.mt < mt.upper) %>%
    rownames()
  
  # Create raw matrices with high quality cells alone
  cell.rna.raw <-  rna[ ,qc_cells]
  cellmd <-  cellmd[qc_cells, ]
  
  ### Create output ------------------------------------------------------------
  
  out <- list(cell.rna.raw, cellmd)
  
  names(out) <- c("rna", "meta")
  
  return(out)
  
}

# Process matrices
matrices <- lapply(samples, matrix.cleaner)

names(matrices) <- samples

### Seurat objects -------------------------------------------------------------

# Create list of Seurat objects with rna data and metadata
seurat.creator <- function(pt.id, matrix.list){
  
  print(pt.id)
  
  # Create Seurat object with RNA assay
  out <- CreateSeuratObject(counts = matrix.list[[pt.id]]$rna,
                            meta.data = matrix.list[[pt.id]]$meta,
                            assay = "RNA",
                            project = pt.id,
                            min.cells = 3)
  
  # Add sample id as orig.ident
  out@meta.data <- mutate(out@meta.data, orig.ident = pt.id)
  
  return(out)
  
}

obj.list <- lapply(samples,
                   function(x) seurat.creator(x, matrices))

names(obj.list) <- samples

## Normalize, scale and reduce RNA assay #######################################

obj.list <- lapply(obj.list, lognorm.and.scale)

obj.list <- lapply(obj.list, PCA.calculator)

## Clustering and visualization ################################################

# Find Neighbors, clusters and run UMAP
clust.and.umap <- function(seur.obj){
  
  # Clusters cells and calculate UMAP embeddings
  
  set.seed(123)
  
  seur.obj <- FindNeighbors(seur.obj, 
                            dims = 1:15, 
                            verbose = F)
  
  seur.obj <- FindClusters(seur.obj, 
                           resolution = 0.6, 
                           algorithm = 4,
                           verbose = F)
  
  seur.obj <- RunUMAP(seur.obj, 
                      dims = 1:15, 
                      verbose = F)
  
  return(seur.obj)
  
}

obj.list <- lapply(obj.list, clust.and.umap)

# Visually inspect
umap.list <- lapply(obj.list, function(x) DimPlot(x,
                                                  label = T) + 
                      NoLegend() +
                      ggtitle(paste(x$orig.ident)))

cowplot::plot_grid(plotlist = umap.list)

# Save object list for ambient RNA quantification
saveRDS(obj.list, paste0(processed_data, "steps/objects/aim_single_objects.rds"))

## Doublet removal #############################################################

### RNA-based doublet identification -------------------------------------------

rna.dbl.marker <- function(seur.obj){
  
  # Returns object where potential doublets (based on gene expression) are 
  # marked
  
  # Run 1 ----------------------------------------------------------------------
  
  print(paste(unique(seur.obj$orig.ident), "run 1"))
  
  # Set RNA as Default assay
  DefaultAssay(seur.obj) <- "RNA"
  
  #Run scDblFinder
  set.seed(1)
  sce.dbl <-  scDblFinder(GetAssayData(seur.obj, slot = "counts"), 
                          clusters = Idents(seur.obj))
  
  # Add resulting scores and classification back to original object
  seur.obj$scDblFinder.score1 <- sce.dbl$scDblFinder.score
  seur.obj$scDblFinder.class1 <- sce.dbl$scDblFinder.class
  
  # Run 2 ----------------------------------------------------------------------
  
  print(paste(unique(seur.obj$orig.ident), "run 2"))
  
  # Set RNA as Default assay
  DefaultAssay(seur.obj) <- "RNA"
  
  #Run scDblFinder
  set.seed(2)
  sce.dbl2 <-  scDblFinder(GetAssayData(seur.obj, slot = "counts"),
                           clusters = Idents(seur.obj))
  
  # Add resulting scores and classification back to original object
  seur.obj$scDblFinder.score2 <- sce.dbl2$scDblFinder.score
  seur.obj$scDblFinder.class2 <- sce.dbl2$scDblFinder.class
  
  # Run 3 ----------------------------------------------------------------------
  
  print(paste(unique(seur.obj$orig.ident), "run 3"))
  
  # Set RNA as Default assay
  DefaultAssay(seur.obj) <- "RNA"
  
  #Run scDblFinder
  set.seed(3)
  sce.dbl3 <-  scDblFinder(GetAssayData(seur.obj, slot = "counts"),
                           clusters = Idents(seur.obj))
  
  # Add resulting scores and classification back to original object
  seur.obj$scDblFinder.score3 <- sce.dbl3$scDblFinder.score
  seur.obj$scDblFinder.class3 <- sce.dbl3$scDblFinder.class
  
  return(seur.obj)
  
}

obj.list <- lapply(obj.list, rna.dbl.marker)

### TCR-based doublet identification -------------------------------------------

# Collect filtered tcr contig files in one list of data frames
AIM.contig.list <- lapply(samples, function(x)
  read.csv(paste0(raw_data,
                  x, 
                  "/filtered_contig_annotations.csv")))

names(AIM.contig.list) <- samples

# Modify barcode to match barcode in object - i.e. add sample info to cell name
tcr.tab <- lapply(samples, function(x){
  
  AIM.contig.list[[x]] %>% 
    tibble() %>% 
    filter(is_cell == "true" &
             productive == "true") %>% 
    select(barcode,
           chain,
           v_gene, d_gene, j_gene, c_gene,
           cdr1, cdr1_nt,
           cdr2, cdr2_nt,
           cdr3, cdr3_nt,
           umis) %>% 
    mutate(barcode = paste0(x,
                            "_",
                            barcode))
  
}) %>% 
  # Bind to one data frame
  do.call(rbind, .)

# Find cells with > 3 productive TCR chains
dbl.tcr <- tcr.tab %>% 
  group_by(barcode) %>% 
  count() %>% 
  filter(n > 3) %>% 
  pull(barcode)

tcr.dbl.remover <- function(seur.obj){
  
  # Mark cells with > 3 TCR chains
  seur.obj@meta.data <- seur.obj@meta.data %>% 
    tibble:: rownames_to_column(var = "bc") %>% 
    mutate(tcr.selection = ifelse(bc %in% dbl.tcr, "remove", "keep")) %>% 
    tibble::column_to_rownames(var = "bc")
  
  return(seur.obj)
  
}

obj.list <- lapply(obj.list, tcr.dbl.remover)

### Remove doublets ------------------------------------------------------------
# Merge to one object
merged.obj <- merge(obj.list[[1]],
                    unlist(obj.list[2:length(obj.list)]),
                    merge.data = TRUE)

# Use scDblFinder class to call doublets
rna.dbl.count <- merged.obj@meta.data[grep("scDblFinder.class", 
                                           colnames(merged.obj@meta.data), 
                                           value = T)]

rna.dbl.id <- rowSums(rna.dbl.count == "doublet") %>% 
  data.frame(doublet.count = .) %>% 
  mutate(RNA.doublet = ifelse(doublet.count == 0, "singlet", "doublet")) %>% 
  select(RNA.doublet)

merged.obj <- AddMetaData(merged.obj, rna.dbl.id)

# Flag cells according to singlet selection results
merged.obj@meta.data <- merged.obj@meta.data %>% 
  mutate(dbl.id = case_when(RNA.doublet == "doublet" & tcr.selection == "keep" ~ "RNA_doublet",
                            RNA.doublet == "singlet" & tcr.selection == "remove" ~ "TCR_doublet",
                            RNA.doublet == "doublet" & tcr.selection == "remove" ~ "RNA+TCR_doublet",
                            RNA.doublet == "singlet" & tcr.selection == "keep" ~ "singlet"))

# Print table for qc
qc.out <- merged.obj@meta.data %>%
  group_by(orig.ident, dbl.id) %>% 
  dplyr::count() %>% 
  group_by(orig.ident) %>% 
  mutate(total_cells = sum(n)) %>% 
  pivot_wider(names_from = dbl.id,
              values_from = n) %>% 
  mutate(total_doublets = sum(RNA_doublet,
                              TCR_doublet,
                              `RNA+TCR_doublet`,
                              na.rm = T)) %>% 
  mutate(dbl.pct = total_doublets/total_cells) %>% 
  select(orig.ident,
         total_cells,
         RNA_doublet,
         TCR_doublet,
         `RNA+TCR_doublet`,
         total_doublets,
         dbl.pct,
         singlet) %>% 
  rename(Sample = orig.ident,
         `Cells after QC, before doublet removal` = total_cells,
         `RNA only doublets` = RNA_doublet,
         `TCR only doublets` = TCR_doublet,
         `RNA & TCR doublets` = `RNA+TCR_doublet`,
         `Total Doublets` = total_doublets,
         `Pct. Doublets` = dbl.pct, 
         Singlets = singlet)

# Remove NAs from number columns
qc.out[2:ncol(qc.out)] <- apply(qc.out[2:ncol(qc.out)], 2, function(x){
  
  ifelse(is.na(x), 0, x)
  
})

write.table(qc.out,
            paste0(processed_data, 
                   "results/qc/aim_singlet_stats.txt"),
            sep = " ",
            dec = ",",
            row.names = F)

# Remove doublets
Idents(merged.obj) <- "dbl.id"

Raw.singlets <- subset(merged.obj, idents = "singlet")

# Save
saveRDS(Raw.singlets, paste0(processed_data, "steps/objects/aim_raw_singlets.rds"))