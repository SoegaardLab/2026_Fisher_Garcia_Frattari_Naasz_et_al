################################################################################
## Title        : Raw objects integration, clustering and cleaning - AIM+ cells
##
## Input        : List of Seurat objects with singlets from AIM+ sorted cells
##
## Output       : Seurat object of AIM+ sorted cells
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading libraries and data:##################################################

library(Seurat)
library(SeuratWrappers)
library(clustree)
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

# Read object
Raw.singlets <- readRDS(paste0(processed_data, "steps/objects/aim_raw_singlets.rds")) 


## Pre-processing: cell-cycle-score regression, dimensionality reduction #######

### Normalize and reduce RNA assay ---------------------------------------------
DefaultAssay(Raw.singlets) <- "RNA"

Raw.singlets <- lognorm.and.scale(Raw.singlets)

Raw.singlets <- PCA.calculator(Raw.singlets, n.var.features = 2000)

### Cell cycle scoring ---------------------------------------------------------

# Join RNA layers
Raw.singlets[["RNA"]] <- JoinLayers(Raw.singlets[["RNA"]])

# Calculate cell scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Raw.singlets <- CellCycleScoring(Raw.singlets, 
                                 s.features = s.genes, 
                                 g2m.features = g2m.genes, 
                                 set.ident = TRUE)

# Visualize
DimPlot(Raw.singlets, 
        reduction = "pca",
        group.by = "Phase") + 
  ggtitle("AIM Raw object - uncorrected")

### Regress out cell cycle score difference ------------------------------------

# Calculate difference
Raw.singlets$CC.Difference <- Raw.singlets$S.Score - Raw.singlets$G2M.Score

# Split RNA asays to allow individual scaling
Raw.singlets[["RNA"]] <- split(Raw.singlets[["RNA"]], f = Raw.singlets$orig.ident)

Raw.singlets <- ScaleData(Raw.singlets, 
                          vars.to.regress = "CC.Difference")

### Calculate PCA on corrected values ------------------------------------------

Raw.singlets <- PCA.calculator(Raw.singlets, n.var.features = 2000)

# Visualize
DimPlot(Raw.singlets, 
        reduction = "pca",
        group.by = "Phase") + 
  ggtitle("AIM Raw object - corrected")


## Integrate RNA with Harmony ##################################################

DefaultAssay(Raw.singlets) <- "RNA"

set.seed(123)

Raw.integrated <- IntegrateLayers(object = Raw.singlets, 
                                  method = HarmonyIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "rna.harmony",
                                  verbose = T)

# Rejoin RNA layer
Raw.integrated[["RNA"]] <- JoinLayers(Raw.integrated[["RNA"]])

## Clustering and visualization ################################################

### Find Neighbors, calculate UMAP ---------------------------------------------

set.seed(123)

Raw.integrated <- FindNeighbors(Raw.integrated,
                                reduction = "rna.harmony",
                                graph.name = c("rna_snn", "rna_nn"),
                                dims = 1:15,
                                verbose = F)

Raw.integrated <- RunUMAP(Raw.integrated,
                          reduction = "rna.harmony",
                          reduction.name = "rna_umap",
                          reduction.key = "rnaUMAP_",
                          dims = 1:15,
                          verbose = F)

### Resolution assessment with clustree ----------------------------------------

# Resolutions to test
res <-seq(0.1, 1.3, 0.1)

# Clustering loop
for(r in res){
  
  print(r)
  
  Raw.integrated <- FindClusters(Raw.integrated,
                                 resolution = r,
                                 graph.name = "rna_snn",
                                 method = "igraph",
                                 algorithm = 4,
                                 verbose = T)
  
}

# Assess tree
cluster.tree <- clustree(Raw.integrated, 
                         prefix = "rna_snn_res.",
                         prop_filter = 0.2,
                         show_axis = T) + 
  scale_edge_colour_gradientn(colours = viridis::viridis(256),
                              limits = c(0, 15000),
                              oob = scales::squish) +
  ggtitle("AIM+ sorted cells")

cluster.tree

# Chosen resolution: 0.8
DimPlot(Raw.integrated, 
        reduction = "rna_umap",
        group.by = "rna_snn_res.0.8", 
        label = T) + 
  NoLegend()

Idents(Raw.integrated) <- "rna_snn_res.0.8"

## Cleaning process ############################################################

### Find RNA markers for clusters ----------------------------------------------

# RNA markers
rna.markers  <- FindAllMarkers(Raw.integrated,
                               assay = "RNA",
                               logfc.threshold = 0.25,
                               min.pct = 0.25,
                               min.cells.group = 1,
                               only.pos = T)

# Export markers for Supplementary file
write.table(rna.markers,
            file = paste0(processed_data, "results/deg/aim_cleaning_genes.txt"),
            dec = ",",
            row.names = F,
            col.names = T)

### QC metrics -----------------------------------------------------------------

#### TCR retrieval -------------------------------------------------------------

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

# Mark cells with TCR data
Raw.integrated@meta.data <- Raw.integrated@meta.data %>% 
  rownames_to_column() %>% 
  mutate(bc = rowname) %>% 
  mutate(tcr = ifelse(bc %in% unique(tcr.tab$barcode), 
                      "tcr", 
                      "no.tcr")) %>% 
  column_to_rownames()

# Collect data for QC
tcr.qc <- Raw.integrated@meta.data %>% 
  group_by(rna_snn_res.0.8, tcr) %>% 
  count() %>% 
  ungroup() %>% 
  pivot_wider(names_from = tcr,
              values_from = n) %>% 
  mutate(`TCR recovery rate` = tcr/(tcr+no.tcr), .keep = "unused")

#### Collect and export QC values ----------------------------------------------

qc.out <- Raw.integrated@meta.data %>% 
  mutate(n = 1) %>% 
  group_by(rna_snn_res.0.8) %>% 
  summarise(n = sum(n),
            med_UMI = median(nCount_RNA),
            med_genes = median(nFeature_RNA),
            med_mt = median(percent.mt),
            iqr_mt = paste0(round(quantile(percent.mt, 0.25), 2),
                            "-",
                            round(quantile(percent.mt, 0.75), 2))) %>% 
  left_join(tcr.qc, .) %>% 
  relocate(n, .after = rna_snn_res.0.8) %>% 
  rename(Cluster = rna_snn_res.0.8,
         `Number of cells` = n,
         `Median number of UMIs` =med_UMI,
         `Median number of genes` = med_genes,
         `Percent mitochondrial genes (median)` = med_mt,
         `Percent mitochondrial genes (IQR)` = iqr_mt)

write.table(qc.out,
            file = paste0(processed_data, "results/qc/aim_cleaning_qc.txt"),
            dec = ",",
            row.names = F,
            col.names = T)

### Annotate cells based on QC criteria ----------------------------------------

Raw.integrated@meta.data <- Raw.integrated@meta.data %>% 
  mutate(qc_selection = case_when(rna_snn_res.0.8 %in% c(1:8, 10, 14) ~ "T-cells",
                                  rna_snn_res.0.8 %in% c(9, 12, 16) ~ "Monocytes",
                                  rna_snn_res.0.8 == 11 ~ "B cells",
                                  # 15 extreme low number of genes, 13 with stress genes
                                  # low TCR retrieval
                                  rna_snn_res.0.8 %in% c(13, 15) ~ "Low quality"))

# Visualize on UMAP
DimPlot(Raw.integrated, 
        group.by = "qc_selection") +
  scale_color_manual(values = rev(scales::hue_pal()(4)))

### Subset Seurat --------------------------------------------------------------

# Adjust meta data
Raw.integrated@meta.data <- Raw.integrated@meta.data %>% 
  select(bc, orig.ident,
         nCount_RNA, nFeature_RNA, log10_nCount_RNA, percent.mt,
         seurat_clusters,
         rna_snn_res.0.3:rna_snn_res.1.3,
         tcr,
         qc_selection)

# Save raw rds as result
saveRDS(Raw.integrated, paste0(processed_data, "results/objects/aim_raw_object.rds"))

# Identify bc of high-quality cells
highq <- Raw.integrated@meta.data %>%
  filter(qc_selection == "T-cells") %>% 
  rownames()

saveRDS(highq, paste0(processed_data, "steps/qc/highq_cells.rds"))
