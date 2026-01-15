################################################################################
## Title        : Reclustering and annotation - AIM+ cells
##
## Input        : Object with high quality, AIM+ sorted cells, no contaminants
##
## Output       : Annotated AIM+ Seurat object for downstream analysis 
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

# Read object
Raw.singlets <- readRDS(paste0(processed_data, "steps/objects/aim_raw_singlets.rds")) 

# Read high-q cells if not in environment

highq <- readRDS(paste0(processed_data, "steps/qc/highq_cells.rds")) 

## Subset object ###############################################################

Clean_1 <- subset(Raw.singlets, cells = highq)

rm(Raw.singlets)

gc()

# Adjust meta data
Clean_1@meta.data <- Clean_1@meta.data %>% 
  select(orig.ident:nFeature_RNA,
         percent.mt)

## Pre-processing: cell-cycle-score regression, dimensionality reduction #######

### Normalize and reduce RNA assay ---------------------------------------------
Clean_1 <- lognorm.and.scale(Clean_1)

Clean_1 <- PCA.calculator(Clean_1, n.var.features = 1500)

### Cell cycle scoring ---------------------------------------------------------

# Join RNA layers
Clean_1[["RNA"]] <- JoinLayers(Clean_1[["RNA"]])

# Calculate cell scores
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

Clean_1 <- CellCycleScoring(Clean_1, 
                            s.features = s.genes, 
                            g2m.features = g2m.genes, 
                            set.ident = TRUE)

# Visualize
DimPlot(Clean_1, 
        reduction = "pca",
        group.by = "Phase") + 
  ggtitle("AIM clean object - uncorrected")

### Regress out cell cycle score difference ------------------------------------

# Calculate difference
Clean_1$CC.Difference <- Clean_1$S.Score - Clean_1$G2M.Score

# Split RNA asays to allow individual scaling
Clean_1[["RNA"]] <- split(Clean_1[["RNA"]], f = Clean_1$orig.ident)

Clean_1 <- ScaleData(Clean_1, 
                     vars.to.regress = "CC.Difference")

### Calculate PCA on corrected values ------------------------------------------

Clean_1 <- PCA.calculator(Clean_1, n.var.features = 1500)

# Visualize
DimPlot(Clean_1, 
        reduction = "pca",
        group.by = "Phase") + 
  ggtitle("AIM Clean object - corrected")

## Integrate RNA with Harmony ##################################################

DefaultAssay(Clean_1) <- "RNA"

set.seed(123)

AIM <- IntegrateLayers(object = Clean_1, 
                       method = HarmonyIntegration,
                       orig.reduction = "pca",
                       new.reduction = "rna.harmony",
                       verbose = T)

# Rejoin RNA layer
AIM[["RNA"]] <- JoinLayers(AIM[["RNA"]])

## Clustering and visualization ################################################

### Find Neighbors, calculate UMAP ---------------------------------------------

set.seed(123)

AIM <- FindNeighbors(AIM,
                     reduction = "rna.harmony",
                     graph.name = c("rna_snn", "rna_nn"),
                     dims = 1:20,
                     verbose = F)

AIM <- RunUMAP(AIM,
               reduction = "rna.harmony",
               reduction.name = "rna_umap",
               reduction.key = "rnaUMAP_",
               dims = 1:20,
               verbose = F)

### Resolution assessment with clustree ----------------------------------------

# Resolutions to test
res <- seq(0.1, 1.3, 0.1)

# Clustering loop
for(r in res){
  
  print(r)
  
  AIM <- FindClusters(AIM,
                      resolution = r,
                      graph.name = "rna_snn",
                      method = "igraph",
                      algorithm = 4,
                      verbose = T)
                                 
  
}

# Assess tree
cluster.tree <- clustree(AIM, 
                         prefix = "rna_snn_res.",
                         prop_filter = 0.15,
                         show_axis = T) + 
  scale_edge_colour_gradientn(colours = viridis::viridis(256),
                              limits = c(0, 10000),
                              oob = scales::squish)

cluster.tree

# Chosen resolution: 0.7
DimPlot(AIM, 
        reduction = "rna_umap",
        group.by = "rna_snn_res.0.7", 
        label = T) + 
  NoLegend()

DimPlot(AIM, 
        reduction = "rna_umap",
        group.by = "rna_snn_res.0.7",
        split.by = "orig.ident",
        ncol = 4,
        label = T) + 
  NoLegend()

Idents(AIM) <- "rna_snn_res.0.7"

## Adapt metadata ##############################################################

### Note current clusters ------------------------------------------------------

AIM@meta.data <- AIM@meta.data %>% 
  mutate(current_clusters = as.character(rna_snn_res.0.7)) %>% 
  mutate(current_clusters = factor(current_clusters,
                                   levels = c(1:max(as.numeric(current_clusters))))) %>% 
  mutate(seurat_clusters = current_clusters)

### Donor-related metadata -----------------------------------------------------

# Split orig ident information into different variables
AIM@meta.data <- AIM@meta.data %>% 
  separate(orig.ident,
           into = c("cells",
                    "donor",
                    "visit"),
           sep = "_",
           remove = FALSE)

# Add clinical information
AIM@meta.data <- AIM@meta.data %>% 
  mutate(clinical_time_point = case_when(orig.ident %in% c("aim_ID104_v13",
                                                           "aim_ID107_v13",
                                                           "aim_ID112_v13",
                                                           "aim_ID142_v02") ~ "ART",
                                         orig.ident %in% c("aim_ID104_v24",
                                                           "aim_ID107_v24",
                                                           "aim_ID107_v28",
                                                           "aim_ID107_v34",
                                                           "aim_ID107_v37",
                                                           "aim_ID107_v40",
                                                           "aim_ID112_v24",
                                                           "aim_ID142_v11b",
                                                           "aim_ID142_v11c") ~ "ATI",
                                         orig.ident %in% c("aim_ID104_v01",
                                                           "aim_ID107_v01",
                                                           "aim_ID112_v01") ~ "viremic"))

## Annotate clusters ###########################################################

annotation <- c("1. CD8 Differentiating",
                "2. Regulatory T cells",
                "3. CD8 Early Effectors",
                "4. Cytotoxic TM",
                "5. Mixed Terminal\nEffectors/\u03b3\u03b4 T cells",
                "6. CD4 Differentiating",
                "7. CD4 CM",
                "8. CD4 Effectors",
                "9. CD8 Late Effectors",
                "10. Proliferating",
                "11. MAIT cells")

AIM$annotation <- factor(annotation[as.numeric(as.character(AIM$current_clusters))],
                         levels = annotation)

# Define w_clusters (working clusters = index of cluster annotations)
AIM$w_clusters <- factor(as.numeric(gsub("\\. .+", "", AIM$annotation)),
                         levels = c(1:length(annotation)))

# Set w_clusters as default identity
Idents(AIM) <- "w_clusters"

# Visualize
DimPlot(AIM, group.by = "annotation")
DimPlot(AIM, group.by = "w_clusters")

saveRDS(AIM, paste0(processed_data, "steps/objects/aim_annotated.rds"))

## Find RNA markers for clusters ###############################################

# !!! NOTE: this step can take several hours to run !!!

all.markers <- FindAllMarkers(AIM,
                              assay = "RNA",
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              only.pos = T,
                              test.use = "MAST",
                              latent.vars = "donor")

# Save DEG results to use in downstream scripts
saveRDS(all.markers, paste0(processed_data, "steps/deg/aim_markers.rds"))

# Save markers for Supplementary file
write.table(all.markers,
            paste0(processed_data, "results/deg/cluster_markers_aim.txt"),
            dec = ",",
            row.names = F)
