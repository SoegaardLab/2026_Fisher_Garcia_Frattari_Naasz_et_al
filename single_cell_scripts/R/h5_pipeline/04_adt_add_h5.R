################################################################################
## Title        : ADT add - H5 version
##
## Input        : Annotated AIM object with RNA assay
##
## Output       : Final AIM object with RNA and ADT assay (where available)
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading required libraries:##################################################

library(Seurat)
library(tidyverse)
library(dsb)

## Environment setup:###########################################################

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

# Define samples
samples <- c(
  paste("aim", "ID107", c("v28", "v34", "v37"), sep = "_"),
  paste("aim", "ID142", c("v02", "v11b", "v11c"), sep = "_")
)

# Read object
AIM <- readRDS(paste0(processed_data, "steps/objects/aim_annotated.rds"))

## Read and process matrices ###################################################

matrix.cleaner <- function(sample.id){
  
  # Takes feature-bc 10x matrices from a sample, returns list of 4 matrices:
  # 1. rna counts
  # 2. positive (cells) and 3. negative (background) adt counts to use for adt
  # count correction
  # 4. metadata
  # where low-quality droplets are removed .
  # These matrices are then used to generate a seurat object.
  
  print(sample.id)
  
  data.slot <- paste0(raw_data, sample.id, "/")
  
  ## Matrix generation ---------------------------------------------------------
  
  # Read data in
  raw <-  Read10X_h5(paste0(data.slot, "raw_feature_bc_matrix.h5"))
  
  cells <- Read10X_h5(paste0(data.slot, "sample_filtered_feature_bc_matrix.h5"))
  
  # Modify cell names to later match add.tcr
  name.changer <- function(mtx){
    
    colnames(mtx) <- paste(sample.id, colnames(mtx), sep= "_")
    
    return(mtx)
    
  }
  
  raw <- lapply(raw, name.changer)
  
  cells <- lapply(cells, name.changer)
  
  # Divide barcodes in those that have been called as cells from 10x 
  # and background
  called_cells <-  colnames(cells$`Gene Expression`)
  
  background <-  setdiff(colnames(raw$`Gene Expression`), called_cells)
  
  # Create raw matrices
  adt <-  raw$`Antibody Capture`
  
  rna <-  raw$`Gene Expression`
  
  # Correct names of ADT features
  rownames(adt) <- gsub("Hu(Ms)*(R)*(t)*\\.", 
                        "", 
                        rownames(adt))
  
  rownames(adt) <- gsub("([[:digit:]])(_.+)", 
                        "\\1", 
                        rownames(adt))
  
  ## Create metadata -----------------------------------------------------------
  
  # Define mitochondrial genes
  mtgene <-  grep(pattern = "^MT-", rownames(rna), value = TRUE)
  
  # Create metadata
  md <-  data.frame(
    log10_nCount_RNA = log10(Matrix::colSums(rna)), 
    nFeature_RNA = Matrix::colSums(rna > 0),
    log10_nCount_ADT = log10(Matrix::colSums(adt)), 
    percent.mt = (Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna))*100
  )
  
  md$drop.class <-  ifelse(rownames(md) %in% called_cells, 'cell', 'background')
  
  # Filter barcodes with no counts
  md <- md %>% filter (log10_nCount_RNA > 0 & log10_nCount_ADT > 0)
  
  # Set limits to define background ADT matrix
  # ADT cutoff
  if(sample.id %in% paste0("aim_", c("ID142_v02",
                                     "ID142_v11b",
                                     "ID142_v11c"))){
    
    bkg.adt.min <- 1.5
    bkg.adt.max <- 3.5
    
  } else {
    
    bkg.adt.min <- 1.5
    bkg.adt.max <- 3
    
  }
  
  bkg.adt.rna <- 2.5
  
  # Visualize
  bkg.adt.plot <- ggplot(md,
                         aes(x = log10_nCount_ADT, 
                             y = log10_nCount_RNA)) +
    theme_bw() + 
    geom_bin2d(bins = 300) + 
    geom_vline(xintercept = bkg.adt.min, color = "red") +
    geom_vline(xintercept = bkg.adt.max, color = "red") +
    geom_hline(yintercept = bkg.adt.rna, color = "red") +
    scale_fill_viridis_c(option = "C",
                         limits = c(0, 500),
                         oob = scales::squish) + 
    facet_wrap(~drop.class)
  
  ggsave(filename = paste0(processed_data, "steps/adt/dsb/", sample.id, "_adt_background.png"),
         plot = bkg.adt.plot,
         width = 25,
         height = 15,
         units = "cm",
         dpi = 300)
  
  # Create background adt matrix
  background_drops <- md %>% 
    filter (log10_nCount_ADT > bkg.adt.min &
              log10_nCount_ADT < bkg.adt.max & 
              log10_nCount_RNA < bkg.adt.rna) %>% 
    rownames()
  
  background.adt.mtx <-  as.matrix(adt[ , background_drops])
  
  # Filter metadata
  cellmd <- md %>% filter(drop.class == "cell")
  
  # Export csv's for quality check
  
  write.table(cellmd,
              paste0(processed_data, "steps/qc/adt_", sample.id, "_cells_features.csv"),
              quote = F,
              sep = ",")
  
  ## Filter matrix -------------------------------------------------------------
  
  # Define filtering limits
  adt.lower <- median(cellmd$log10_nCount_ADT) - 3*mad(cellmd$log10_nCount_ADT)
  adt.upper <- median(cellmd$log10_nCount_ADT) + 3*mad(cellmd$log10_nCount_ADT)
  rna.lower <- log10(1000)
  gen.lower <- 500
  mt.upper <- 7.5
  
  # Find high quality cells
  qc_cells <- cellmd %>% 
    filter(log10_nCount_ADT > adt.lower & 
             log10_nCount_ADT < adt.upper &
             log10_nCount_RNA > rna.lower &
             cellmd$nFeature_RNA > gen.lower & 
             cellmd$percent.mt < mt.upper) %>%
    rownames()
  
  # Create raw matrices with high quality cells alone
  cell.adt.raw <-  as.matrix(adt[ , qc_cells])
  
  cell.rna.raw <-  rna[ ,qc_cells]
  
  cellmd <-  cellmd[qc_cells, ]
  
  ## Create output -------------------------------------------------------------
  
  out <- list(cell.rna.raw, cell.adt.raw, background.adt.mtx, cellmd)
  
  names(out) <- c("rna", "pos_adt", "neg_adt", "meta")
  
  return(out)
  
}

# Create raw (target cells + impurities) object
matrices <- lapply(samples, matrix.cleaner)

names(matrices) <- samples

## Create DSB assay ############################################################

# Identify Isotype controls
isotype.controls <- grep("Isotype", rownames(matrices[[1]]$pos_adt), value = T)

# Create DSB matrices
dsb.matrix.list <- lapply(matrices, function(m){
  
  # Positive matrix to clean
  adt_cells <- m$pos_adt
  
  # Negative matrix to calculate noise
  adt_bg <- m$neg_adt
  
  # Run DSB
  dsb.mat <- DSBNormalizeProtein(cell_protein_matrix = adt_cells, 
                                 empty_drop_matrix = adt_bg, 
                                 denoise.counts = TRUE, 
                                 use.isotype.control = TRUE, 
                                 isotype.control.name.vec = isotype.controls)
  
})

names(dsb.matrix.list) <- samples

# Collect DSB data matrices into one
dsb_common <- do.call(cbind, dsb.matrix.list)

# Keep only cells that are in the final object
pre_dsb <- dsb_common[,intersect(colnames(dsb_common), colnames(AIM))]

# Add DSB assay to AIM object
AIM[["DSB"]] <-  CreateAssay5Object(data = pre_dsb)

## Create ADT (CLR-normalized) assay ###########################################

# Create one assay for each matrix
adt.seur.list <- lapply(matrices, function(x){
  
  m <- x$pos_adt
  
  # Create assay while removing isotypes (keep off CLR normalization)
  CreateAssay5Object(counts = m[!rownames(m) %in% c(isotype.controls), ])
  
})

# Merge into one assay obect
merged.adt <- merge(adt.seur.list[[1]],
                    unlist(adt.seur.list[2:length(adt.seur.list)]),
                    merge.data = TRUE)

# Keep only cells that are in the final object
common.adt <- subset(merged.adt, cells = intersect(colnames(merged.adt),
                                                   colnames(AIM)))

# CLR normalize individually
common.adt <- NormalizeData(common.adt, 
                            normalization.method = "CLR", 
                            margin = 2)

# Join layers into one single assay
common.adt <- JoinLayers(common.adt)

# Add assay to AIM
AIM[["ADT"]] <-  CreateAssay5Object(data = common.adt$data)

## ADT differential expression #################################################

#### Markers pre-filtering -----------------------------------------------------

# Keep only features with mean expression > mean expression of specific isotype
# with Z > 2 - following Wei et al., Immunity 2025

# Read isotype data in
totalseq <- read.csv2("additional_data/totalseq_isotypes.txt", sep="") %>% 
  tibble()

# Extract assay data
adt.mat <- GetAssayData(AIM, 
                        assay = "DSB",
                        layer = "data") %>%
  data.frame()

# Summarize assay data across fatures
adt.mat$mean <- rowMeans(adt.mat, na.rm=TRUE)
adt.mat$sd <- apply(adt.mat, 1, sd, na.rm=TRUE)
adt.mat$marker <- rownames(adt.mat)

# Remove cell barcodes
adt.mat <- adt.mat[c("marker", "mean", "sd")] %>% tibble()

# Keep isotypes from totalseq data sheet
iso <- totalseq %>% 
  mutate(marker = feature) %>% 
  select(marker, Isotype)

# Attach isotype mean and sd
iso.master <- full_join(adt.mat, iso) %>% 
  filter(grepl("Isotype", Isotype)) %>% 
  mutate(Isotype = sub("Isotype ", "", Isotype)) %>% 
  rename(iso.mean = mean,
         iso.sd = sd) %>% 
  select(Isotype, iso.mean, iso.sd)

# Attach isotype information to markers
adt.master <- full_join(adt.mat, iso) %>% 
  filter(!grepl("Isotype", Isotype))

# Join markers and isotype means and sd
final.tab <- full_join(adt.master, iso.master)

# Calculate z-score
final.tab$x2.sample.z.score <- (final.tab$mean - final.tab$iso.mean)/sqrt(final.tab$sd^2/ncol(AIM[["DSB"]]) +
                                                                            final.tab$iso.sd^2/ncol(AIM[["DSB"]]))

# Pull markers with z-score > 2 
to.test <- final.tab %>% 
  filter(x2.sample.z.score > 2) %>% 
  pull(marker)

### Markers testing ------------------------------------------------------------

# Wilcoxon for significance
adt.markers <- FindAllMarkers(AIM,
                              assay = "ADT",
                              features = to.test,
                              logfc.threshold = 0.25,
                              only.pos = T,
                              test.use = "wilcox")

# ROC for characterizing markers
adt.roc <- FindAllMarkers(AIM,
                          assay = "ADT",
                          features = to.test,
                          logfc.threshold = 0.25,
                          only.pos = T,
                          test.use = "roc")

# Pool informations: keep markers with AUC > 0.7 and adjusted p < 0.01
sig.markers <- adt.roc[c("avg_log2FC", "cluster", "gene", "myAUC")] %>% 
  left_join(adt.markers[c("cluster", "gene", "p_val_adj")]) %>% 
  filter(p_val_adj < 0.01) %>% 
  pull(gene)

## Significant markers expression clustering and scaling #######################

### Prepare basis data ---------------------------------------------------------

# Extract cluster info
cluster.info <- AIM@meta.data["w_clusters"] %>% 
  rownames_to_column(var = "bc")

# Get DSB data
adt.data <- GetAssayData(AIM,
                         assay = "DSB",
                         layer = "data") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "bc") %>%
  pivot_longer(!bc,
               names_to = "marker",
               values_to = "expression") %>% 
  left_join(cluster.info, .)

# Define positive cells and summarize data
plot.data <- adt.data %>% 
  # Keep only relevant features
  filter(marker %in% sig.markers) %>% 
  # Define threshold value to DSB = 4
  mutate(thr_expression = ifelse(expression > 4,
                                 expression,
                                 NA)) %>% 
  # Summarize by cluster for each marker
  group_by(marker, w_clusters) %>% 
  summarise(avg_expr = mean(expression),
            n = n(),
            pct_expr = sum(is.na(thr_expression) == F)/n,
            .groups = "keep")

### Cluster features and clusters based on expression values -------------------

# Transform data to matrix
avg.mat <- plot.data %>% 
  select(marker, w_clusters, avg_expr) %>%
  pivot_wider(names_from = w_clusters,
              values_from = avg_expr) %>% 
  column_to_rownames(var = "marker") %>% 
  as.matrix()

# Calculate distance matrix
marker.dist <- dist(avg.mat)

# Cluster markers
marker.clust <- hclust(marker.dist)

# Pull marker order
marker.order <- rownames(avg.mat[marker.clust$order,])

# Calculate distance for w_clusters (transpose matrix, since dist works on rows)
cluster.dist <- dist(t(avg.mat))

# Cluster w_clusters
cluster.clust <- hclust(cluster.dist)

# Pull w_clusters order
cluster.order <- colnames(avg.mat[,cluster.clust$order])

# Assign factor levels for clusters
plot.data$w_clusters <- factor(as.character(plot.data$w_clusters),
                               levels = cluster.order)

### Scale expression by feature ------------------------------------------------

scaled.data <- plot.data %>% 
  # Transform to matrix with markers as columns
  select(marker, w_clusters, avg_expr) %>% 
  pivot_wider(names_from = marker,
              values_from = avg_expr) %>% 
  column_to_rownames(var = "w_clusters") %>% 
  as.matrix() %>% 
  # Scale by column
  scale() %>% 
  # Back to original format
  as.data.frame() %>% 
  rownames_to_column(var = "w_clusters") %>% 
  pivot_longer(!w_clusters,
               values_to = "scaled_expr",
               names_to = "marker")

# Define factor levels for markers (do for both scaled.data and plot.data)
scaled.data$marker <- factor(scaled.data$marker,
                             levels = marker.order)

plot.data$marker <- factor(plot.data$marker,
                           levels = marker.order)

### Collect in significant markers df ------------------------------------------

adt.markers.plot.tab <- plot.data %>% 
  left_join(scaled.data)

adt.markers.plot.tab$w_clusters <- factor(as.character(plot.data$w_clusters),
                                          levels = cluster.order)

## Save outputs ################################################################

# Save Seurat object with ADT assay
saveRDS(AIM, paste0(processed_data, "steps/objects/aim_annotated_adt.rds"))

# Save rds objects
adt_results <- list(adt_z_score_selection = final.tab,
                    adt_test = list(adt_wilcox = adt.markers,
                                    adt_roc = adt.roc),
                    adt_selected = list(adt_plot_tab = adt.markers.plot.tab,
                                        adt_marker_clust = cluster.clust))

saveRDS(adt_results, paste0(processed_data, "steps/adt/adt_results.rds"))

# Export test results for supplementary
write.table(adt.markers,
            file = paste0(processed_data, "results/adt/adt_markers.txt"),
            col.names = T,
            row.names = F)

write.table(adt.roc,
            file = paste0(processed_data, "results/adt/adt_markers_roc.txt"),
            col.names = T,
            row.names = F)
