################################################################################
## Title        : Gene expression analysis
##
## Input        : Annotated AIM object with RNA and ADT (not used)
##
## Output       : 1. Differential gene testing results 
##                2. Gene ontology enrichment results
##                3. Antiviral module scoring results
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading required libraries:##################################################

library(clusterProfiler)
library(enrichplot)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(Seurat)
library(tidyverse)

## Environment setup:###########################################################

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

AIM <- readRDS(paste0(processed_data, "steps/objects/aim_annotated_adt.rds"))

Idents(AIM) <- "w_clusters"

## Differential expression testing #############################################

### Cluster 1 vs 10 ------------------------------------------------------------

# Testing with MAST, adjusting for sample origin
cl1v10 <- FindMarkers(AIM,
                      ident.1 = 1,
                      ident.2 = 10,
                      logfc.threshold = 0.5,
                      min.pct = 0.1,
                      only.pos = T,
                      test.use = "MAST",
                      latent.vars = c("donor"))

# Extract all upregulated genes
cl1v10.genes <- cl1v10 %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>%
  rownames()

# Top 25 genes
top25.10 <- cl1v10 %>% 
  filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>% 
  top_n(n = 25, wt = avg_log2FC) %>% 
  rownames()

### Cluster 1 vs 3 -------------------------------------------------------------

# Testing with MAST, adjusting for sample origin
cl1v3 <- FindMarkers(AIM,
                     ident.1 = 1,
                     ident.2 = 3,
                     logfc.threshold = 0.5,
                     min.pct = 0.1,
                     only.pos = T,
                     test.use = "MAST",
                     latent.vars = "donor")

# Extract all upregulated genes
cl1v3.genes <- cl1v3 %>%
  filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>%
  rownames()

# Top 25 genes
top25.3 <- cl1v3 %>% 
  filter(p_val_adj < 0.01 & avg_log2FC > 0.5) %>% 
  top_n(n = 25, wt = avg_log2FC) %>% 
  rownames()

### Top 25 + 25 expression heatmap build ---------------------------------------

# Average expression per cluster 
avg.aim <- AverageExpression(subset(AIM, idents = c(1,3,10)),
                             assays = "RNA",
                             features = union(top25.10,
                                              top25.3),
                             return.seurat = T)

# Extract data matrix
rna.mat <- GetAssayData(avg.aim,
                        assay = "RNA",
                        layer = "scale.data")

# Split into v3 genes and v10 genes
v3.mat <- rna.mat[top25.3,]
v10.mat <- rna.mat[top25.10,]

# Cluster genes by expression
order.finder <- function(m){
  
  # Calculate distance
  mat_dist <- dist(m, method = "euclidean")
  
  # Cluster
  mat_clust <- hclust(mat_dist, method = "complete")
  
  # Return gene order
  mat_clust$order
  
}

# Collect into one result table
top25.data <- as.data.frame(v3.mat[order.finder(v3.mat), ]) %>%
  # Bind the two data frames by row
  rbind(as.data.frame(v10.mat[order.finder(v10.mat), ])) %>% 
  rownames_to_column(var = "gene") %>% 
  tibble() %>% 
  # Factorize genes to define order
  mutate(gene = factor(gene,
                       levels = gene)) %>%
  # one column for each value to plot
  pivot_longer(!gene,
               names_to = "cluster",
               values_to = "expression") %>% 
  # Clean cluster name to be shown
  mutate(w_cluster = factor(sub("g", "", cluster),
                            levels = c(3, 1, 10)))

## Gene ontology enrichment analysis ###########################################

### Cluster 1 vs 10 ------------------------------------------------------------

# Run enrich GO on genes enriched in 1 vs 10
cl1v10.go <- enrichGO(cl1v10.genes,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "SYMBOL",
                      ont = "BP",
                      pvalueCutoff = 0.01,
                      qvalueCutoff = 0.05,
                      minGSSize = 15,
                      maxGSSize = 500,
                      universe = rownames(AIM[["RNA"]]))

# Simplify results
cl1v10.go <- pairwise_termsim(cl1v10.go)
cl1v10.simply <- clusterProfiler::simplify(cl1v10.go)

# Collect top 5
cl1v10.top5.tab <- cl1v10.simply@result %>% 
  top_n(n = 5, wt = -log10(p.adjust)) %>% 
  select(ID, Description, p.adjust, Count, FoldEnrichment) %>% 
  mutate(upregulated = "cl10") %>% 
  arrange(desc(FoldEnrichment)) %>% 
  tibble()

### Cluster 1 vs 3 -------------------------------------------------------------

# Run enrich GO on genes enriched in 1 vs 3
cl1v3.go <- enrichGO(cl1v3.genes,
                     OrgDb = "org.Hs.eg.db",
                     keyType = "SYMBOL",
                     ont = "BP",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05,
                     minGSSize = 15,
                     maxGSSize = 500,
                     universe = rownames(AIM[["RNA"]]))

# Simplify results
cl1v3.go <- pairwise_termsim(cl1v3.go)
cl1v3.simply <- clusterProfiler::simplify(cl1v3.go)

# Collect top 5
cl1v3.top5.tab <- cl1v3.simply@result %>% 
  top_n(n = 5, wt = -log10(p.adjust)) %>%  
  select(ID, Description, p.adjust, Count, FoldEnrichment) %>% 
  mutate(upregulated = "cl3") %>% 
  arrange(desc(FoldEnrichment)) %>% 
  tibble()

## Antiviral module score ######################################################

### Calculate module score -----------------------------------------------------

# Define module (from Fuchs et al, Front Immunol, 2019 + Collins et al., 
# Sci Immunology, 2023)
cd8_specific <- list(c("ACTB",
                       "ACTG1",
                       "CD82",
                       "CRTAM",
                       "CTNNA1",
                       "EGR2",
                       "GNG4",
                       "GZMB",
                       "HSP90AB1",
                       "HSPA8",
                       "IFNG",
                       "IL2RA",
                       "MIR155HG",
                       "NR4A2",
                       "PKM",
                       "TAGAP",
                       "TNFRSF18",
                       "TNFRSF9", 
                       "XCL1",
                       "XCL2"))

# Calculate score
AIM <- AddModuleScore(AIM,
                      features = cd8_specific,
                      assay = "RNA",
                      name = "antigen_specific_cd8",
                      seed = 123)

# Extract module scores (save as result now, add again later, avoid saving too
# many Seurat object, preserve modularity)

antiviral_module <- AIM$antigen_specific_cd81

### Antiviral module heatmap build ---------------------------------------------

# Define cluster plot order based on median module expression
cluster.order <- AIM@meta.data %>% 
  # Keep only relevant columns
  select(antigen_specific_cd81, w_clusters) %>% 
  group_by(w_clusters) %>% 
  summarise(median_score = median(antigen_specific_cd81)) %>% 
  # Arrange by score
  arrange(desc(median_score)) %>% 
  pull(w_clusters) %>% 
  as.character()

# Create object with average expression
avg.cd8.antigen <- AverageExpression(AIM,
                                     assays = "RNA",
                                     features = cd8_specific[[1]],
                                     group.by = "w_clusters",
                                     return.seurat = T)

# Extract data matrix
cd8.antigen.matrix <- GetAssayData(avg.cd8.antigen,
                                   assay = "RNA",
                                   layer = "scale.data")

# Use order finder function to cluster genes as defined above
cd8.antigen.order <- order.finder(cd8.antigen.matrix)

# Create tibble for plot
plot.cd8.agen.clust <- cd8.antigen.matrix[cd8.antigen.order, ] %>%  
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "gene") %>% 
  tibble() %>%
  mutate(gene = factor(gene, levels = gene)) %>%  
  pivot_longer(!gene,
               names_to = "w_cluster",
               values_to = "expr") %>% 
  # Remove prefix g from cluster names
  mutate(w_cluster = sub("g", "", w_cluster)) %>% 
  # Factorize to ensure correct clustered order
  mutate(w_cluster = factor(w_cluster, 
                            levels = cluster.order))

# Create tibble with median scores
antigen.median.scores <- AIM@meta.data %>% 
  select(antigen_specific_cd81, w_clusters) %>% 
  group_by(w_clusters) %>% 
  summarise(median_score = median(antigen_specific_cd81))

## Collect and export results ##################################################

# To be added to Seurat object later
gene_expr <- list(cluster1_markers = list(vs3_all = cl1v3,
                                          vs10_all = cl1v10,
                                          top25 = top25.data),
                  cluster1_go = list(vs3_go = cl1v3.simply,
                                     vs3_top5 = cl1v3.top5.tab,
                                     vs10_go = cl1v10.simply,
                                     vs10_top5 = cl1v10.top5.tab),
                  antiviral_cd8_module = list(module_score = antiviral_module,
                                              heatmap_data = plot.cd8.agen.clust,
                                              median_scores = antigen.median.scores))

saveRDS(gene_expr, paste0(processed_data, "steps/deg/gene_expression.rds"))

# For supplementary
write.table(rownames_to_column(cl1v10, var = "gene"),
            file = paste0(processed_data, "results/deg/cluster_1v10.txt"),
            dec = ",",
            col.names = T,
            row.names = F)

write.table(rownames_to_column(cl1v3, var = "gene"),
            file = paste0(processed_data, "results/deg/cluster_1v3.txt"),
            dec = ",",
            col.names = T,
            row.names = F)

write.table(cl1v10.simply,
            file = paste0(processed_data, "results/deg/cluster_1v10_go.txt"),
            col.names = T,
            row.names = F)

write.table(cl1v3.simply,
            file = paste0(processed_data, "results/deg/cluster_1v3_go.txt"),
            col.names = T,
            row.names = F)
