################################################################################
## Title        : Figures
##
## Input        : Annotated AIM+ object with all data collected
##
## Output       : AIM+ sorted, single-cell results figures
##
## Author       : Giacomo S Frattari
################################################################################
##
## Environment setup:###########################################################

print("--- Setting up the environment ---")

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

# Figure utilities
source("R/utilities/figure_utilities.R")

### AIM object -----------------------------------------------------------------

print("--- Loading data ---")

AIM <- readRDS(paste0(processed_data, "results/objects/aim_final.rds"))

# Make sure that working clusters are the identities
Idents(AIM) <- "w_clusters"

# Make cluster labels to a separate data frame
cluster.info <- AIM@meta.data[c("barcode", "w_clusters", "donor")]

print("--- Generating figure panels ---")

################################################################################
## Figure 5 ####################################################################
################################################################################

### Figure 5A: AIM+ UMAP #######################################################

fig.5a <- DimPlot(
  AIM,
  group.by = "w_clusters",
  label = T,
  label.size = 2.4,
  raster = T,
  raster.dpi = c(250, 250)) +
  scale_color_manual(values = w_cluster_color, 
                     labels = levels(AIM$annotation)) +
  theme_minimal() +
  basic_theme +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.height = unit(1, "lines")) +
  guides(color = guide_legend(override.aes = list(size = 1.75))) +
  labs(x = "UMAP 1",
       y = "UMAP 2",
       title = "AIM-Sorted T-Cells")

fig.5a

# Export data
write.table(rownames_to_column(fig.5a$data, var = "barcode"),
            file = paste0(processed_data, "results/figures/data_figure5a.txt"),
            col.names = T,
            row.names = F)

# Save figure
ggsave(paste0(processed_data, "results/figures/figure_5a.svg"),
       plot = fig.5a,
       width = 10.5,
       height = 6.5,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5a.png"),
       plot = fig.5a,
       width = 10.5,
       height = 6.5,
       units = "cm",
       dpi = 300)

### Figure 5B: Virus-Specific heatmap ##########################################

#### Create heatmap ------------------------------------------------------------

gene.heat <- AIM@misc$rna$antiviral_cd8_module$heatmap_data %>% 
  ggplot(aes(x = w_cluster,
             y = gene,
             fill = expr)) +
  geom_tile(color = "gray50",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradientn(colors = colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 7,
                                 name = "RdYlBu")))(100)) +
  labs(fill = "Mean scaled\nRNA\nexpression",
       x = "Clusters",
       y = "Module genes") +
  theme_classic () +
  basic_theme +
  theme(axis.text.y = element_text(hjust = 1, 
                                   size = 6),
        axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5,
                                   size = 6,
                                   face = "bold"),
        axis.line = element_blank(),
        axis.ticks = element_blank()) +
  guides(fill = guide_colorbar(barwidth = 0.4))

# Export data
write.table(gene.heat$data,
            file = paste0(processed_data, "results/figures/data_figure5b_heatmap.txt"),
            col.names = T,
            row.names = F)

#### Create top bar ------------------------------------------------------------

median.scores <- AIM@misc$rna$antiviral_cd8_module$median_scores

cluster.order <- levels(AIM@misc$rna$antiviral_cd8_module$heatmap_data$w_cluster)

# Create top bar with median scores
top.bar <- median.scores %>% 
  ggplot(aes(x = factor(w_clusters, 
                        levels = cluster.order),
             y = median_score + abs(min(median_score)))) +
  geom_col(color = "black",
           fill = "gray") +
  theme_void() +
  scale_y_continuous(name = NULL,
                     sec.axis = dup_axis(name = "Median\nmodule\nscore")) +
  coord_cartesian(ylim = c(min(median.scores$median_score),
                           max(median.scores$median_score))) +
  theme(legend.position = "none",
        axis.title.y.right = element_text(angle = 90,
                                          size = 6))

# Export data
write.table(top.bar$data,
            file = paste0(processed_data, "results/figures/data_figure5b_topbar.txt"),
            col.names = T,
            row.names = F)

#### Put together with patchwork and export ------------------------------------

fig.5b <- top.bar + gene.heat + 
  plot_layout(ncol = 1,
              heights = c(.2, 1)) +
  plot_annotation(title = "Virus-Specific CD8 Activation Module") &
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 7,
                                  face = "bold"),
        plot.margin = margin())

fig.5b

ggsave(paste0(processed_data, "results/figures/figure_5b.svg"),
       plot = fig.5b,
       width = 8,
       height = 7.2,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5b.png"),
       plot = fig.5b,
       width = 8,
       height = 7.2,
       units = "cm",
       dpi = 300)

### Figure 5C: Virus-Specific UMAP #############################################

fig.5c <- FeaturePlot(AIM,
                      features = "antigen_specific_cd8",
                      raster = T,
                      raster.dpi = c(250, 250)) +
  scale_color_viridis_c(option = "B", breaks = c(-0.5, 0.75, 2.0)) +
  theme_minimal() +
  basic_theme +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_colorbar(title.position = "left",
                                barheight = 0.4)) +
  labs(color = "Module score",
       title = "Virus-Specific CD8 Activation Module",
       x = "UMAP 1",
       y = "UMAP2")

fig.5c

# Export data
write.table(rownames_to_column(fig.5c$data, var = "barcode"),
            file = paste0(processed_data, "results/figures/data_figure5c.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5c.svg"),
       plot = fig.5c,
       width = 7.5,
       height = 7.5,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5c.png"),
       plot = fig.5c,
       width = 7.5,
       height = 7.5,
       units = "cm",
       dpi = 300)

### Figure 5D: ADT Feature Plot ################################################

DefaultAssay(AIM) <- "DSB"

fig.5d <- FeaturePlot(AIM, 
                      features = c("CD69", "CD107a", "CD25", "CD71"), 
                      min.cutoff = 4, 
                      max.cutoff = 20,
                      order = T) +
  plot_layout(guides = "collect") &
  labs(color = "DSB-normalized\nADT expression") &
  theme_void() &
  basic_theme &
  theme(plot.title = element_text(size = 6, 
                                  hjust = 0.5, 
                                  face = "bold"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) &
  guides(color = guide_colorbar(barheight = 0.4)) &
  scale_color_gradient(low = "white", 
                       high = heat.colors(100)[1], 
                       na.value = "white",
                       breaks = seq(5, 20, 5),
                       limits = c(4, 20),
                       oob = scales::squish) &
  plot_annotation(title = "Surface Protein Profiling",
                  theme = theme(plot.title = element_text(size = 7,
                                                          face = "bold")))

fig.5d

fig.5d.data <- FeaturePlot(AIM, 
                          features = c("CD69", "CD107a", "CD25", "CD71"), 
                          min.cutoff = 4, 
                          max.cutoff = 20,
                          order = T,
                          combine = F) %>%
  lapply(., function(x){
    
    x$data %>% rownames_to_column(var = "barcode")
    
  }) %>%
  reduce(full_join)

# Export data
write.table(fig.5d.data,
            file = paste0(processed_data, "results/figures/data_figure5d.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5d.svg"),
       plot = fig.5d,
       width = 6.8,
       height = 7.2,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5d.png"),
       plot = fig.5d,
       width = 6.8,
       height = 7.2,
       units = "cm",
       dpi = 300)

### Figure 5E: Cluster abundance at ATI start ##################################

fig.5e <- AIM@misc$dab$cluster_proportions_art_group %>% 
  ggplot(aes(x = status)) +
  geom_col(aes(y = freq,
               fill = cluster),
           color = "black",
           position = "stack",
           linewidth = 0.25,
           alpha = 0.6) +
  geom_text(aes(label = lab,
                y = label.position),
            size = 2.4) +
  scale_fill_manual(values = rev(scales::hue_pal()(11))) +
  theme_classic() +
  basic_theme +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 1),
        axis.text.x = element_text(size = 6, face = "bold"),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  labs (x = "",
        y = "Proportion of cluster\nin AIM-sorted T cells",
        fill = "",
        title = "Cluster Abundance\nat ATI Start")

fig.5e

# Export data
write.table(fig.5e$data,
            file = paste0(processed_data, "results/figures/data_figure5e.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5e.svg"),
       plot = fig.5e,
       width = 3.5,
       height = 7,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5e.png"),
       plot = fig.5e,
       width = 3.5,
       height = 7,
       units = "cm",
       dpi = 300)

### Figure 5F: Cluster 1 heatmap ###############################################

fig.5f <- AIM@misc$rna$cluster1_markers$top25 %>% 
  mutate(rank = row_number()) %>% 
  mutate(group = factor(ifelse(rank <76, # Rank 76 because there are 25 genes x 3 clusters entries before
                               "vs Cluster 3", 
                               "vs Cluster 10"),
                        levels = c("vs Cluster 10",
                                   "vs Cluster 3"))) %>%
  ggplot(aes(x = w_cluster,
             y = gene,
             fill = expression)) +
  geom_tile(color = "gray50",
            lwd = 0.5,
            linetype = 1) +
  scale_fill_gradientn(colors = colorRampPalette(
    rev(RColorBrewer::brewer.pal(n = 7,
                                 name = "RdYlBu")))(100)) +
  facet_wrap(~group, scales = "free_y", ncol = 1) +
  labs(fill = "Mean scaled\nRNA expression",
       x = "Clusters",
       title = "Top 25 Upregulated Genes in Cluster 1") +
  theme_classic() +
  basic_theme +
  theme(axis.text.y = element_text(hjust = 1),
        axis.text.x = element_text(face = "bold",
                                   size = 6),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold",
                                  margin = margin(b = 2))) +
  guides(fill = guide_colorbar(barheight = 0.4))

# Export data
write.table(fig.5f$data,
            file = paste0(processed_data, "results/figures/data_figure5f.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5f.svg"),
       plot = fig.5f,
       width = 6.5,
       height = 12,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5f.png"),
       plot = fig.5f,
       width = 6.5,
       height = 12,
       units = "cm",
       dpi = 300)

### Figure 5G: Cluster 1 GO terms ##############################################

fig.5g <- rbind(AIM@misc$rna$cluster1_go$vs10_top5,
                AIM@misc$rna$cluster1_go$vs3_top5) %>% 
  mutate(Description = ifelse(ID == "GO:0002460",
                              "TCR-based immune response*",
                              Description)) %>%
  mutate(Description = factor(Description,
                              levels = rev(Description))) %>% 
  mutate(upregulated = factor(ifelse(upregulated == "cl10",
                                     "vs Cluster 10",
                                     "vs Cluster 3"),
                              levels = c("vs Cluster 10",
                                         "vs Cluster 3"))) %>% 
  ggplot(aes(y = Description,
             fill = upregulated)) +
  geom_col(aes(x = FoldEnrichment),
           show.legend = F) +
  geom_text(aes(x = 0,
                label = Description),
            hjust = 0,
            size = 2) +
  scale_fill_manual(values = alpha(scales::hue_pal()(11)[c(3, 10)], 0.2)) +
  facet_wrap(~upregulated,
             scales = "free_y",
             ncol = 1) +
  labs(x = "Fold Enrichment",
       title = "Top Enriched GO Terms in Cluster 1",
       fill = "") +
  coord_cartesian(xlim = c(-0.2, 6.4)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_classic() +
  basic_theme +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold",
                                  margin = margin(b = 2))) 

fig.5g

# Export data
write.table(fig.5g$data,
            file = paste0(processed_data, "results/figures/data_figure5g.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5g.svg"),
       plot = fig.5g,
       width = 6,
       height = 11.2,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5g.png"),
       plot = fig.5g,
       width = 6,
       height =11.2,
       units = "cm",
       dpi = 300)

### Figure 5H: Top 10 Clonotypes over Time #####################################

# Pseudonymize clonotypes (avoid CDR3 nt)
alluvial.data <- AIM@misc$tcr$tcr_10_development %>%
  # Keep one entry per clonotype per donor
  distinct(donor, CT) %>%
  # Create pseudonym
  mutate(plot.ct = paste("ct", row_number(), sep = "_")) %>%
  # Rejoin to original table
  right_join(AIM@misc$tcr$tcr_10_development) %>%
  # Swap original CT with pseudonyms
  select(!CT) %>%
  rename(CT = plot.ct)

fig.5h <- alluvial.data %>% 
  left_join(v_to_add) %>% 
  mutate(week = factor(as.character(week),
                       levels = sort(unique(v_to_add$week)))) %>%
  ggplot(aes(x = week,
             stratum = CT,
             alluvium = CT,
             fill = rank,
             y = frequency)) +
  geom_rect(xmin = -Inf, xmax = 3.5,
            ymin = -Inf, ymax = Inf,
            fill = "#EAEAEA") +
  geom_alluvium(decreasing = T,
                knot.pos = 0.05,
                alpha = 0.3) +
  geom_stratum (decreasing = T,
                width = 0.5,
                color = "black",
                linewidth = 0.1) +
  geom_vline(xintercept = 3.5,
             linetype = "longdash",
             linewidth = 0.25,
             color = "gray20") +
  theme_classic() +
  scale_fill_manual(values = top.colors) +
  labs(y = "Clonotype proportion in visit repertoire",
       x = "ATI week",
       fill = "Within-participant clone ranking",
       title = "Top 10 Clonotypes Over Time") +
  theme(panel.border = element_rect(color = "black",
                                    fill = "transparent"),
        axis.line = element_line(linewidth = 0.25)) +
  scale_y_continuous(n.breaks = 4) +
  facet_nested_wrap(~control + donor,
                    strip.position = "right",
                    nest_line = element_line(color = "black"),
                    ncol = 1) +
  basic_theme +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 6,
                                  face = "bold",
                                  margin = margin(l = 3, r = 3)),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.size = unit(0.2, "cm")) +
  guides(fill = guide_legend(nrow = 2, 
                             byrow = TRUE,
                             title.position = "top"))

fig.5h

# Export data
write.table(fig.5h$data,
            file = paste0(processed_data, "results/figures/data_figure5h.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/figure_5h.svg"),
       plot = fig.5h,
       width = 6,
       height = 11.5,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/figure_5h.png"),
       plot = fig.5h,
       width = 6,
       height = 11.5,
       units = "cm",
       dpi = 300)

print("--- Done! ---")
