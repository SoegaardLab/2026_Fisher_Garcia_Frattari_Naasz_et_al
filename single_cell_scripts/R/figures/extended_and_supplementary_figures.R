################################################################################
## Title        : Extended Data and Supplementary Figures
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
## Extended Data Figure 5 ######################################################
################################################################################

### Extended Data Figure 5A: Marker Genes ######################################

AIM@meta.data <- AIM@meta.data %>% 
  mutate(dot_clusters = factor(as.character(w_clusters),
                               levels = c(11:1)))
efig.5a <- DotPlot(AIM, 
                   assay = "RNA",
                   features = rna.feats,
                   group.by = "dot_clusters",
                   dot.scale = 2) +
  theme_classic() +
  basic_theme +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1),
        axis.title.x = element_blank(),
        axis.line = element_line(linewidth = 0.5),
        legend.title = element_text(vjust = 0.5),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(rep(5, 4))) +
  scale_color_gradientn(colors = c("white", rev(heat.colors(100)))) +
  guides(color = guide_colorbar(title.position = "left",
                                barheight = 0.4)) +
  labs(color = "Average Expression",
       size = "Percent expressed",
       y = "Cluster")

efig.5a

# Export data
write.table(efig.5a$data,
            file = paste0(processed_data, "results/figures/data_efigure5a.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_5a.svg"),
       plot = efig.5a,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_5a.png"),
       plot = efig.5a,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

### Extended Data Figure 5B: Marker proteins ###################################

cluster.order <- as.character(AIM@misc$adt$adt_selected$adt_marker_clust$order)

adt.dots <- AIM@misc$adt$adt_selected$adt_plot_tab %>% 
  ggplot(aes(y = w_clusters, 
             x = marker)) +
  geom_point(aes(size = pct_expr *100,
                 color = scaled_expr)) +
  theme_classic() +
  scale_color_gradientn(colors = c("white", rev(heat.colors(100)))) +
  scale_y_discrete(limits = cluster.order) +
  scale_size_continuous(range = c(0, 2)) +
  basic_theme +
  labs(color = "Average Expression", size = "Percent Expressed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.title = element_blank(),
        legend.title = element_text(vjust = 0.5),
        legend.position = "bottom",
        axis.line = element_line(linewidth = 0.5)) +
  guides(color = guide_colorbar(title.position = "left",
                                barheight = 0.4,
                                order = 1),
         size = guide_legend(order = 2))

dend <- as.dendrogram(AIM@misc$adt$adt_selected$adt_marker_clust)

dend_data <- dendro_data(dend)

dend_plot <- ggdendrogram(dend,
                          rotate = TRUE,      # rotate so it's vertical
                          theme_dendro = FALSE) +
  scale_y_reverse() + # Flip left to right
  scale_x_discrete(limits = cluster.order) +
  theme_void() +
  theme(plot.margin = margin(r = 0))

efig.5b <- dend_plot + adt.dots + plot_layout(widths = c(0.075, 1))

efig.5b

# Export data
write.table(efig.5b$data,
            file = paste0(processed_data, "results/figures/data_efigure5b.txt"),
            col.names = T,
            row.names = F)

write.table(dend_data$segments,
            file = paste0(processed_data, "results/figures/data_efigure5b_dendrogram.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_5b.svg"),
       plot = efig.5b,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_5b.png"),
       plot = efig.5b,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

### Extended Data Figure 5C: ADT expression ####################################

adt.expr <- GetAssayData(AIM,
                         assay = "DSB")[adt.feats,] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "barcode") %>% 
  left_join(cluster.info) %>% 
  tibble() %>% 
  pivot_longer(all_of(adt.feats),
               names_to = "feature",
               values_to = "expression") %>% 
  mutate(expr_color = factor(ifelse(expression > 4, "Positive", "Negative"),
                             levels = c("Positive", "Negative")))

adt.rank <- adt.expr %>% 
  group_by(feature, w_clusters) %>% 
  summarise(avg_expr = mean(expression)) %>% 
  group_by(feature) %>% 
  arrange(desc(avg_expr), .by_group = T) %>% 
  mutate(feat_rank = row_number()) %>% 
  select(feature, w_clusters, feat_rank)

marker.plotter <- function(marker){
  
  adt.expr %>% 
    filter(feature == marker) %>% 
    left_join(adt.rank) %>% 
    arrange(feat_rank) %>% 
    mutate(w_clusters = factor(w_clusters,
                               levels = unique(w_clusters))) %>% 
    ggplot(aes(x = w_clusters,
               y = expression)) +
    ggbeeswarm::geom_quasirandom(size = AutoPointSize(AIM@meta.data),
                                 aes(color = expr_color)) +
    geom_boxplot(aes(fill = w_clusters),
                 width = 0.4,
                 outliers = F,
                 show.legend = F) +
    geom_hline(yintercept = 4,
               linetype = "dashed") +
    scale_fill_manual(values = w_cluster_color) +
    scale_color_manual(values = c("red", "gray"),
                       breaks = c("Positive", "Negative")) +
    labs(title = marker,
         color = "Positivity cutoff\n(Expression > 4)",
         x = "Clusters\n(ranked by average expression)",
         y = "DSB-normalized ADT expression")
  
}

adt.violins <- lapply(adt.feats, marker.plotter)

efig.5c <- wrap_plots(
  adt.violins, 
  nrow = 1, 
  guides = "collect", 
  axes = "collect") &
  theme_classic() &
  basic_theme &
  theme(plot.margin = margin(rep(5, 4)),
        legend.position = "bottom")

efig.5c

efig.5c.data <- efig.5c %>%
  lapply(., function(x) x$data) %>%
  do.call(rbind, .)

# Export data
write.table(efig.5c.data,
            file = paste0(processed_data, "results/figures/data_efigure5c.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_5c.svg"),
       plot = efig.5c,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_5c.png"),
       plot = efig.5c,
       width = 20.6,
       height = 6.8,
       units = "cm",
       dpi = 300)

### Extended Data Figure 5D: DAB results #######################################

pre.ati.clusters.tab <- AIM@misc$dab$cluster_proportions_art

ati.clusters.tab <- AIM@misc$dab$cluster_proportions_ati

ylim <- range(c(pre.ati.clusters.tab$freq, ati.clusters.tab$freq))

avg.pre.ati <- pre.ati.clusters.tab %>% 
  group_by(cluster_label, status) %>% 
  summarise(freq = mean(freq)) %>% 
  mutate(min = ifelse(status == "PICs", 0.75, 1.75),
         max = ifelse(status == "PICs", 1.25, 2.25)) %>%  
  mutate(status = factor(status,
                         levels = c("PICs", "NCs")))

avg.ati <- ati.clusters.tab %>% 
  group_by(cluster_label, status) %>% 
  summarise(freq = mean(freq)) %>% 
  mutate(min = ifelse(status == "PICs", 0.75, 1.75),
         max = ifelse(status == "PICs", 1.25, 2.25)) %>% 
  mutate(status = factor(status,
                         levels = c("PICs", "NCs")))

pre.ati.clusters <- pre.ati.clusters.tab %>%
  ggplot(aes(x = status,
             y = freq)) +
  ggbeeswarm::geom_beeswarm(aes(color = donor),
                            cex = 5,
                            size = 1.5,
                            shape = 1,
                            stroke = 1) + 
  geom_linerange(data = avg.pre.ati, 
                 aes(xmin = min, xmax = max)) +
  annotate(geom = 'segment',
           y = Inf,
           yend = Inf,
           x = 0.5,
           xend = 2.5) +
  scale_color_manual(values = id.cols) +
  theme_minimal() +
  facet_grid(clinical_time_point ~ cluster_label) +
  theme(legend.position = "none") +
  coord_cartesian(ylim = ylim) +
  labs(y = "Proportion of cluster in AIM-sorted cells",
       color = "Participant ID")

ati.clusters <- ati.clusters.tab %>% 
  ggplot(aes(x = status,
             y = freq)) +
  ggbeeswarm::geom_beeswarm(aes(color = donor),
                            cex = 5,
                            size = 1.5,
                            shape = 1,
                            stroke = 1) +
  geom_linerange(data = avg.ati, 
                 aes(xmin = min, xmax = max)) +
  annotate(geom = 'segment',
           y = Inf,
           yend = Inf,
           x = 0.5,
           xend = 2.5) +
  scale_color_manual(values = id.cols) +
  theme_minimal() +
  facet_grid(clinical_time_point ~ cluster_label) +
  coord_cartesian(ylim = ylim) +
  theme(legend.direction = "horizontal",
        legend.position = "bottom") +
  labs(y = "Proportion of cluster in AIM-sorted cells",
       color = "Participant ID")

efig.5d <- pre.ati.clusters + ati.clusters +
  plot_layout(ncol = 1,
              axes = "collect") +
  plot_annotation(title = "Comparison of Cluster Abundances Between Groups") &
  theme(legend.margin = margin(b = 0, t = 0.1),
        strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y = element_part_rect(side = "l"),
        strip.text = element_text(face = "bold", size = 8),
        legend.title = element_text(size = 6, face = "bold"),
        plot.title = element_text(size = 7,
                                  face = "bold",
                                  margin = margin()),
        axis.line.y = element_line(linewidth = 0.25),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(face = "bold",
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 6),
        axis.title.y = element_text(size = 6, face = "bold"),
        axis.text.y = element_text(size = 5),
        legend.box.margin = unit(2, "mm"),
        plot.margin = margin(l = 5, r = 5))

efig.5d

efig.5d.data <- efig.5d %>%
  lapply(., function(x) x$data) %>%
  do.call(rbind, .)

efig.5d.avg.data <- rbind(avg.pre.ati, avg.ati)

# Export data
write.table(efig.5d.data,
            file = paste0(processed_data, "results/figures/data_efigure5d.txt"),
            col.names = T,
            row.names = F)

write.table(efig.5d.avg.data,
            file = paste0(processed_data, "results/figures/data_efigure5d_averaged.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_5d.svg"),
       plot = efig.5d,
       width = 20,
       height = 8,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_5d.png"),
       plot = efig.5d,
       width = 20,
       height = 8,
       units = "cm",
       dpi = 300)

################################################################################
## Extended Data Figure 6 ######################################################
################################################################################

### Extended Data Figure 6A: Clonal ranks ######################################

# Pseudonymize TCR data
pseudo_tcr_w_ranks <- AIM@misc$tcr$tcr_w_ranks %>%
  # Keep one entry per clonotype per donor
  distinct(donor, CT) %>%
  # Create pseudonym
  mutate(plot.ct = paste("ct", row_number(), sep = "_")) %>%
  # Rejoin to original table
  right_join(AIM@misc$tcr$tcr_w_ranks) %>%
  # Swap original CT with pseudonyms
  select(!CT) %>%
  rename(CT = plot.ct) %>%
  ungroup() %>%
  # Keep minimal set
  select(barcode, donor, visit, CT, prevalence, rank_group, w_clusters, control)

efig.6a <- pseudo_tcr_w_ranks %>% 
  mutate(w_clusters = factor(w_clusters,
                             levels = c(11:1))) %>% 
  ggplot(aes(y = w_clusters)) +
  geom_bar(aes(fill = forcats::fct_rev(rank_group)),
           position = "fill",
           color = "gray30",
           alpha = 0.8,
           width = 0.75,
           linewidth = 0.25) +
  geom_vline(xintercept = 0.5, 
             linetype = "dashed",
             linewidth = 0.35) +
  facet_nested_wrap(~control + donor, 
                    nest_line = element_line(linewidth = 0.5),
                    nrow = 1,
                    strip = strip_nested(text_x = elem_list_text(face = c("bold",
                                                                          "plain")),
                                         by_layer_x = T)) +
  scale_fill_manual(values = rank.colors) +
  scale_x_continuous(labels = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 6),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(linewidth = 0.25),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.title = element_text(size = 7, face = "bold"),
        axis.text = element_text(size = 5),
        axis.title = element_text(size = 6,
                                  face = "bold"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6, face = "bold"),
        legend.margin = margin(b = 0, t = 0.1)) +
  labs(fill = "Within-participant\nclonotype indices",
       y = "Cluster",
       x = "Proportion of rank among T cells with productive TCR-\u03b1 and TCR-\u03b2",
       title = "Clonal Rank Distribution by Cluster and Participant") +
  guides(fill = guide_legend(reverse = T,
                             title.position = "left"))

efig.6a

# Export data
write.table(efig.6a$data,
            file = paste0(processed_data, "results/figures/data_efigure6a.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_6a.svg"),
       plot = efig.6a,
       width = 20.2,
       height = 7.6,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_6a.png"),
       plot = efig.6a,
       width = 20.2,
       height = 7.6,
       units = "cm",
       dpi = 300)

### Extended Data Figure 6B: Clonality indices #################################

# Plot
efig.6b <- AIM@misc$tcr$tcr_clonality %>% 
  ggplot(aes(x = control,
             y = gini)) +
  ggbeeswarm::geom_beeswarm(aes(color = donor),
                            cex = 5,
                            size = 1.5,
                            shape = 1,
                            stroke = 1) +
  annotate(geom = 'segment',
           y = Inf,
           yend = Inf,
           x = 0.5,
           xend = 2.5) +
  scale_color_manual(values = id.cols) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line()) +
  facet_grid(~w_clusters) +
  theme(strip.placement = "outside",
        strip.background.x = element_blank(),
        strip.background.y = element_rect(),
        strip.text = element_text(face = "bold", size = 6),
        axis.text.x = element_text(face = "bold",
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 6),
        axis.text.y = element_text(size = 5),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom") +
  labs(y = "Gini coefficient",
       color = "Participant ID",
       title = "Median Cluster Clonality Across Visits") +
  theme(legend.margin = margin(b = 0, t = 0.1),
        plot.title = element_text(size = 7,
                                  face = "bold"),
        legend.title = element_text(size = 6, face = "bold"),
        axis.line.y = element_line(linewidth = 0.25),
        axis.ticks.y = element_line(linewidth = 0.25),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 6, face = "bold"),
        axis.line.x = element_line(linewidth = 0.2),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(face = "bold",
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 6),
        axis.text.y = element_text(size = 5))

efig.6b

# Export data
write.table(efig.6b$data,
            file = paste0(processed_data, "results/figures/data_efigure6b.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_6b.svg"),
       plot = efig.6b,
       width = 20,
       height = 7.4,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_6b.png"),
       plot = efig.6b,
       width = 20,
       height = 7.4,
       units = "cm",
       dpi = 300)

### Extended Data Figure 6C: Repertoire overlap ################################

# Heatmap
efig.6c <- AIM@misc$tcr$tcr_overlap %>% 
  ggplot(aes(x = factor(annotation2,
                        levels = levels(AIM$w_clusters)),
             y = factor(annotation1,
                        levels = rev(levels(AIM$w_clusters))),
             fill = morisita)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(morisita <= 0.05, "", round(morisita, 2))),
            size.unit = "pt",
            size = 6) +
  scale_fill_gradient(low =  "gray90",
                      high = "red",
                      limits = c(0.2, 1), 
                      na.value = "gray90") +
  theme_minimal() +
  facet_nested_wrap(~control + donor1, 
                    nest_line = element_line(linewidth = 0.5),
                    nrow = 2,
                    strip = strip_nested(text_x = elem_list_text(face = c("bold",
                                                                          "plain")),
                                         by_layer_x = T)) +
  theme(strip.background = element_blank(),
        legend.position = "right",
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 7, 
                                  face = "bold",
                                  margin = margin(t = 0,
                                                  b = 0.5)),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 6,
                                  face = "bold"),
        legend.title = element_text(size = 6, 
                                    face = "bold"),
        legend.text = element_text(size = 6),
        legend.margin = margin(b = 0, t = 0.1)) +
  guides(fill = guide_colorbar(title.position = "top",
                               barwidth = 0.4)) +
  labs(x = "Cluster",
       y = "Cluster",
       fill = "Morisita-Horn\nindex",
       title = "TCR Repertoire Overlap Between Clusters")

efig.6c

# Export data
write.table(efig.6c$data,
            file = paste0(processed_data, "results/figures/data_efigure6c.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/extended_figure_6c.svg"),
       plot = efig.6c,
       width = 20.2,
       height = 8.5,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/extended_figure_6c.png"),
       plot = efig.6c,
       width = 20.2,
       height = 8.5,
       units = "cm",
       dpi = 300)

### Supplementary Figure 5: Visit Overview #####################################

visit.plot.data <- visit %>% 
  mutate(study = ifelse(ID == "ID142", "TITAN", "eCLEAR")) %>% 
  mutate(control = factor(ifelse(ID %in% c("ID112", "ID104"), 
                                 "NCs",
                                 "PICs"),
                          levels = c("PICs",
                                     "NCs")))

sfig.5 <- visit.plot.data %>%
  left_join(rgb.v.colors) %>% 
  filter(original %in% AIM$visit) %>%
  filter(ID %in% AIM$donor) %>% 
  mutate(ID = factor(ID,
                     levels = rev(c("ID107",
                                    "ID142",
                                    "ID104",
                                    "ID112")))) %>% 
  ggplot(aes(x = week,
             y = ID)) +
  geom_rect(xmin = -Inf, xmax = 0,
            ymin = -Inf, ymax = Inf,
            fill = "#EAEAEA") +
  geom_linerange(xmin = -Inf, 
                 xmax = Inf,
                 color = "gray50") +
  geom_point(aes(color = visit_colors),
             size = 2) +
  geom_vline(xintercept = 0,
             linetype = "longdash",
             linewidth = 0.25,
             color = "gray20") +
  scale_x_continuous(breaks = c(-52, 0, 50, 100, 150, 200, 250, 300)) +
  scale_color_manual(values = visit_colors, guide = "none") +
  theme_classic() +
  facet_wrap(~control,
             strip.position = "left",
             ncol = 1,
             scales = "free_y") +
  basic_theme +
  theme(strip.placement = "outside",
        strip.background = element_part_rect("r"),
        strip.text = element_text(face = "bold",
                                  size = 6),
        axis.line = element_line(linewidth = 0.25),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_blank()) +
  labs(x = "ATI week",
       title = "Overview of AIM-Sorted Single-Cell-Sequenced Samples")

sfig.5 

# Export data
write.table(sfig.5$data,
            file = paste0(processed_data, "results/figures/data_sfigure5.txt"),
            col.names = T,
            row.names = F)

ggsave(paste0(processed_data, "results/figures/supplementary_figure_5.svg"),
       plot = sfig.5 ,
       width = 10.4,
       height = 6.2,
       units = "cm",
       dpi = 300)

ggsave(paste0(processed_data, "results/figures/supplementary_figure_5.png"),
       plot = sfig.5,
       width = 10.4,
       height = 6.2,
       units = "cm",
       dpi = 300)

print("--- Done! ---")
