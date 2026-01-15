################################################################################
## Title        : Figures utilities
##
## Description  : Common set of values (mostly colors) used across figures
##
## Author       : Giacomo S Frattari
################################################################################

## Loading required libraries:##################################################

library(ggh4x)
library(ggalluvial)
library(ggdendro)
library(Seurat)
library(tidyverse)
library(patchwork)

### Visit metadata -------------------------------------------------------------
visit <- readRDS("additional_data/visit_names.rds")

v_to_add <- visit %>% 
  mutate(week_label = paste0("Week ", week)) %>% 
  mutate(week_label = case_when(week == -57 ~ "pre-ART",
                                week %in% c(-2, 0) ~ "pre-ATI",
                                TRUE ~ week_label)) %>% 
  arrange(week) %>% 
  mutate(week_label = factor(week_label,
                             levels = unique(week_label))) %>% 
  rename(donor = ID,
         visit = original)

### Cluster colors -------------------------------------------------------------

w_cluster_color <- scales::hue_pal()(11)
names(w_cluster_color) <- c(1:11)

### Top 10 clonotypes colors ---------------------------------------------------

top.colors <- viridis::viridis(n = 10, option = "D", direction = -1)

names(top.colors) <- c(1:10)

### Rank colors ----------------------------------------------------------------

rank.colors <- viridis::viridis(5, 
                                option = "C",
                                direction = -1)

names(rank.colors) <- c("1:10",
                        "11:50",
                        "51:200",
                        "201:500",
                        "> 500")

### Visit colors ---------------------------------------------------------------

# Follow scheme from th rest of the paper
rgb.v.colors <- data.frame(
  week = c(
    -57, 0, 14, 29, 95, 130, 181, 227, 281,
    -2, 51, 78,
    -57, 0, 5,
    -57, 0, 3),
  ID = c(
    rep("ID107", 9),
    rep("ID142", 3),
    rep("ID104", 3),
    rep("ID112", 3)),
  r = c(160, 0, 232, 60, 29, 255, 25, 178, 0, 247, 139, 205, 76, 153, 230, 140, 192, 221),
  g = c(32, 197, 148, 179, 187, 105, 62, 34, 0, 6, 0, 105, 76, 153, 230, 127, 177, 212),
  b = c(240, 205, 190, 113, 229, 180, 130, 34, 205, 6, 0, 201, 76, 153, 230, 0, 0, 90)
)

rgb.v.colors$visit_colors <- rgb(rgb.v.colors[c("r", "g", "b")], maxColorValue = 255)

visit_colors <- rgb.v.colors$visit_colors
names(visit_colors) <- rgb.v.colors$visit_colors

### Donors' colors -------------------------------------------------------------

id.cols <- scales::hue_pal()(4)[c(2, 4, 1, 3)]
names(id.cols) <- c("ID107", "ID142", "ID104", "ID112")

### Standard theme settings ----------------------------------------------------

basic_theme <- theme(
  text = element_text(family = "Helvetica"),
  axis.text = element_text(size = 5),
  axis.title = element_text(size = 6),
  legend.text = element_text(size = 6),
  legend.title = element_text(size = 6,
                              vjust = 1),
  plot.title = element_text(size = 7,
                            face = "bold",
                            hjust = 0.5),
  legend.box.margin = margin(),
  legend.box.spacing = unit(2, "mm"),
  legend.margin = margin(),
  plot.margin = margin()
) 

### Features in RNA DotPlot ----------------------------------------------------

rna.feats <- c("CD3D",
               "CD8A",
               "CD4",
               "EGR2",
               "NR4A2",
               "NFKB1",
               "CRTAM",
               "TNFRSF9",
               "CD38",
               "SLAMF7",
               "XCL1",
               "XCL2",
               "GZMB",
               "GZMH",
               "IFNG",
               "CCL4",
               "CCL4L2",
               "ZBED2",
               "GNG4",
               "CCL5",
               "APOBEC3G",
               "IL2RA",
               "FOXP3",
               "TCF7",
               "IL7R",
               "KLRG1",
               "KIR2DL1",
               "FCGR3A",
               "TYROBP",
               "GNLY",
               "TRDC",
               "CX3CR1",
               "KLRB1",
               "TOX",
               "CD40LG",
               "TNFSF8",
               "DPP4",
               "LEF1",
               "ID3",
               "TNFRSF4",
               "GNA15",
               "CCL20",
               "CCR6",
               "CTSH",
               "OSM",
               "TNF",
               "SELL",
               "MKI67",
               "RRM2",
               "SPC25",
               "CCNA2",
               "MCM7",
               "LMNA",
               "TYMS",
               "SLC4A10",
               "TRAV1-2",
               "IL23R",
               "ZBTB16",
               "IL18RAP")

### Features in ADT Feature Plot -----------------------------------------------

adt.feats <- c("CD69", "CD107a", "CD25", "CD71")
