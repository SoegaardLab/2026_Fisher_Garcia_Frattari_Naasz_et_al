################################################################################
## Title        : Cluster Abundance Analysis - AIM+ cells
##
## Input        : Annotated AIM object
##
## Output       : Results of scanpro differential abundance analysis
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading libraries and data:##################################################

library(Seurat)
library(reticulate)
library(tidyverse)

# Activate Python environment
reticulate::use_virtualenv(Sys.getenv("PYTHON_ENV"))

## Environment setup:###########################################################

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

# Set seed
set.seed(123)

# Read object
AIM <- readRDS(paste0(processed_data, "steps/objects/aim_annotated.rds"))

# Summarize cluster abundance in a df
cluster.tab <- AIM@meta.data %>% 
  select(orig.ident, w_clusters, 
         donor, visit,
         clinical_time_point) %>% 
  tibble() %>% 
  rename(sample = orig.ident,
         cluster = w_clusters) %>% 
  mutate(status = factor(ifelse(donor %in% c("ID104", "ID112"),
                                "non_controller",
                                "controller"),
                         levels = c("controller", "non_controller")))

# Save as csv file for results and downstream
write.csv(cluster.tab,
          paste0(processed_data, "steps/dab/clustertab.csv"),
          row.names = F)

## Test differences in clusters with scanpro ###################################

# Run the python script
py_run_file("python/cluster_prop_aim.py")

### Inspect and summarize pre-ATI ----------------------------------------------

# Save results with full data
write.table(py$prop_art$results,
            file = paste0(processed_data, "results/dab/art_scanpro_full.txt"),
            col.names = T,
            row.names = F)

# Summarize bootstrapped results
art_results <- lapply(names(py$art_scanpro), function(x){
  
  py$art_scanpro[[x]]$results %>%
    rownames_to_column(var = "Cluster") %>%
    mutate(Cluster = factor(Cluster, levels = c(1:11))) %>%
    mutate(boot = as.numeric(x) + 1)
  
  
}) %>%
  do.call(rbind, .) %>%
  tibble()

art_summary <- art_results %>%
  group_by(Cluster) %>%
  summarise(`Baseline Proportion (Mean)` = mean(baseline_props),
            `Baseline Proportion (95% CI)` = paste(round(t.test(baseline_props)$conf.int[1], 3),
                                                   round(t.test(baseline_props)$conf.int[2], 3),
                                                   sep = "-"),
            `NC Proportion (Mean)` = mean(mean_props_non_controller),
            `NC Proportion (95% CI)` = paste(round(t.test(mean_props_non_controller)$conf.int[1], 3),
                                             round(t.test(mean_props_non_controller)$conf.int[2], 3),
                                             sep = "-"),
            `PIC Proportion (Mean)` = mean(mean_props_controller),
            `PIC Proportion (95% CI)` = paste(round(t.test(mean_props_controller)$conf.int[1], 3),
                                              round(t.test(mean_props_controller)$conf.int[2], 3),
                                              sep = "-"),
            `Proportion Ratio (Mean)` = mean(prop_ratio),
            `Proportion Ratio (95% CI)` = paste(round(t.test(prop_ratio)$conf.int[1], 3),
                                                round(t.test(prop_ratio)$conf.int[2], 3),
                                                sep = "-"),
            `Adjusted p-value (Median)` = median(adjusted_p_values),
            `Runs With p < 0.05 (Frequency)` = sum(adjusted_p_values < 0.05)/max(row_number()))

# Save to results
write.table(art_summary,
            file = paste0(processed_data, "results/dab/art_scanpro_boot.txt"),
            col.names = T,
            row.names = F)

# Summarize bootstrapped proportions
art_simulations <- lapply(names(py$art_boot), function(x){
  
  print(paste("Working on boot", x))
  
  df <- py$art_boot[[x]]
  
  df %>%
    count(cluster, status, donor) %>%
    mutate(boot = x)
  
}) %>%
  do.call(rbind, .)

art_boot_summary <- art_simulations %>%
  filter(donor != "ID112") %>%
  group_by(status, donor, boot) %>%
  mutate(prop = n/sum(n)) %>%
  group_by(cluster, status, donor) %>%
  summarise(`Proportion (Mean)` = mean(prop),
            `Proportion (95% CI)` = paste(round(t.test(prop)$conf.int[1], 3),
                                          round(t.test(prop)$conf.int[2], 3),
                                          sep = "-")) %>%
  ungroup() %>%
  mutate(Cluster = cluster,
         Group = ifelse(status == "controller", "PIC", "NC"),
         Donor = donor,
         .before = 3,
         .keep = "unused")

# Save to results
write.table(art_boot_summary,
            file = paste0(processed_data, "results/dab/art_scanpro_boot_proportions.txt"),
            col.names = T,
            row.names = F)

### Summarize ATI --------------------------------------------------------------

# Save results with full data
write.table(py$prop_ati$results,
            file = paste0(processed_data, "results/dab/ati_scanpro_full.txt"),
            col.names = T,
            row.names = F)

# Summarize bootstrapped results
ati_results <- lapply(names(py$ati_scanpro), function(x){
  
  py$ati_scanpro[[x]]$results %>%
    rownames_to_column(var = "Cluster") %>%
    mutate(Cluster = factor(Cluster, levels = c(1:11))) %>%
    mutate(boot = as.numeric(x) + 1)
  
  
}) %>%
  do.call(rbind, .) %>%
  tibble()

# Summarize bootstrapped proportions
ati_summary <- ati_results %>%
  group_by(Cluster) %>%
  summarise(`Baseline Proportion (Mean)` = mean(baseline_props),
            `Baseline Proportion (95% CI)` = paste(round(t.test(baseline_props)$conf.int[1], 3),
                                                   round(t.test(baseline_props)$conf.int[2], 3),
                                                   sep = "-"),
            `NC Proportion (Mean)` = mean(mean_props_non_controller),
            `NC Proportion (95% CI)` = paste(round(t.test(mean_props_non_controller)$conf.int[1], 3),
                                             round(t.test(mean_props_non_controller)$conf.int[2], 3),
                                             sep = "-"),
            `PIC Proportion (Mean)` = mean(mean_props_controller),
            `PIC Proportion (95% CI)` = paste(round(t.test(mean_props_controller)$conf.int[1], 3),
                                              round(t.test(mean_props_controller)$conf.int[2], 3),
                                              sep = "-"),
            `Proportion Ratio (Mean)` = mean(prop_ratio),
            `Proportion Ratio (95% CI)` = paste(round(t.test(prop_ratio)$conf.int[1], 3),
                                                round(t.test(prop_ratio)$conf.int[2], 3),
                                                sep = "-"),
            `Adjusted p-value (Median)` = median(adjusted_p_values),
            `Runs With p < 0.05 (Frequency)` = sum(adjusted_p_values < 0.05)/max(row_number()))

# Save to results
write.table(ati_summary,
            file = paste0(processed_data, "results/dab/ati_scanpro_boot.txt"),
            col.names = T,
            row.names = F)

# Summarize bootstrapped proportions
ati_simulations <- lapply(names(py$ati_boot), function(x){
  
  print(paste("Working on boot", x))
  
  df <- py$ati_boot[[x]]
  
  df %>%
    count(cluster, status, donor) %>%
    mutate(boot = x)
  
}) %>%
  do.call(rbind, .)

ati_boot_summary <- ati_simulations %>%
  filter(donor != "ID142") %>%
  group_by(status, donor, boot) %>%
  mutate(prop = n/sum(n)) %>%
  group_by(cluster, status, donor) %>%
  summarise(`Proportion (Mean)` = mean(prop),
            `Proportion (95% CI)` = paste(round(t.test(prop)$conf.int[1], 3),
                                          round(t.test(prop)$conf.int[2], 3),
                                          sep = "-")) %>%
  ungroup() %>%
  mutate(Cluster = cluster,
         Group = ifelse(status == "controller", "PIC", "NC"),
         Donor = donor,
         .before = 3,
         .keep = "unused")

# Save to results
write.table(ati_boot_summary,
            file = paste0(processed_data, "results/dab/ati_scanpro_boot_proportions.txt"),
            col.names = T,
            row.names = F)

## Summarized data for plots and supplementary #################################

### pre-ATI --------------------------------------------------------------------

non.control.high <- art_summary %>% 
  filter(`Adjusted p-value (Median)` < 0.05 & `Proportion Ratio (Mean)` > 1) %>% 
  pull(Cluster)

control.high <- art_summary %>% 
  filter(`Adjusted p-value (Median)` < 0.05 & `Proportion Ratio (Mean)` < 1) %>% 
  pull(Cluster)

art.clusters <- cluster.tab %>% 
  filter(visit %in% c("v13", "v02")) %>% 
  mutate(cluster = factor(cluster,
                          levels = c(11:1))) %>% 
  count(status,
        cluster) %>% 
  group_by(status) %>% 
  arrange(desc(cluster), .by_group = T) %>% 
  mutate(tot.cells = sum(n)) %>% 
  mutate(freq = n/tot.cells) %>% 
  mutate(label.position = cumsum(freq)-freq/2) %>% 
  mutate(lab = ifelse(freq < 0.04, 
                      "", 
                      forcats::fct_rev(cluster))) %>%
  mutate(lab = case_when(status == "controller" & cluster %in% control.high ~ paste0(lab, "*"),
                         status == "non_controller" & cluster %in% non.control.high ~ paste0(lab, "*"),
                         TRUE ~ lab)) %>% 
  ungroup() %>%
  mutate(status = ifelse(status == "controller",
                         "PICs",
                         "NCs")) %>% 
  mutate(status = factor(status,
                         levels = c("PICs", "NCs")))

art.clusters.by.donor <- cluster.tab %>% 
  filter(visit %in% c("v13", "v02")) %>%  
  count(sample, 
        donor,
        clinical_time_point,
        status,
        cluster) %>%
  group_by(sample) %>% 
  mutate(tot.cells = sum(n)) %>% 
  mutate(freq = n/tot.cells) %>% 
  ungroup() %>%
  mutate(status = ifelse(status == "controller",
                         "PICs",
                         "NCs")) %>% 
  mutate(status = factor(status,
                         levels = c("PICs", "NCs"))) %>% 
  mutate(cluster_label = ifelse(cluster == 1,
                                paste0(cluster, "*"),
                                cluster)) %>%
  mutate(cluster_label = factor(cluster_label,
                                levels = c(paste0(1, "*"),
                                           2:11))) %>% 
  mutate(cluster = factor(cluster,
                          levels = c(1:11))) %>% 
  mutate(clinical_time_point = "pre-ATI") %>% 
  mutate(donor = factor(donor,
                        levels = c("ID107", "ID142", "ID104", "ID112")))

### ATI ------------------------------------------------------------------------

# ATI cluster proportions by group
ati.clusters <- cluster.tab %>% 
  filter(visit %in% c("v24", "v11b")) %>% 
  mutate(cluster = factor(cluster,
                          levels = c(11:1))) %>% 
  count(status,
        cluster) %>% 
  group_by(status) %>% 
  arrange(desc(cluster), .by_group = T) %>% 
  mutate(tot.cells = sum(n)) %>% 
  mutate(freq = n/tot.cells) %>% 
  mutate(label.position = cumsum(freq)-freq/2) %>% 
  mutate(lab = ifelse(freq < 0.04, 
                      "", 
                      forcats::fct_rev(cluster))) %>%
  ungroup() %>%
  mutate(status = ifelse(status == "controller",
                         "PICs",
                         "NCs")) %>% 
  mutate(status = factor(status,
                         levels = c("PICs", "NCs")))

# ATI cluster proportions by donor
ati.clusters.by.donor <- cluster.tab %>% 
  filter(visit %in% c("v24", "v11b")) %>%  
  count(sample, 
        donor,
        clinical_time_point,
        status,
        cluster) %>%
  group_by(sample) %>% 
  mutate(tot.cells = sum(n)) %>% 
  mutate(freq = n/tot.cells) %>% 
  ungroup() %>%
  mutate(status = ifelse(status == "controller",
                         "PICs",
                         "NCs")) %>% 
  mutate(status = factor(status,
                         levels = c("PICs", "NCs"))) %>% 
  mutate(cluster_label = factor(cluster,
                                levels = c(1:11))) %>% 
  mutate(cluster = factor(cluster,
                          levels = c(1:11))) %>% 
  mutate(clinical_time_point = "ATI") %>% 
  mutate(donor = factor(donor,
                        levels = c("ID107", "ID142", "ID104", "ID112")))

### Export data for supplementary file -----------------------------------------
art.clusters.by.donor %>% 
  rbind(ati.clusters.by.donor) %>% 
  mutate(cluster = sub("\\*", "", cluster)) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  select(clinical_time_point, cluster, status, donor, n, tot.cells, freq) %>%
  arrange(desc(clinical_time_point), cluster, status, donor, n, tot.cells, freq) %>% 
  rename(`Time point` = clinical_time_point,
         Cluster = cluster,
         Group = status,
         Donor = donor,
         `Number of cells in Cluster` = n,
         `Number of cells in Sample` = tot.cells,
         `Frequency of Cluster in Sample` = freq) %>% 
  write.table(file = paste0(processed_data, "results/dab/dab_analysis_summary_data.txt"),
              dec = ",",
              row.names = F)

## Collect dab results in one rds to add to Seurat #############################

dab.results <- list(cluster_proportions_art = art.clusters.by.donor,
                    cluster_proportions_art_group = art.clusters,
                    cluster_proportions_ati = ati.clusters.by.donor,
                    cluster_proportions_ati_group = ati.clusters,
                    scanpro_results_full_data = list(art = py$prop_art$results,
                                                     ati = py$prop_ati$results),
                    cluster_proportions_bootstrap = list(art = art_boot_summary,
                                                         ati = ati_boot_summary),
                    scanpro_results_bootstrap = list(art = art_summary,
                                                     ati = ati_summary))

saveRDS(dab.results, paste0(processed_data, "steps/dab/dab_results.rds"))
