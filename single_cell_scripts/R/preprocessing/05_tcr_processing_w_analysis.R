################################################################################
## Title        : TCR pre-processing - AIM+ cells
##
## Input        : 10x filtered contig csv file
##
## Output       : CSV file with TCR genes and cdr3 info at single cell level
##
## Author       : Giacomo S Frattari
################################################################################
##
## Loading required libraries:##################################################

library(tidyverse)

## Environment setup:###########################################################

# Get data locations from environment
processed_data <- Sys.getenv("PROCESSED_DATA")
raw_data <- Sys.getenv("RAW_DATA")

samples <- c(paste("aim", "ID104", c("v01", "v13", "v24"), sep = "_"),
             paste("aim", "ID107", c("v01", "v13", "v24", "v28", "v34", "v37", "v40"), sep = "_"),
             paste("aim", "ID112", c("v01", "v13", "v24"), sep = "_"),
             paste("aim", "ID142", c("v02", "v11b", "v11c"), sep = "_"))

AIM <- readRDS(paste0(processed_data, "steps/objects/aim_annotated.rds"))

################################################################################
# TCR pre-processing ###########################################################
################################################################################

## Read TCR data ###############################################################

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
                            barcode)) %>% 
    mutate(donor = sub(".+(ID...).+", 
                       "\\1", 
                       barcode),
           visit = sub(".+(v.+)_.+", 
                       "\\1", 
                       barcode),
           .after = 1)
  
}) %>% 
  # Bind to one data frame
  do.call(rbind, .)

# Filter barcodes: keep only barcodes that are in the final object
tcr.tab <- tcr.tab %>% 
  filter(barcode %in% colnames(AIM))

## Identify clonotypes #########################################################

### Identify TRA and TRB cts ---------------------------------------------------

tra <- tcr.tab %>%
  filter(chain == "TRA") %>% 
  select(barcode, v_gene, j_gene, cdr3_nt) %>% 
  # Pool genes to name genetic TRA CT
  mutate(ct = paste0(v_gene, ".", j_gene, "_", cdr3_nt),
         .keep = "unused") %>% 
  # Mark TRA1 and TRA2 if multiple chains available
  group_by(barcode) %>% 
  # Arrange according to CT name to ensure systeatic ordering - important for 
  # later matching
  arrange(ct, .by_group = T) %>% 
  mutate(ct.tra = paste0("TRA", row_number())) %>%
  # Transform to wide format - collect chains that belong to same cell
  # All cells will have one TRA1, some cells will have no TRA2
  pivot_wider(names_from = ct.tra,
              values_from = ct) %>% 
  ungroup()

# Same procedure as tra
trb <- tcr.tab %>%
  filter(chain == "TRB") %>% 
  select(barcode, v_gene, j_gene, cdr3_nt) %>% 
  mutate(ct = paste0(v_gene, ".", j_gene, "_", cdr3_nt),
         .keep = "unused") %>% 
  group_by(barcode) %>% 
  arrange(ct, .by_group = T) %>% 
  mutate(ct.trb = paste0("TRB", row_number())) %>%
  pivot_wider(names_from = ct.trb,
              values_from = ct)  %>% 
  ungroup()

# Match TRA and TRB to cells
trab.tab <- full_join(tra, trb)

### Create CT identifiers: pool TRA and TRB together ---------------------------

ct.tab <- trab.tab %>% 
  mutate(CTa = case_when(is.na(TRA1) ~ NA, # If no TRA available, return NA
                         is.na(TRA2) ~ TRA1, # If no TRA2 available, returns only TRA1
                         TRUE ~ paste0(TRA1, ";", TRA2))) %>% # If multiple TRAs available, paste them
  mutate(CTb = case_when(is.na(TRB1) ~ NA,
                         is.na(TRB2) ~ TRB1,
                         TRUE ~ paste0(TRB1, ";", TRB2))) %>% 
  mutate(CT = case_when(is.na(CTa) ~ CTb, # If no CTa available, return CTb only
                        is.na(CTb) ~ CTa, # If no CTb available, return CTa only
                        TRUE ~ paste0(CTa, "+", CTb)), # If both available, paste them
         .keep = "unused")

### Process cdr3aa chain information -------------------------------------------

cdr3.a <- tcr.tab %>% 
  filter(chain == "TRA") %>% 
  select(barcode, cdr3) %>% 
  # Mark TRA1 and TRA2 if multiple chains available
  group_by(barcode) %>% 
  # Arrange according to cdr3
  arrange(cdr3, .by_group = T) %>% 
  mutate(cdr3.tra = paste0("TRA", row_number(), "_cdr3aa")) %>%
  # Transform to wide format - collect chains that belong to same cell
  # All cells will have one TRA1, some cells will have no TRA2
  pivot_wider(names_from = cdr3.tra,
              values_from = cdr3) %>% 
  ungroup() %>% 
  mutate(length_a1 = nchar(TRA1_cdr3aa)) %>% 
  mutate(length_a2 = nchar(TRA2_cdr3aa))

cdr3.b <- tcr.tab %>% 
  filter(chain == "TRB") %>% 
  select(barcode, cdr3) %>% 
  # Mark TRA1 and TRA2 if multiple chains available
  group_by(barcode) %>% 
  # Arrange according to cdr3
  arrange(cdr3, .by_group = T) %>% 
  mutate(cdr3.trb = paste0("TRB", row_number(), "_cdr3aa")) %>%
  # Transform to wide format - collect chains that belong to same cell
  # All cells will have one TRA1, some cells will have no TRA2
  pivot_wider(names_from = cdr3.trb,
              values_from = cdr3) %>% 
  ungroup() %>% 
  mutate(length_b1= nchar(TRB1_cdr3aa)) %>% 
  mutate(length_b2 = nchar(TRB2_cdr3aa))

# Add is not NA function to calculate n. of chain. If the chain exists, returns
# 1, otherwise 0

is.not.na <- function(x){
  
  abs(is.na(x) -1)
  
}

# Merge cdr3a and cdr3b and calculate number of chains
cdr3.ab <- full_join(cdr3.a, cdr3.b) %>% 
  mutate(n.tra = is.not.na(length_a1) + is.not.na(length_a2)) %>%
  mutate(n.trb = is.not.na(length_b1) + is.not.na(length_b2)) %>%
  mutate(n.chains = n.tra + n.trb) %>%
  mutate(chain = case_when (
    n.tra > n.trb ~ "TRA > TRB",
    n.trb > n.tra ~ "TRB > TRA",
    n.tra == n.trb ~ "pair"))

### Identify TRA and TRB cts ---------------------------------------------------

cdr3.tra <- tcr.tab %>%
  filter(chain == "TRA") %>% 
  select(barcode, v_gene, j_gene, cdr3) %>% 
  # Pool genes to name genetic TRA CT
  mutate(ct = paste0(v_gene, ".", j_gene, "_", cdr3),
         .keep = "unused") %>% 
  # Mark TRA1 and TRA2 if multiple chains available
  group_by(barcode) %>% 
  # Arrange according to CT name to ensure systeatic ordering - important for 
  # later matching
  arrange(ct, .by_group = T) %>% 
  mutate(ct.tra = paste0("TRA", row_number())) %>%
  # Transform to wide format - collect chains that belong to same cell
  # All cells will have one TRA1, some cells will have no TRA2
  pivot_wider(names_from = ct.tra,
              values_from = ct) %>% 
  ungroup()

# Same procedure as tra
cdr3.trb <- tcr.tab %>%
  filter(chain == "TRB") %>% 
  select(barcode, v_gene, j_gene, cdr3) %>% 
  mutate(ct = paste0(v_gene, ".", j_gene, "_", cdr3),
         .keep = "unused") %>% 
  group_by(barcode) %>% 
  arrange(ct, .by_group = T) %>% 
  mutate(ct.trb = paste0("TRB", row_number())) %>%
  pivot_wider(names_from = ct.trb,
              values_from = ct)  %>% 
  ungroup() 

# Match TRA and TRB to cells
cdr3.trab.tab <- full_join(cdr3.tra, cdr3.trb)

### Create CT identifiers: pool TRA and TRB together ---------------------------

cdr3.ct.tab <- cdr3.trab.tab %>% 
  mutate(CTa = case_when(is.na(TRA1) ~ NA, # If no TRA available, return NA
                         is.na(TRA2) ~ TRA1, # If no TRA2 available, returns only TRA1
                         TRUE ~ paste0(TRA1, ";", TRA2))) %>% # If multiple TRAs available, paste them
  mutate(CTb = case_when(is.na(TRB1) ~ NA,
                         is.na(TRB2) ~ TRB1,
                         TRUE ~ paste0(TRB1, ";", TRB2))) %>% 
  mutate(CT = case_when(is.na(CTa) ~ CTb, # If no CTa available, return CTb only
                        is.na(CTb) ~ CTa, # If no CTb available, return CTa only
                        TRUE ~ paste0(CTa, "+", CTb)), # If both available, paste them
         .keep = "unused") %>% 
  select(barcode, CT) %>% 
  rename(CTaa = CT)

## Create master tab with cts, cdr3s, genes and n.chains #######################

tcr.master <- full_join(cdr3.ab, ct.tab) %>% 
  full_join(cdr3.ct.tab) %>% 
  # Split CT informations back to cdr3nt and TCR genes
  mutate(CTdata = CT) %>% 
  # Split TRA and TRB
  separate(CTdata,
           into = c("TRA", "TRB"),
           sep = "\\+",
           fill = "right") %>%
  # Move TRB info back to TRB (moved to TRA if only one chain present)
  mutate(TRB = ifelse(n.chains == 1 & chain == "TRB > TRA",
                      TRA,
                      TRB)) %>%
  mutate(TRA = ifelse(n.chains == 1 & chain == "TRB > TRA",
                      NA,
                      TRA)) %>%
  # Split into two chains if two chains are present
  separate(TRA,
           into = c("TRA1", "TRA2"),
           sep = ";") %>% 
  separate(TRB,
           into = c("TRB1", "TRB2"),
           sep = ";") %>%
  # Split into gene and cdr3 sequences
  separate(TRA1,
           into = c("TRA1_gene", "TRA1_cdr3"),
           sep = "_") %>%   
  separate(TRB1,
           into = c("TRB1_gene", "TRB1_cdr3"),
           sep = "_") %>% 
  # Split genes into V and J genes
  separate(TRA1_gene,
           into = c("TRAV_1", "TRAJ_1"),
           sep = "\\.") %>%
  separate(TRA2,
           into = c("TRA2_gene", "TRA2_cdr3"),
           sep = "_") %>% 
  separate(TRA2_gene,
           into = c("TRAV_2", "TRAJ_2"),
           sep = "\\.") %>%
  separate(TRB1_gene,
           into = c("TRBV_1", "TRBJ_1"),
           sep = "\\.") %>% 
  separate(TRB2,
           into = c("TRB2_gene", "TRB2_cdr3"),
           sep = "_") %>% 
  separate(TRB2_gene,
           into = c("TRBV_2", "TRBJ_2"),
           sep = "\\.") %>% 
  # Add donor and visit information from barcodes for later merging
  mutate(donor = sub(".+(ID...).+", 
                     "\\1", 
                     barcode),
         visit = sub(".+(v.+)_.+", 
                     "\\1", 
                     barcode),
         .after = 1)

## Create output ###############################################################

### For downstream use ---------------------------------------------------------
ct.final <- tcr.master %>% 
  # Define order
  select(barcode, donor, visit,
         TRAV_1, TRAJ_1, TRA1_cdr3,
         TRAV_2, TRAJ_2, TRA2_cdr3,
         TRBV_1, TRBJ_1, TRB1_cdr3,
         TRBV_2, TRBJ_2, TRB2_cdr3,
         TRA1_cdr3aa, TRA2_cdr3aa, length_a1, length_a2,
         TRB1_cdr3aa, TRB2_cdr3aa, length_b1, length_b2,
         n.tra, n.trb, n.chains, chain,
         CT, CTaa)

################################################################################
# TCR analysis #################################################################
################################################################################

## Clonal ranks ################################################################

# Make cluster labels to a separate data frame
cluster.info <- AIM@meta.data[c("w_clusters", "donor")] %>% 
  rownames_to_column(var = "barcode")

ranked.ct <- ct.final %>%
  # Keep only CT with paired chains
  filter(n.chains > 1) %>% 
  count(donor, CT, name = "prevalence") %>%
  # Define rank
  group_by(donor) %>% 
  distinct(donor, CT, .keep_all = T) %>% 
  arrange(desc(prevalence), .by_group = T) %>% 
  mutate(rank = row_number()) %>% 
  mutate(rank_group = case_when(rank >= 1 & rank <= 10 ~ "1:10",
                                rank >= 11 & rank <= 50 ~ "11:50",
                                rank >= 51 & rank <= 200 ~ "51:200",
                                rank >= 201 & rank <= 500 ~ "201:500",
                                rank >= 501 ~ "> 500")) %>% 
  mutate(rank_group = factor(rank_group,
                             levels = c("1:10",
                                        "11:50",
                                        "51:200",
                                        "201:500",
                                        "> 500"))) %>%
  # Rejoin full table with tcr data and cluster
  left_join(ct.final) %>% 
  left_join(cluster.info) %>% 
  mutate(control = factor(ifelse(donor %in% c("ID112", "ID104"), 
                                 "NCs",
                                 "PICs"),
                          levels = rev(c("NCs",
                                         "PICs"))))

ranked.summarized <- ranked.ct %>% 
  count(control, donor, w_clusters, rank_group) %>% 
  group_by(control, donor, w_clusters) %>% 
  mutate(frequency = (n/sum(n)))

### Top 10 ranked across visits ------------------------------------------------

top.10.cts <- ct.final %>% 
  # Keep only clones with > 2 chains
  filter(n.chains > 1) %>%
  # Count
  group_by(donor, CT) %>% 
  count() %>%
  # Ransk
  group_by(donor) %>% 
  arrange(desc(n), .by_group = T) %>% 
  mutate(rank = row_number()) %>% 
  # Select top10
  slice_min(n = 10,
            order_by = rank,
            with_ties = F) %>%  
  select(donor, CT, rank)

# Basisplot til alluvial
top10.development <- ct.final %>% 
  # Calculate frequencies
  group_by(donor, visit, CT) %>% 
  count() %>% 
  group_by(donor, visit) %>% 
  mutate(total = sum(n)) %>% 
  mutate(frequency = n/total) %>% 
  ungroup() %>% 
  # Select only top 10
  left_join(top.10.cts) %>% 
  filter(is.na(rank) == F) %>%
  # Add control info
  mutate(control = ifelse(donor %in% c("ID104", "ID112"),
                          "NCs",
                          "PICs")) %>% 
  mutate(control = factor(control,
                          levels = c("PICs",
                                     "NCs"))) %>% 
  mutate(rank = factor(rank,
                       levels = c(1:10)))

## Clonality calculation #######################################################

#### Calculate Gni coefficients with bootstrapping -----------------------------

# Calculate prevalence of each CTaa in each cluster for each sample, filter out
# too small repertoires, find size of smallest remaining repertoires
ct.min <- ct.final %>% 
  left_join(cluster.info) %>% 
  # Remove cells with only one chain
  filter(n.chains > 1) %>% 
  # Count number per cluster per sample
  count(CT, donor, w_clusters, visit) %>% 
  # Count total number of cells pr visit for each cluster
  group_by(donor, w_clusters, visit) %>% 
  mutate(n_visit = sum(n)) %>% 
  # Remove too small clusters
  filter(n_visit >= 50) %>% 
  # Find the minimum number of cells in each cluster (sample with lowest output)
  # for later bootstrapping
  group_by(w_clusters) %>% 
  mutate(n_anno = min(n_visit)) %>% 
  ungroup()

# Create table with one entry = one cell to be used in calculations
ct.groups <- ct.final %>% 
  left_join(cluster.info) %>% 
  # Keep only CTs selected in ct.min
  filter(CT %in% ct.min$CT) %>% 
  # Add ct.min info
  left_join(ct.min) %>% 
  # Remove cells with CT/cluster combinations that were below limit
  filter(is.na(n_anno) == F) %>% 
  # Group for later split
  group_by(donor, w_clusters, visit) %>% 
  # Keep only relevant columns
  select(CT, donor, w_clusters, visit, n_anno)

# Split the data frame according to the grouping: one data frame pr cluster pr
# sample
ct.list <- group_split(ct.groups)

# Extract grouping names to create unique id's for the ct list
group.names <- ct.groups %>% 
  group_keys() %>% 
  mutate(name = paste(donor, visit, w_clusters, sep = "_")) %>% 
  pull(name)

names(ct.list) <- group.names

# Calculate gini index in each cluster with bootstrapping (minimize effect of
# size on Gini index - at least to compare between visits within same cluster)
ct.ginis <- lapply(ct.list, function(x) {
  
  print(paste(unique(x$donor), 
              "cluster", unique(x$w_clusters), 
              unique(x$visit)))
  
  # Bootstrap 100 times
  gini.values <- vector(length = 100)
  
  set.seed(123)
  
  for (i in 1:100){
    
    to.test <- slice_sample(x,
                            n = as.numeric(unique(x$n_anno)))
    
    test.counts <- to.test %>% 
      group_by(CT) %>% 
      count()
    
    gini.values[i] <- DescTools::Gini(test.counts$n)
    
  }
  
  out <- mean(gini.values) %>% as.data.frame()
  
  colnames(out) <- "gini"
  
  return(out)
  
})

# Create final tab for gini values
gini.results <- do.call(rbind, ct.ginis) %>% 
  rownames_to_column(var = "id") %>% 
  separate(id,
           into = c("donor", "visit", "w_clusters"),
           sep = "_") %>% 
  group_by(donor, w_clusters) %>% 
  # Calculate median value for each cluster in each donor
  summarise(gini = median(gini, na.rm = T), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(w_clusters = factor(w_clusters,
                             levels = c(1:11))) %>% 
  mutate(donor = factor(donor,
                        levels = c("ID107", "ID142", "ID104", "ID112"))) %>% 
  mutate(control = factor(ifelse(donor %in% c("ID112", "ID104"), 
                                 "NCs",
                                 "PICs"),
                          levels = c("PICs",
                                     "NCs")))

## Repertoire overlap ##########################################################

#### Calculate Morisita-Horn indices to quantify overlap -----------------------

# Morisita-Horn index of similarity
mh.index <- function(tab, id.1, id.2){
  
  X <- colSums(tab[id.1])
  Y <- colSums(tab[id.2])
  
  num <- 2 * sum(tab[,id.1] * tab[,id.2])
  den <- ((sum(tab[,id.1]^2)/(X^2)) + (sum(tab[,id.2]^2)/(Y^2))) * X * Y
  
  return(num/den)
  
}

# Overlaps between clusters (1 cluster = 1 annotation in 1 donor)

# Create frequency tab for clonotypes in each donor_cluster combination
anno.overlap.tab <- ct.final %>%
  # Filter to contain only cells with paired chains
  filter(n.chains > 1) %>% 
  # Add cluster info
  left_join(cluster.info) %>% 
  # Count CT frequency per donor/cluster
  count(donor, w_clusters, CT) %>%
  # Create unique identifier
  mutate(annori = paste(donor, w_clusters, sep = "_")) %>% 
  # Prepare for calculation
  select(CT, annori, n) %>% 
  pivot_wider(names_from = annori,
              values_from = n,
              values_fill = 0)

# Retrieve identifier to test (i.e. all column names except the first with CT id)
id.to.test <- grep("CT", colnames(anno.overlap.tab), value = T, invert = T)

# Create comparison tab: unique combinations of overlaps tested
comparisons <- t(combn(id.to.test, 2)) %>% data.frame()

colnames(comparisons) <- c("id1", "id2")

# Initialize and calculate morisita index column
comparisons$morisita <- 0

for(i in 1:nrow(comparisons)){
  
  comparisons$morisita[i] <- mh.index(anno.overlap.tab , 
                                      comparisons[i,"id1"], 
                                      comparisons[i,"id2"])
  
}

# Wrangle to table
anno.overlap.tab <- comparisons %>% 
  # Extract first cluster id
  mutate(annotation1 = sub(".+_", "", id1),
         donor1 = sub("(.+)_.+", "\\1", id1),
         .keep = "unused") %>% 
  # Extract second cluster id
  mutate(annotation2 = sub(".+_", "", id2),
         donor2 = sub("(.+)_.+", "\\1", id2),
         .keep = "unused") %>%
  # Select only relevant columns
  select (annotation1, donor1,
          annotation2, donor2,
          morisita) %>% 
  # Keep only comparisons within same donor and between different clusters
  filter(donor1 == donor2 & annotation1 != annotation2) %>%
  # Add control metadata
  mutate(control = factor(ifelse(donor1 %in% c("ID112", "ID104"), 
                                 "NCs",
                                 "PICs"),
                          levels = rev(c("NCs",
                                         "PICs"))))

################################################################################
## Export data #################################################################
################################################################################

### Collected in one RDS to add to final object --------------------------------

tcr_results <- list(tcr_processed = ct.final,
                    tcr_w_ranks = ranked.ct,
                    tcr_w_ranks_summarized = ranked.summarized,
                    tcr_10_ranks = top.10.cts,
                    tcr_10_development = top10.development,
                    tcr_clonality = gini.results,
                    tcr_overlap = anno.overlap.tab)

saveRDS(tcr_results, paste0(processed_data, "steps/tcr/tcr_data.rds"))
