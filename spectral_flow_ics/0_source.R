# LOAD LIBRARIES
library(ggh4x)
library(data.table)
library(formulaic)
library(flowCore)
library(flowViz)
library(PeacoQC)
library(CATALYST)
library(SingleCellExperiment)
library(uwot)
library(knitr)
library(writexl)
library(stringr)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(glue)
library(cowplot)
library(scater)
library(patchwork)
library(ComplexHeatmap)
library(flowWorkspace)
library(scales)
library(ggcyto)
library(CytoExploreR)
library(MetaCyto)
library(cyCombine)
library(outliers)
library(clustree)
library(ggrepel)
library(cytoqc)
library(gridExtra)
library(grid)
library(tibble)
library(RColorBrewer)
library(ggridges)
library(ggtext)


# SET AND WRITE DIRECTORIES
dir_project = dirname(getwd())
scriptname = str_replace(basename(rstudioapi::getSourceEditorContext()$path),".Rmd|.R","")
dir_pipeline = glue("{dir_project}/2_pipeline/{scriptname}")
for (folder in c("tmp","store","out")) dir.create(glue("{dir_pipeline}/{folder}"), showWarnings = FALSE, recursive = TRUE)
dir_tmp = glue("{dir_pipeline}/tmp")
dir_store = glue("{dir_pipeline}/store")
dir_out = glue("{dir_pipeline}/out")
rm(list = c("folder","scriptname"))

# LOAD FUNCTIONS

q_AddExperimentInfo = function(fs){
  
  if(any(grepl(" ",sampleNames(fs)))){
    sample_names_clean = str_replace_all(sampleNames(fs)," ","_")
    sampleNames(fs) = sample_names_clean
    pData(fs)$name = rownames(pData(fs))
  }
  
  
  # find filenames
  filenames = sampleNames(fs)
  
  # find well id
  well_ids = well_ids = sapply(str_split(filenames,"_"),"[[",1)
  
  # find dates
  dates = fsApply(fs, keyword, "$DATE")
  filedates = matrix(unlist(dates), ncol=1, byrow=FALSE)
  filedates = lubridate::dmy(filedates)
  # sorted dates
  filedates_sorted = unique(sort(filedates))
  
  # define batch
  batches = factor(filedates, levels = filedates_sorted)
  levels(batches) = LETTERS[1:length(unique(batches))]
  
  # make df
  ei = data.frame(name = filenames,
                  well_id = well_ids,
                  batch = batches,
                  date = filedates)
  
  rownames(ei) = ei$name
  
  pData(fs) = ei
  
  pData(fs)$sample_id = str_replace(q_SplitExtractString(q_SplitExtractString(pData(fs)$name, "_WLSM",1),(glue("{pData(fs)$well_id}_")),2),"_MC","")
  
  pData(fs) = pData(fs)[,c("date","batch","sample_id","well_id")]
  
  return(fs)
}

q_DefinePanel = function(x,list_type,list_state){
  
  # get column names
  fcs_colname = colnames(x)
  
  # get antigen from description
  antigen = pData(parameters(x[[1]]))$desc
  antigen = sapply(str_split(antigen, " - ",2),"[",1)
  antigen[c(grep(" : ",antigen,value=FALSE,invert=TRUE))]=NA
  antigen = sapply(str_split(antigen, " : ",2),"[",1)
  
  # table of panel
  panel = data.frame(fcs_colname, antigen, row.names = NULL)
  
  # create dict
  values_type = rep("type", length(list_type))
  values_state = rep("state", length(list_state))
  keys = c(list_type,list_state)
  values = c(values_type,values_state)
  names(values) = keys
  
  # add column of marker_class
  panel = panel %>% 
    mutate(marker_class = recode(antigen, !!!values, .default="none")) %>% 
    replace_na(list(marker_class = "none"))
  
  return(panel)
}

q_densityplot = function(parameter, data_ff_fs, ...){
  
  if (is(data_ff_fs, "flowFrame")){
    
    minRange_maxRange_df = range(data_ff_fs, type = "data")
    minRange_maxRange_df_t = data.frame(t(minRange_maxRange_df))
    rownames(minRange_maxRange_df_t) = colnames(minRange_maxRange_df)
    colnames(minRange_maxRange_df_t) = c("minRange","maxRange")
    
    data_ff_fs@parameters@data$minRange = minRange_maxRange_df_t$minRange
    data_ff_fs@parameters@data$maxRange = minRange_maxRange_df_t$maxRange
    
  }
  
  if (is(data_ff_fs, "flowSet")) {
    
    n_flowframes = nrow(pData(data_ff_fs))
    
    for (flowframe in (1:n_flowframes)){
      
      minRange_maxRange_df = range(data_ff_fs[[flowframe]], type = "data")
      minRange_maxRange_df_t = data.frame(t(minRange_maxRange_df))
      rownames(minRange_maxRange_df_t) = colnames(minRange_maxRange_df)
      colnames(minRange_maxRange_df_t) = c("minRange","maxRange")
      
      data_ff_fs[[flowframe]]@parameters@data$minRange = minRange_maxRange_df_t$minRange
      data_ff_fs[[flowframe]]@parameters@data$maxRange = minRange_maxRange_df_t$maxRange
      
    }
    
  }
  
  return(flowViz::densityplot(parameter, data_ff_fs, ...))
  
}

q_TransformEstimateLogicle = function(fs, ref_samplename, panel = panel, seed = 1234, plot = F, add_suffix = "_tf"){
  
  # seed for reproducibility
  set.seed(seed)
  
  # define markers from panel
  markerstotransform = panel[!is.na(panel$antigen), ]$fcs_colname
  
  # calculate transformation
  max_value = max(fsApply(fs, each_col, max))
  trans = estimateLogicle(fs[(pData(fs)$name==ref_samplename)][[1]],
                          channels = markerstotransform)
  
  # transform flowset
  fs_transform = transform(fs, trans)
  
  # save plots
  if(plot == TRUE){
    
    p1 = q_densityplot(~`.`,
                       channels = markerstotransform,
                       fs[(pData(fs)$name==ref_samplename)][[1]],
                       main = "Before transformation")
    p2 = q_densityplot(~`.`,
                       channels = markerstotransform,
                       fs_transform[(pData(fs_transform)$name==ref_samplename)][[1]],
                       main = "After transformation")
    
    pdf(file = glue("{dir_store}/Transformation_before.pdf"), width = 15, height = 12)
    print(p1)
    dev.off()
    
    pdf(file = glue("{dir_store}/Transformation_after.pdf"), width = 15, height = 12)
    print(p2)
    dev.off()
    
    print(p1)
    print(p2)
    
  }
  
  if(is.null(add_suffix)==FALSE){
    fs_transform = q_PrefixSuffixSampleNames(fs_transform, "suffix", add_suffix)
  }
  
  
  return(fs_transform)
  
}

q_SplitExtractString = function(string, split_by, number_element){
  
  extracted_string = sapply(strsplit(string,
                                     split = split_by),
                            "[[",
                            number_element)
  
  return(extracted_string)
}

q_PrefixSuffixSampleNames = function(fs, add = c("prefix","suffix"), to_add){
  if(add == "prefix"){
    filenames_new = glue("{to_add}{sampleNames(fs)}.fcs")
  }
  
  if(add == "suffix"){
    filenames = str_replace(sampleNames(fs),".fcs","")
    filenames_new = glue("{filenames}{to_add}.fcs")
  }
  
  sampleNames(fs) = filenames_new
  return(fs)
}

q_MergeMetadataFlowSet = function(fs, dataframe, merge_col){
  
  pData(fs)$name = rownames(pData(fs))
  
  pd = merge(x = pData(fs),
             y = dataframe,
             by = merge_col,
             all.x = TRUE)
  
  rownames(pd) = pd$name
  pData(fs) = pd
  
  return(fs)
  
}

q_PeacoQC_FindParameters = function(fs, sample_name = pData(fs)$name[1], panel = panel, MAD_val, IT_val, outdir){
  
  peacoqc_res = PeacoQC(ff = fs[(pData(fs)$name == sample_name)][[1]],
                        channels = panel[panel$marker_class!="none",]$fcs_colname,
                        determine_good_cells = "all",
                        save_fcs = FALSE,
                        plot = TRUE,
                        output_directory = outdir,
                        IT_limit = IT_val,
                        MAD = MAD_val,
                        time_channel_parameter = "TIME",
                        time_unit = .1)
}

q_PeacoQC_RunQC = function(fs, panel = panel, MAD_val = 6, IT_val = .6, outdir){
  
  set.seed(1234)
  for(i in 1:nrow(pData(fs))){
    ff = fs[[i]]
    
    peacoqc_res = PeacoQC(ff,
                          channels = panel[panel$marker_class!="none",]$fcs_colname,
                          determine_good_cells = "all",
                          save_fcs = TRUE,
                          plot = TRUE,
                          output_directory = outdir,
                          IT_limit = IT_val,
                          MAD = MAD_val,
                          time_channel_parameter = "TIME",
                          time_unit = .1)
  }
  
  # read QC files
  fcs_dir = glue("{outdir}/PeacoQC_results/fcs_files")
  fs_clean = read.flowSet(path=fcs_dir,
                          pattern="*.fcs",
                          transformation = FALSE,
                          truncate_max_range = FALSE)
  
  fs_clean = q_AddExperimentInfo(fs_clean)
  
  return(fs_clean)
  
}


q_Gating_GetBoolMatrix = function(flowSet, gatingSet = NULL, gatingTemplate = NULL, gatenames = NULL, get_matrix = T){
  
  if(is.null(gatingSet)){
    
    gatingSet = GatingSet(flowSet)
    gatingSet = cyto_gatingTemplate_apply(gatingSet, gatingTemplate)
    
  }else{
    gatingSet = gatingSet
  }
  
  files_in_wsp = flowWorkspace::sampleNames(gatingSet)
  # get counts manually (these are normally as suffix when using the original function)
  counts2 = fsApply(flowSet, keyword, "$TOT")
  # check the order of files:
  stopifnot(names(counts2) == files_in_wsp)
  counts = as.numeric(unlist(counts2, use.names=FALSE))
  names(counts) = names(counts2)
  files_in_wsp = gsub("_[0-9]*$", "", files_in_wsp)
  result = list()
  for(file in files_in_wsp){
    
    if(is.null(gatenames)){
      
      cellTypes_tmp = sapply(str_split(gs_get_leaf_nodes(gatingSet),"/"),tail,1) 
      
    } else {
      
      cellTypes_tmp = gatenames
    }
    
    gate_names = cellTypes_tmp
    
    gatingMatrix = matrix(NA,
                          nrow = counts[file],
                          ncol = length(gate_names),
                          dimnames = list(NULL,
                                          gate_names))
    for(gate in gate_names){
      gatingMatrix[, gate] =
        flowWorkspace::gh_pop_get_indices(gatingSet[[file]], gate)
    }
    
    manual = FlowSOM::ManualVector(gatingMatrix, gate_names)
    
    result[[file]] = list("matrix" = gatingMatrix,
                          "manual" = manual)
    
  }
  
  
  gating = result
  
  if(get_matrix == T){
    # get boolean matrix
    matrices = sapply(gating, "[[", 1)
    combined_matrix = do.call(rbind,matrices)
    
    stopifnot(nrow(combined_matrix)==sum(counts))
    
    final_matrix = data.frame(combined_matrix)
    
    names(final_matrix) = gsub("[.]","_pos",names(final_matrix))
    
    return(final_matrix)
    
  }else{
    
    return(gating)
    
  }
  
}

q_DownsampleSCE = function(sce, samplesize){
  
  # split cell indices by sample
  cs = split(seq(ncol(sce)), sce$sample_id)
  # downsample to at most samplesize per sample
  cs = lapply(cs, \(.) sample(., min(samplesize, length(.))))
  sce_down = sce[, unlist(cs)]
  
  return(sce_down)
}

q_AnnotateClusters = function(sce, clusterList, new_id = "cell_type"){
  
  merging_table = data.frame(stack(clusterList)[2:1])
  merging_table = merging_table %>% rename(old_cluster = values,
                                           new_cluster = ind)
  merging_table = merging_table[,c("old_cluster","new_cluster")]
  
  # check validity
  repeated_metaclusters = c(unique(which(table(merging_table$old_cluster)>1)))
  for(repeated_metacluster in repeated_metaclusters){
    stop(paste(glue("Warning: metacluster {repeated_metacluster} is repeated in input list")))
  }
  
  meta_name = glue("meta{nrow(merging_table)}")
  print(glue("Renaming {meta_name} using following table:"))
  print(kable(merging_table))
  
  
  sce = mergeClusters(sce,
                      k = meta_name,
                      table = merging_table,
                      id = new_id,
                      overwrite = TRUE)
  
  return(sce)
  
}

q_plotDR = function(x, color_by, facet_by=NULL, labelsize=3, highlight_clust = NULL, background_color = "grey90", label_repel = TRUE,...){
  
  
  # create plot and df
  plt = plotDR(x,"UMAP",color_by = color_by,facet_by=facet_by, ...)
  df = plt$data
  
  # without facet_by
  if (is.null(facet_by)){
    df = aggregate(x = df[c("x", "y")],
                   by = list(c = get(color_by,df)),
                   FUN = median)
    
    no_clusters = length(unique(df$c))
    palette = plt[["plot_env"]][["k_pal"]][1:no_clusters]
    
    if(!is.null(highlight_clust)){
      
      all_clust = unique(df$c)
      grey_clusters = which(!all_clust %in% highlight_clust)
      palette[grey_clusters] = background_color
      
      if(!is.null(labelsize)){
        df$c = ifelse(df$c %in% highlight_clust, as.character(df$c), NA)
      }else{
        df$c = NA
        labelsize=0
      }
      
      
    }
    
    if(label_repel == TRUE){
      plt +
        geom_label_repel(aes(x, y, label = c),
                         size=labelsize,
                         col=palette,
                         data = df,
                         show.legend = FALSE,
                         inherit.aes = FALSE,
                         seed = 1234) +
        scale_colour_manual(values=palette)
    } else {
      plt +
        geom_label(aes(x, y, label = c),
                   size=labelsize,
                   col=palette,
                   data = df,
                   show.legend = FALSE,
                   inherit.aes = FALSE) +
        scale_colour_manual(values=palette)
    }
    
    
    
    # with facet_by
  } else {
    df = aggregate(x = df[c("x", "y")],
                   by = list(c = get(color_by,df),f=get(facet_by,df)),
                   FUN = median)
    
    no_clusters = length(unique(df$c))
    no_facets = length(unique(df$f))
    palette = plt[["plot_env"]][["k_pal"]][1:no_clusters]
    
    if(!is.null(highlight_clust)){
      
      all_clust = unique(df$c)
      grey_clusters = which(!all_clust %in% highlight_clust)
      palette[grey_clusters] = background_color
      
      if(!is.null(labelsize)){
        df$c = ifelse(df$c %in% highlight_clust, as.character(df$c), NA)
      }else{
        df$c = NA
        labelsize = 0
      }
      
    }
    
    df = df %>% rename(!!facet_by := "f")
    
    if(label_repel == TRUE){
      
      plt +
        geom_label_repel(data = df,
                         mapping = aes(x = x, y = y, label = c),
                         size=labelsize,
                         col=rep(palette,no_facets),
                         show.legend = FALSE,
                         inherit.aes = FALSE,
                         seed = 1234) +
        scale_colour_manual(values=palette)
    } else {
      
      plt +
        geom_label(data = df,
                   mapping = aes(x = x, y = y, label = c),
                   size=labelsize,
                   col=rep(palette,no_facets),
                   show.legend = FALSE,
                   inherit.aes = FALSE) +
        scale_colour_manual(values=palette)
    }
    
    
  }
  
  
}


q_SCEAddBoolMatrix = function(sce,bool_matrix){
  
  keep_cols = setdiff(names(colData(sce)), names(bool_matrix))
  new_coldata = cbind(colData(sce)[,keep_cols],bool_matrix)
  colData(sce) = new_coldata
  
  return(sce)
  
}

q_SCEBooleanCombinations = function(sce, colnames_combine, palette_legend1 = NULL){
  
  BoolMatrix = data.frame(colData(sce)[,colnames_combine])
  
  # get marker names from sce
  markers_gated = c()
  for(marker in rownames(sce)){
    if(any(grepl(marker,colnames_combine))){
      markers_gated = c(markers_gated,marker)
    }
  }
  print(markers_gated)
  markers_gated = setdiff(markers_gated,"CD3") 
  print(markers_gated)
  # rename columns in gated data to fit those in sce
  new_colnames = c()
  for(colname in names(BoolMatrix)){
    for(marker_gated in markers_gated){
      if(grepl(marker_gated,colname)){
        new_colnames = c(new_colnames,marker_gated)
      }
    }
  }
  names(BoolMatrix) = new_colnames
  
  # get order of combinations
  list_combinations = list()
  for(n in c(1:length(new_colnames))){
    combination_n = paste0(combn(new_colnames, n, paste0, collapse = "+"),"+")
    list_combinations = append(list_combinations, combination_n)
  }
  levels_combinations = unlist(list_combinations,recursive = F)
  
  # assign label to each cell
  df_long <- BoolMatrix %>%
    dplyr::mutate(Cell = row_number()) %>%  # Assign unique IDs for each cell
    pivot_longer(cols = -Cell,
                 names_to = "Marker", 
                 values_to = "Positive") %>%
    dplyr::group_by(Cell) %>%
    dplyr::summarize(Marker_combination = if (any(Positive)) {
      paste(Marker[Positive], collapse = "+")%>% paste0("+")  # Concatenate marker names for positives
    } else {
      "All negative"  # Assign "All negative" if no markers are TRUE
    }, .groups = "drop")
  
  # make factor 
  df_long$Marker_combination = factor(df_long$Marker_combination, levels = c("All negative",levels_combinations))
  
  # Add polyfunctionality columns to sce
  colData(sce)$Polyfunctionality_individual = df_long$Marker_combination
  colData(sce)$Polyfunctionality_grouped = str_count(df_long$Marker_combination,pattern = "[+]")
  
  levels_polyfunctionality_grouped = str_sort(as.character(unique(colData(sce)$Polyfunctionality_grouped)), numeric = TRUE)
  
  # make factor
  colData(sce)$Polyfunctionality_grouped = factor(colData(sce)$Polyfunctionality_grouped, levels = levels_polyfunctionality_grouped)
  
  # get counts
  counts_individual = data.frame(colData(sce)[,"Polyfunctionality_individual",drop=F]) %>% group_by(Polyfunctionality_individual)%>%count()
  counts_grouped = data.frame(colData(sce)[,"Polyfunctionality_grouped",drop=F]) %>% group_by(Polyfunctionality_grouped)%>%count()
  
  # make legend with marker combinations
  plus_matrix = data.frame()
  for(column_BoolComb in levels(colData(sce)$Polyfunctionality_individual)){
    for(marker in new_colnames){
      plus_matrix[marker,column_BoolComb] = grepl(marker,column_BoolComb)
    }
  }
  
  # plot legend with marker combinations
  plus_matrix[plus_matrix==TRUE] = 1
  plus_matrix[plus_matrix==FALSE] = 0
  plus_matrix = as.matrix(plus_matrix)
  
  col_annotation <- data.frame(Polyfunctionality_grouped = str_count(colnames(plus_matrix), "[+]"))
  unique_col_annos = unique(col_annotation$Polyfunctionality_grouped)
  if(is.null(palette_legend1)){
    col_palette = c("grey85",brewer_pal(type="qual")(length(unique_col_annos)-1))
  }else{
    col_palette = palette_legend1
  }
  names(col_palette) <- unique_col_annos
  
  col_annot <- HeatmapAnnotation(df = col_annotation, 
                                 col = list(Polyfunctionality_grouped = col_palette),
                                 show_legend = T,
                                 show_annotation_name = F)
  
  heatmap <- Heatmap(plus_matrix, 
                     row_names_side = "left", 
                     top_annotation = col_annot,
                     cluster_rows = F,
                     cluster_columns = F,
                     col = c("grey96","black"),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x, y, width, height,
                                 gp = gpar(fill = fill, col = "white", lwd = 1))},
                     show_heatmap_legend = F,
                     column_split = col_annotation,
                     column_title=NULL,
                     row_title=NULL,
                     show_column_names = F,
                     row_names_gp = gpar(fontsize = 8))
  
  p = heatmap
  
  # save results
  list_result = list()
  list_result$sce = sce
  list_result$counts_individual = counts_individual
  list_result$counts_grouped = counts_grouped
  list_result$legend = p
  
  return(list_result)
  
}



q_PsuedobulkPercPosMarker = function(dataframe, polyfunctionality = c("Polyfunctionality_grouped","Polyfunctionality_individual"), background_col = "Stimulation", background_value = "Neg"){
  
  all_polyfunctionality = unique(dataframe[,polyfunctionality])  # Extract all unique values of the polyfunctionality column
  
  group_cols = setdiff(names(dataframe), polyfunctionality) # Identify all grouping columns except the polyfunctionality column
  
  # Calculate frequency 
  df_freq = dataframe %>% count(!!!syms(group_cols), !!sym(polyfunctionality)) %>%  # Count dynamically using group columns
    complete(!!!syms(group_cols), !!sym(polyfunctionality) := all_polyfunctionality, fill = list(n = 0)) %>%  # Ensure all polyfunctionality values exist
    group_by(!!!syms(group_cols)) %>%
    mutate(
      Total_cells = sum(n),  # Sum of all cells in the group
      Percentage = (n / Total_cells) * 100  # Compute percentage
    ) %>%
    ungroup() %>%
    dplyr::filter(Total_cells > 0) # Remove rows with no cells
  
  # Subtrack background
  df_freq_background_subtract = df_freq %>% 
    group_by(across(-c(!!sym(background_col), n, Total_cells, Percentage))) %>%
    mutate(Percentage_backgroundSubtracted = pmax(Percentage - Percentage[!!sym(background_col) == background_value],0))
  
  # Check all backgrounds are zero after background subtraction 
  stopifnot(all(df_freq_background_subtract[df_freq_background_subtract[,background_col]==background_value,"Percentage_backgroundSubtracted"]==0))
  
  # Filter out the 0 and All negative polyfunctionality group
  df_final = df_freq_background_subtract %>% dplyr::filter(!(!!sym(polyfunctionality) %in% c("All negative","0")))
  
  return(df_final)
  
}

q_plotPeptideSpecificResponse_NC_PTC = function(dataframe, x_value = "Visit", fill_value = "Marker", facet_y = c("cell_type","cell_type_merge"),subtract_background = T, legend_title = NULL){
  
  if(subtract_background == T){
    dataframe = dataframe[dataframe$Stimulation!="Neg",]
    dataframe[,"y_value"] = dataframe[,"percent_pos_minusBackground"]
  }else{
    dataframe[,"y_value"] = dataframe[,"percent_pos"]
  }
  
  if(is.null(legend_title)){
    legend_title = "Marker positivity"
  }
  
  # data cleanup
  n_stims = length(unique(dataframe$Stimulation))
  n_PIDs = length(unique(dataframe$PID))
  dataframe[,"facet_y"] = dataframe[,facet_y]
  dataframe[,"x_value"] = dataframe[,x_value]
  dataframe[,"fill_value"] = dataframe[,fill_value]
  
  dataframe$Stimulation = factor(dataframe$Stimulation, levels = c("Neg", "Gag", "Pol", "Nef", "Env"))
  
  p = dataframe %>% ggplot(aes(x = x_value,
                               y = y_value,
                               fill = fill_value))+
    facet_nested(facet_y~Condition+PID+Stimulation,
                 scales = "free",
                 space = "free_x",
                 switch = "y",
                 strip = strip_nested(background_y = list(element_rect(fill = "grey90",  colour = "black")),
                                      background_x = list(element_rect(colour = "black", fill = "grey90"), 
                                                          element_part_rect(side = "b", colour = "black", fill = NA),                   
                                                          element_rect(colour = NA, fill = NA)),
                                      by_layer_y = T,
                                      by_layer_x = T))+
    geom_bar(stat="identity") +
    ylab(label = "Percent positive in population (%)")+
    xlab(label = "Timepoint")+
    scale_fill_manual(values = CATALYST:::.cluster_cols[1:length(unique(dataframe$fill_value))])+
    guides(fill=guide_legend(title=legend_title,
                             ncol = 1))+
    theme(legend.title=element_text(size = 14),
          plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size = 14,
                                   angle = 90,
                                   hjust = 1, 
                                   vjust = 0.5),       
          axis.text.y=element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 10),
          axis.title=element_text(size = 18, face="bold"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_rect(fill = 'grey95'),
          panel.spacing.x = unit(c(rep(c(rep(.5,n_stims-1),2),n_PIDs-1),rep(.5,n_stims-1)), "lines"),
          panel.spacing.y = unit(1, "lines"),
    )
  
  return(p)
  
}


