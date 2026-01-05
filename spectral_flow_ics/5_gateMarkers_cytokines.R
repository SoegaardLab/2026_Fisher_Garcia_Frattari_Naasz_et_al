#___________________#
#### Load source ####
#___________________#

rm(list=ls())
source("0_source.R")

#____________________________________#
#### Load data and prepare gating ####
#____________________________________#

# Panel
panel = read.csv(file=glue("{dir_project}/2_pipeline/{folder_transform}/out/Panel.csv"))[,-1]

# Clean flowSet
fs_norm = readRDS(file = glue("{dir_project}/2_pipeline/4_cyCombine/out/FlowSet_norm_new.rds"))

# Create gatingSet
gs = GatingSet(fs_norm)

# Split samples for IFN gating
gs_rest = gs[-which(cyto_details(gs)$PID == "112"&cyto_details(gs)$Visit == "1")]
gs_112v1 = gs[which(cyto_details(gs)$PID == "112"&cyto_details(gs)$Visit == "1")]

my_gatingTemplate_rest = "gatingTemplate_cytokine+rest.csv"
my_gatingTemplate_112v1 = "gatingTemplate_cytokine+112v1.csv"

cyto_gatingTemplate_apply(gs_rest, my_gatingTemplate_rest)
cyto_gatingTemplate_apply(gs_112v1, my_gatingTemplate_112v1)


#__________________#
#### Draw gates ####
#__________________#

channel_CD8 = panel[(!is.na(panel$antigen) & panel$antigen=="CD8"),"fcs_colname"]
channel_IFN = panel[(!is.na(panel$antigen) & panel$antigen=="IFN"),"fcs_colname"]
channel_TNFa = panel[(!is.na(panel$antigen) & panel$antigen=="TNFa"),"fcs_colname"]
channel_IL2 = panel[(!is.na(panel$antigen) & panel$antigen=="IL2"),"fcs_colname"]
channel_CD107a = panel[(!is.na(panel$antigen) & panel$antigen=="CD107a"),"fcs_colname"]

#
# cyto_gate_draw(gs,
#                parent = "root",
#                select = which(cyto_details(gs)$IFN_gatingGroup=="1"),
#                overlay = gs_rest,
#                point_col = c("red","blue"),
#                channels = c("BV421-A","SparkBlue-574-A"),
#                alias = "IFN+",
#                axes_limits = "instrument",
#                type = "rectangle",
#                display = 10000,
#                merge_by = "IFN_gatingGroup",
#                gatingTemplate = my_gatingTemplate)
# cyto_gate_edit(gs_112v1,
#                select = list("Stimulation" == "Gag"),
#                parent = "root",
#                channels = c(channel_IL2,"SparkBlue-574-A"),
#                alias = "IL2+",
#                axes_limits = "instrument",
#                type = "rectangle",
#                display = 50000,
#                gatingTemplate = my_gatingTemplate_112v1)
# cyto_gate_edit(gs_rest,
#                select = list("PID" == "142",
#                              "Stimulation"=="Gag"),
#                parent = "root",
#                channels = c(channel_IL2,"SparkBlue-574-A"),
#                alias = "IL2+",
#                axes_limits = "instrument",
#                type = "rectangle",
#                display = 50000,
#                gatingTemplate = my_gatingTemplate_rest)

#_________________________________________#
#### Reorder samples to original order ####
#_________________________________________#

# merge the two gatingSets to one
gs_merge = merge_list_to_gs(list(gs_rest,gs_112v1))

# get back original order of samples
df = cyto_details(gs_merge)[,"sample_id",drop=F]
df$index = c(1:nrow(df))
df$sample_id = factor(df$sample_id, levels = pData(fs_norm)$sample_id)
order_original = df[order(df$sample_id),"index"]
gs_final = gs_merge[order_original]

# check sampleNames
all(sampleNames(gs_final)==rownames(pData(fs_norm)))

cyto_details(gs_final)$name 
rownames(cyto_details(gs_final))
cyto_details(gs_final)

#__________________#
#### Plot gates ####
#__________________#

sample_CD107a = data.frame(cyto_stats_compute(gs,
                   alias = "CD107a+",
                   parent = "root",
                   stat = "freq")) %>% arrange(desc(Frequency))%>%slice(1) %>%pull(sample_id)

sample_IFN = data.frame(cyto_stats_compute(gs,
                                              alias = "IFN+",
                                              parent = "root",
                                              stat = "freq")) %>% arrange(desc(Frequency))%>%slice(5) %>%pull(sample_id)

sample_TNF = data.frame(cyto_stats_compute(gs,
                                           alias = "TNFa+",
                                           parent = "root",
                                           stat = "freq")) %>% arrange(desc(Frequency))%>%slice(1) %>%pull(sample_id)

sample_IL2 = data.frame(cyto_stats_compute(gs,
                                           alias = "IL2+",
                                           parent = "root",
                                           stat = "freq")) %>% arrange(desc(Frequency))%>%slice(1) %>%pull(sample_id)

# Create layout matrix (number indicate plot no.)
mat <- matrix(c(1,2,3,4), ncol = 2, byrow = T)

# Call cyto_plot_save to save layout to png file
#cyto_plot_save(glue("{dir_store}/GatingStrategy2.png"),height = 10, width = 10, units = "in",res = 300)
#cyto_plot_save(glue("{dir_store}/SuppFig17.png"),height = 10, width = 10, units = "in",res = 300)
#cyto_plot_save(glue("{dir_store}/SuppFig17.svg"),height = 10, width = 10, units = "in",res = 300)
cyto_plot_save(glue("{dir_store}/SuppFig17.tiff"),height = 10, width = 10, units = "in",res = 300)

# Set custom layout using cyto_plot_custom
cyto_plot_custom(mat)

cyto_plot(gs[which(cyto_details(gs)$sample_id == sample_CD107a)],
          parent = "root",
          alias = "",
          channels = c(panel[(panel$marker_class=="state"&panel$antigen=="CD107a"),"fcs_colname"],"SparkBlue-574-A"),
          label_text_y = 4.2,
          axes_label_text_size = 1.5,
          point_size = 1.5,
          xlim = c(-1,5),
          ylim = c(-1,4.5),
          display = 150000,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "A                                                                                  ",
          title_text_size = 1.6
)
cyto_plot(gs[which(cyto_details(gs)$sample_id == sample_IFN)],
          parent = "root",
               alias = "",
               channels = c(panel[(panel$marker_class=="state"&panel$antigen=="IFN"),"fcs_colname"],"SparkBlue-574-A"),
               label_text_y = 4.2,
               axes_label_text_size = 1.5,
               point_size = 1.5,
               xlim = c(-1,5),
               ylim = c(-1,4.5),
               display = 150000,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "B                                                                                  ",
          title_text_size = 1.6
)
cyto_plot(gs[which(cyto_details(gs)$sample_id == sample_TNF)],
          parent = "root",
               alias = "",
               channels = c(panel[(panel$marker_class=="state"&panel$antigen=="TNFa"),"fcs_colname"],"SparkBlue-574-A"),
               label_text_y = 4.2,
               axes_label_text_size = 1.5,
               point_size = 1.5,
               xlim = c(-1,5),
               ylim = c(-1,4.5),
               display = 150000,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "C                                                                                  ",
          title_text_size = 1.6
)
cyto_plot(gs[which(cyto_details(gs)$sample_id == sample_IL2)],
          parent = "root",
          alias = "",
          channels = c(panel[(panel$marker_class=="state"&panel$antigen=="IL2"),"fcs_colname"],"SparkBlue-574-A"),
          label_text_y = 4.2,
          axes_label_text_size = 1.5,
          point_size = 1.5,
          xlim = c(-1,5),
          ylim = c(-1,4.5),
          display = 150000,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "D                                                                                  ",
          title_text_size = 1.6
)

# Call cyto_plot_complete once the figure is ready to be saved
cyto_plot_complete()


#___________________#
#### Save matrix ####
#___________________#

boolean_matrix = q_Gating_GetBoolMatrix(flowSet = fs_norm, gatingSet = gs_final,gatenames = c("IFN+","TNFa+","IL2+","CD107a+"))
saveRDS(boolean_matrix, file = glue("{dir_out}/boolean_cytokines.rds"))


