#____________________________#
#### Source and libraries ####
#____________________________#

rm(list=ls())
source("0_source.R")

#____________________________________#
#### Load data and prepare gating ####
#____________________________________#

# read transformed flowSet
fs_transform = readRDS(file = glue("{dir_project}/2_pipeline/{folder_transform}/out/FlowSet_transform.rds"))

# create gatingSet from transformed flowSet
gs = GatingSet(fs_transform)

# make sure rownames and names column are the same (otherwise gating will not work)
cyto_details(gs)$name = rownames(cyto_details(gs))

# gatingTemplate
gatingTemplate_name = "gatingTemp_pregating.csv"
# cyto_gatingTemplate_apply(gs, gatingTemplate_name) (run if gates are already drawn)

#__________________#
#### Draw gates #### 
#__________________#
# (run if no gates have been drawn yet)

# Boundary
cyto_gate_edit(gs,
               select = 1,
               parent = "root",
               alias =  "boundary",
               channels = c("FSC-A", "SSC-A"),
               type = "boundary",
               axes_limits = "data",
               display = .8,
               gatingTemplate = gatingTemplate_name)

# Cells
cyto_gate_draw(gs,
               select = 1,
               parent = "boundary",
               alias = "Cells",
               channels =  c("FSC-A","SSC-A"),
               type = "polygon",
               axes_limits = "data",
               display = 50000,
               gatingTemplate = gatingTemplate_name)

recompute(gs)
cyto_plot(gs[sample(c(1:nrow(cyto_details(gs))),20)],parent = "boundary", alias = "", channels = c("FSC-A","SSC-A"),display = 50000, group_by = "sample_id")


# Live
cyto_gate_edit(gs,
               select = 1,
               parent = "Cells",
               alias = "Live",
               channels =  c("LiveDeadFixableBlue-A","SSC-A"),
               type = "polygon",
               axes_limits = "data",
               display = 50000,
               gatingTemplate = gatingTemplate_name)

recompute(gs)
cyto_plot(gs[sample(c(1:nrow(cyto_details(gs))),20)],parent = "Cells", alias = "", channels = c("LiveDeadFixableBlue-A","SSC-A"),display = 50000, group_by = "sample_id")


# Singlets
cyto_gate_draw(gs,
               select = 1,
               parent = "Live",
               alias = "Singlets",
               channels =  c("FSC-A","FSC-H"),
               type = "polygon",
               axes_limits = "data",
               display = 50000,
               gatingTemplate = gatingTemplate_name)

recompute(gs)
cyto_plot(gs[sample(c(1:nrow(cyto_details(gs))),10)],parent = "Live", alias = "", channels = c("FSC-A","FSC-H"),display = 50000, group_by = "sample_id")


# Lymphocytes
cyto_gate_edit(gs,
               select = 1,
               parent = "Singlets",
               alias = "Lymphocytes",
               channels =  c("FSC-A","SSC-A"),
               type = "polygon",
               axes_limits = "data",
               display = 50000,
               gatingTemplate = gatingTemplate_name)

recompute(gs)
cyto_plot(gs[sample(c(1:nrow(cyto_details(gs))),10)],parent = "Singlets", alias = "", channels = c("FSC-A","SSC-A"),display = 50000, group_by = "sample_id")

#__________________#
#### Plot gates ####
#__________________#

display_cells = 40000
mat <- matrix(c(1,2,3,4), ncol = 4, byrow = T)

# Call cyto_plot_save to save layout to png file
cyto_plot_save(glue("{dir_store}/SuppFig14.svg"),height = 4, width = 16, units = "in",res = 300)
#cyto_plot_save(glue("{dir_store}/SuppFig14.tiff"),height = 4, width = 16, units = "in",res = 300)
#cyto_plot_save(glue("{dir_store}/SuppFig14.png"),height = 4, width = 16, units = "in",res = 300)

# Set custom layout using cyto_plot_custom
cyto_plot_custom(mat)

# Cells
cyto_plot(gs[[1]],
          parent =  "boundary",
          alias = "",
          channels = c("FSC-A","SSC-A"),
          axes_limits = "data",
          label_text_x = 1000000,
          label_text_y = 1000000,
          label_text_size = .8,
          display = display_cells, 
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "root")

# Live
cyto_plot(gs[[1]],
          parent = "Cells",
          alias = "",
          channels =  c("LiveDeadFixableBlue-A","SSC-A"),
          axes_limits = "instrument",
          label_text_x = 0,
          label_text_y = 1600000,
          label_text_size = .8,
          gate_line_width = 1.5,
          gate_line_col = "black",
          display = display_cells,
          title = "Cells")


# Singlets
cyto_plot(gs[[1]],
          parent = "Live",
          alias = "",
          channels = c("FSC-A","FSC-H"),
          label_text_x = 200000,
          label_text_y = 300000,
          axes_limits = "data",
          label_text_size = .8,
          display = display_cells,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "Live")

# Lymphocytes
cyto_plot(gs[[1]],
          parent =  "Singlets",
          alias = "",
          channels = c("FSC-A","SSC-A"),
          axes_limits = "data",
          label_text_x = 600000,
          label_text_y = 1000000,
          label_text_size = .8,
          display = display_cells,
          gate_line_width = 1.5,
          gate_line_col = "black",
          title = "Singlets")

cyto_plot_complete()

#_________________________#
#### Extract viability ####
#_________________________#

df_viability = cyto_stats_compute(gs,
                                  alias = "Live",
                                  parent = "Cells",
                                  stat = "freq")
df_viability = df_viability[,c("sample_id","batch","PID","Visit","Frequency")]
colnames(df_viability) = c("sample_id","batch", "PID", "Visit","Viability")

#_______________________#
#### Save gated data ####
#_______________________#

end_pop = "Lymphocytes"

fs_pregated = cytoset_to_flowSet(gs_pop_get_data(gs,end_pop))
fs_pregated = q_PrefixSuffixSampleNames(fs_pregated, add = "suffix", glue("_lymph_pregate"))

# edit name column
pData(fs_pregated)$name = rownames(pData(fs_pregated))

# add viability
fs_pregated = q_MergeMetadataFlowSet(fs = fs_pregated,
                                     dataframe = df_viability[,c("sample_id","Viability")],
                                     merge_col = "sample_id")

# change col order pData
pData(fs_pregated) = pData(fs_pregated)[,c("name","sample_id","sample_id2","well_id","date","batch","Trial","Condition","PID","Visit","Stimulation","Viability")]

# save RDS
saveRDS(fs_pregated, file = glue("{dir_out}/FlowSet_pregated_lymph.rds"))


