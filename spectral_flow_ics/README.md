# Spectral flow cytometry-based intracellular cytokine staining

This document outlines a detailed data analysis pipeline for spectral flow cytometry-based intracellular cytokine staining. The workflow draws inspiration from the following key resources:

-   *den Braanker H, Bongenaar M and Lubberts E (2021) ["How to Prepare Spectral Flow Cytometry Datasets for High Dimensional Data Analysis: A Practical Workflow."](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.768113/full) Front. Immunol. 12:768113. doi: 10.3389/fimmu.2021.768113*

-   *Van Gassen S, Callebaut B, Van Helden M, Lambrecht B, Demeester P, Dhaene T, Saeys Y (2015). [“FlowSOM: Using self-organizing maps for visualization and interpretation of cytometry data.”](<https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.22625>) Cytometry Part A, 87(7), 636-645.*

Additional guidance on best practices for high-dimensional cytometry data analysis can be found [here](https://www.nature.com/articles/s41590-021-01006-z). 

### Source code

The file [0_source](/spectral_flow_ics/0_source.R), contains custom wrapper functions used throughout the pipeline. Each data analysis step and its corresponding .R or .Rmd file is described below.

### Step 1: Data transformation

**File to use**: [1_transform](/spectral_flow_ics/1_transform.Rmd)

As the physical process of exciting, emitting and detecting fluorescence signals on a flow cytometer results in increased variance of fluorescence signal when mean fluorescence intensity increases (signal variance is inhomogenous), transformation of data is needed to stabilise the variance. 

Data is transformed using *logicle transformation* from the **flowWorkspace** package by using the wrapper function `q_TransformEstimateLogicle` using a reference sample of your choosing. If no reference sample is specified, the first sample of the dataset will be used. Ideally, the reference file should include both a negative and a positive population for each marker, though this can be challenging for markers that show continuous and low expression patterns, for example, exhaustion markers such as
PD-1. 

### Step 2: Pregating

**File to use**: [2_pregate](/spectral_flow_ics/2_pregate.R)

The second step of data analysis is pregating to exclude debris, doublets, and dead cells. Gating is done using the **CytoExploreR** package with interactive pop-up windows to draw you gates. CytoExploreR needs a **gatingSet** object as input, which can be created directly from the transformed flowSet. Drawn gates are saved in a **gatingTemplate** (a csv file containing information on parent population, gate name, gate coordinates, etc.) which can be reused on other samples.

After gating, the end gate population is saved as a flowSet and used for downstream analysis.

### Step 3: Automated quality control

**File to use**: [3_qualitycontrol](/spectral_flow_ics/3_qualitycontrol.Rmd)

Automatic quality control is performed using the **PeacoQC** package for removal of low-quality events, e.g.

-   Temporary shift in signal (clogs/slow uptake)
-   Permanent shift in signal (speed/flow rate change)
-   Monotonic decrease in signal (contamination)
-   Monotonic increase in signal (contamination)

For each file, a report is created. The report shows events sorted by acquisition time and the signal for each marker. Removed events are highlighted in color codes based on which part of the algorithm decided to remove it. The PeacoQC algorithm creates QC fcs files, which can be read back into RStudio as a flowSet for downstream analysis.

### Step 4: Batch correction/normalization

**File to use**: [4_normalize](/spectral_flow_ics/4_normalize.Rmd)

After removal of low quality events, batch correction is performed using **cyCombine**. This is essential for ensuring downstream clustering and visualization reflects biological differences rather than technical noise. A replicated sample is used as anchor, and performance of the batch corrections is evaluated by EMD and MAD, plotting the UMAP before and after batch correction colored by batch label, and lastly a visual inspection of marker distribution in the replicate sample before and after batch correction. 

### Step 5: Gate cytokine positivity

**File to use**: [5_gateMarkers_cytokines](/spectral_flow_ics/5_gateMarkers_cytokines.R)

Like with pregating, gates for cytokine positivity are drawn using CytoExploreR. A gatingSet is again created from the QC flowSet and events are gated for positivity of the cytokines IFN, TNF, and IL2 as well as the degranulation marker CD107a.

Using the function `q_Gating_matrix_vector` we can get a **boolean matrix of TRUE/FALSE for marker positivity** for all cells in the flowSet/gatingSet. This matrix can later be added to the data for visualization and quantification of polyfunctional T cells.

### Step 6a: Clustering 1 – Major Lymphocyte Populations

**File to use**: [6a_clustering1](/spectral_flow_ics/6a_clustering1.Rmd)

To cluster the normalized data, we first need to convert it to a Single Cell Experiment (SCE). This is done with the `prepData` function from **CATALYST**. Next, FlowSOM clustering is performed with the canonical lineage markers CD3, CD4, CD8, CD19, CD19 and a 10x10 grid size. A maximum of 15 metaclusters are found via ConsensusPlus metaclustering. Everything is done in one step with the `cluster` wrapper function from CATALYST.
In addition to FlowSOM clustering, a UMAP is built using the same markers as for clustering and downsampling the data to 10,000 cells per sample. 
Cells are annotated based on the expression of the markers used for clustering by inspecting 2D scatterplots, UMAP overlayed with marker expression, heatmaps, etc.

### Step 6b: Clustering 2 – CD4 and CD8 T Cell Subset subclustering

**File to use**: [6b_clustering2](/spectral_flow_ics/6b_clustering2.Rmd)

To subcluster CD4 and CD8 T cells, cells from step 6a are filtered by their cell annotation labels and the markers to use for clustering are determined by the top 10 variable markers in each dataset. 

Clustering and building the UMAP is performed as described above. 

### Step 7: Determine memory CD4 and CD8 T cell polyfunctionality

**File to use**: [7_Polyfunctionality](/spectral_flow_ics/7_Polyfunctionality.Rmd)

Finally, T cell polyfunctionality is determined. This is done in 3 steps:

1. The cytokine positivity matrix created in step 5 is added to the filtered CD4 and CD8 data sets using the in-house `q_SCEAddBoolMatrix` function. 
2. All possible combinations of the markers are then found using the in-house `q_SCEBooleanCombinations` function. 
3. The percentage of cells positive for each of the combinations is computed and background signal from the negative stimulation sample is subtracted using the in-house `q_PsuedobulkPercPosMarker`.

If you have any questions, please do not hesitate to contact me at lisdie@clin.au.dk. 
