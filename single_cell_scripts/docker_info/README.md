## Reproducing the Analysis with Docker

Starting from a `rocker/rstudio` base container, we developed a Docker environment that reproduces the exact setup used for our single-cell analysis. The container enables users to fully reproduce all analyses or adapt the environment for their own datasets.

## Requirements

### System Requirements

The container can be run on a desktop computer using Docker Desktop. OS-specific system requirements are described in the Docker documentation: https://docs.docker.com/desktop/. Importantly, installing Docker Desktop requires Administrator privileges, but running it does not.

The analysis was run using Docker Desktop v28.4.0 on Windows 11 with 64 GB of RAM. While systems with lower memory may also work, at least 64 GB of RAM is recommended to ensure stable performance when running the entire pipeline, given the dataset size and the computational demands of some steps.

### Raw data

Available for download through the European Genome-Phenome Archive (EGA), accession number [to be updated].

Files should be named using the prefix `aim_[IDXXX]_[vXX]` where `IDXXX` denotes the donor identifier and `vXX` the visit number, followed by a suffix indicating the file type: 

-   `aim_[IDXXX]_[vXX]_raw_feature_bc_matrix.h5` for the raw count matrix
-   `aim_[IDXXX]_[vXX]_sample_filtered_feature_bc_matrix.h5` for the filtered count matrix
-   `aim_[IDXXX]_[vXX]_filtered_contig_annotations.csv` for the TCR contig annotations

All files should be placed in a single directory named `sc_data_by_donor`, located within a `main_directory/` that is set as the working directory.

If the `.rds` file has also been downloaded from the EGA, it should be placed in the `main_directory` and named `aim_final.rds`.

### Docker image

A Docker Desktop version compatible with your operating system can be downloaded from https://www.docker.com/products/docker-desktop/. OS-specific installation instructions are available at https://docs.docker.com/desktop/. After installation, open Docker Desktop to start the Docker engine.

The Docker image can be built locally from the Dockerfile in the project's root directory, or pulled directly from Docker Hub:

```{bash}

docker pull gsf70/fisher_garcia_frattari_naasz_26:published

```

## Further instructions

### Start the container

With `main_directory/` set as the current working directory, start the container from a local terminal:

```{bash}

docker run --name fisher_garcia_frattari_naasz -p 8787:8787 -e PASSWORD=Vajolet -e DEFAULT_USER=rstudio -v '[PATH_TO_DATA]:/fisher_garcia_frattari_naasz_sc_analysis/data' -d gsf70/fisher_garcia_frattari_naasz_26:published

```

This command launches a containerized RStudio Server instance locally. No internet connection is required after the image has been pulled.

The container name (`--name`), port (`-p`), and the RStudio Server username (`DEFAULT_USER`) and password (`PASSWORD`) can be customized, but all must be specified.

`[PATH_TO_DATA]` refers to the full local path of the working directory. For example, on Windows:

```{bash}

-v "C:\Users\FGFN\Projects\project_dir:/fisher_garcia_frattari_naasz_sc_analysis/data"

```

If the terminal is already in `main_directory/`, the relative path can be used:

```{bash}

-v ".:/fisher_garcia_frattari_naasz_sc_analysis/data"

```

This command mounts the local project directory to the container, enabling read and write access to all results directly on the host computer.

### Access RStudio Server

Open a web browser (e.g., Microsoft Edge or Firefox), navigate to `http://localhost:8787/` (adjust the port number if customized), and log in using the username and password specified above (username: `rstudio`, password: `Vajolet`).

### Set Up the Folder Structure

The script `R/utilities/structure_setup.R` recreates the subdirectory structure within the  `main_directory/` that was used to run the analysis pipeline. To execute it, open and run the script, or use the following command in the R console of RStudio:

```{r}

source("R/utilities/structure_setup.R")

```

### Ready to Start

At this stage, the Docker image should be running, and the local directory structure should match the one expected by the scripts.

## Instructions for Use

### Option 1: Run the Entire Pipeline

All steps from raw data to final figures can be reproduced using the scripts in the `R`directory. 

Begin with the scripts in `R/preprocessing/`directory, which are numbered in running order (e.g., `R/preprocessing/01_object_generation.R` is the first).

**Note:** 
To ensure full reproducibility and avoid environment conflicts, **empty the environment** and restart the **R session** between scripts. You can use the script `R/utilities/session_restart`to automate this:

```{r}

source("R/utilities/session_restart.R")

```

Running the complete pipeline is expected to take several hours on a standard desktop computer (approximately 2-3 hours if differential gene expression testing with MAST - performed at the end of `R/preprocessing/03_reclustering_and_annotation.R` - is skipped, otherwise longer). 

### Option 2: Inspect the Final Results Directly

If the `aim_final.rds` file has been downloaded, the scripts in `R/figures/` can be executed directly after the structure setup without running the full pipeline. To reproduce the figures, run the following commands in the R console:

```{r}

# Figures in the main text
source("R/figures/figures.R")

# Extended data and supplementary figures
source("R/figures/extended_and_supplementary_figures.R")

```

Running this option typically takes 5-10 minutes on a standard desktop computer. The expected output of the demo is a set of figures identical to those presented in the manuscript, saved in `processed_data/results/figures/`.