## Reproducing the Analysis with Docker

Starting from a `rocker/rstudio` base container, we developed a Docker environment that reproduces the exact setup used for our single-cell analysis. The container enables users to fully reproduce all analyses or adapt the environment for their own datasets.

## Requirements

### Raw data

Available for download through the European Genome-Phenome Archive (EGA), accession number [to be updated]. 

For each donor, create a directory named `aim_[IDXXX]_[vXX]` and containing:

-   a `raw_feature_bc_matrix`subdirectory
-   a `sample_filtered_feature_bc_matrix` subdirectory
-   a `filtered_contig_annotations` file

Collect all donor-specific directories into a single folder named `sc_data_by_donor`.

### Docker image

The Docker image can be built locally from the Dockerfile in the project's root directory, or pulled directly from Docker Hub:

```{bash}

docker pull gsf70/fisher_garcia_frattari_naasz_26:published

```

## Further instructions

### Start the container

Create a working directory for the project and set it as your working directory. Start the container from the local terminal using:

```{bash}

docker run --name fisher_garcia_frattari_naasz -p 8787:8787 -e PASSWORD=Vajolet -e DEFAULT_USER=rstudio -v '[PATH_TO_DATA]:/fisher_garcia_frattari_naasz_sc_analysis/data' -d gsf70/fisher_garcia_frattari_naasz_26:published

```

Once executed, this command launches a containerized RStudio Server instance locally. No internet connection is required.

The container name (`--name`), port (`-p`), and the RStudio Server username (`DEFAULT_USER`) and password (`PASSWORD`) can be customized, but all must be specified.

`[PATH_TO_DATA]` refers to the full local path of the working directory. For example:

```{bash}

-v "C:\Users\FGFN\Projects\project_dir:/fisher_garcia_frattari_naasz_sc_analysis/data"

```

Since your terminal’s current working directory is the project directory (as set above), you can simply use:

```{bash}

-v ".:/fisher_garcia_frattari_naasz_sc_analysis/data"

```

This command mounts the local project directory to the container, enabling read and write access to all results directly on the host computer.

### Access RStudio Server

Open a web browser (e.g., Microsoft Edge or Firefox), navigate to `http://localhost:8787/` (adjust the port number if customized), and log in using the username and password specified above (username: `rstudio`, password: `Vajolet`).

### Set Up the Folder Structure

The script `R/utilities/structure_setup.R` recreates the subdirectory structure within `Fisher_Garcia_Frattari_Naasz_et_al_2025` that was used to run the analysis pipeline. To execute it, open and run the script, or use the following command in the R console of RStudio:

```{r}

source("R/utilities/structure_setup.R")

```

After running this script, the directory `raw_data/single_cell_data/` will be created. Copy the `sc_data_by_donor` folder into this location.

### Ready to Start

At this stage, the Docker image should be running, and the local directory structure should match the one expected by the scripts.

## Instructions for Use

All steps from raw data to final figures can be reproduced using the scripts in the `R`directory. 

Begin with the scripts in `R/preprocessing`directory, which are numbered in running order (e.g., `R/preprocessing/01_object_generation.R` is the first).

**Note:** 
To ensure full reproducibility and avoid environment conflicts, **empty the environment** and restart the **R session** between scripts. You can use the script `R/utilities/session_restart`to automate this:

```{r}

source("R/utilities/session_restart.R")

```

Running the complete pipeline is expected to take several hours on a standard desktop computer (approximately 2-3 hours if differential gene expression testing with MAST - performed at the end of `R/preprocessing/03_reclustering_and_annotation.R` - is skipped, otherwise longer).
