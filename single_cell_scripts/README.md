# Single cell analysis

Scripts used to analyze single-cell sequencing data in the Fisher, Garcia, Frattari, Naasz, `et al` 2026 manuscript.

## Contents

### R

Contains R-based scripts, organized into three subdirectories:

-   **utilities:** common functions used for data processing steps (`protein_coding_genes_extraction.R`, `qc_ambient_rna.R`, `object_processing_functions.R`), figure generation (`figure_utilities`), and Docker-based analysis (`structure_setup.R`, `session_restart.R`)

-   **preprocessing:** scripts for Seurat object generation, processing and analysis

-   **figures:** scripts used to generate the single-cell data figures

### Python
Contains Python-based scripts and information required to recreate the Python environment used in the analysis.

### Additional data

Additional data used in the scripts, included for reproducibility:

-   `protein_coding_genes.rds`: list of protein-coding genes used by the scripts.  The list was generated using the script `R/utilities/ protein_coding_genes_extraction.R` starting from a GTF with human genes (Gencode v32) from the 10x Genomics reference `refdata-gex-GRCh38-2020-A`. The reference can be downloaded from: https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads and should be placed in `additional_data/`

-   `totalseq_isotypes.txt`: Extract from the ADT reference provided in Supplementary File 1, containing marker names and corresponding isotypes

-   `visit_names.rds`: List of standardized visit names across samples

### renv

Files used by `renv` to recreate the same R environment.

### docker_info

Instructions for setting up and running the analysis within the project's Docker container. The Dockerfile used to build the image can be found in the home directory and may also serve as a template for future builds.

## Contact

If you have any questions, would like further details, or wish to discuss any aspects of the analysis, feel free to start an issue or reach out at gs.frattari\@clin.au.dk
