################################################################################
## Title        : Pipy env create
##
## Input        : Pipy requirements text file with packages needed. Path to
##                python installation (Python 3.11 here)
##
## Output       : Pipy venv to be used for clustering and scanpro analysis 
##
## Author       : Giacomo S Frattari
################################################################################

library(reticulate)

virtualenv_create(envname = "python/python_env/pipy_env",
                  requirements = "python/python_env/pipy_requirements.txt",
                  python = Sys.getenv("PYTHON_PATH")) # Retrieve python path from .Renviron file - alternatively,
                                                      # provide the path to python here

# Check instalation
reticulate::use_virtualenv("python/python_env/pipy_env", required = T)

py_module_available("leidenalg")
py_module_available("scanpro")
