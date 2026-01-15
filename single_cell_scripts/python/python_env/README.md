### Pipy: Python environment used by single-cell scripts

The single-cell analysis requires a dedicated Python environment:

-   `FindClusters` uses the Python package **leidenalg**.
-   `scanpro` is implemented as a Python module.

To support these, a Python virtual environment `(venv)` was created with the necessary packages and dependencies.

### Recreating the environment

-   **Using Docker:**

    The environment is automatically created at `/opt/pipy_env` during image build.

    -   This path is stored in `.Renviron` as `PYTHON_ENV`.
    -   Scripts call it through `Sys.getenv("PYTHON_ENV")`.
    -   **No further action is needed** when running the pipeline inside Docker.

-   **Using R outside Docker**:

    1.  Use `python_env/pipy_requirements.txt` together with `python_env/pipy_env_creation.R` to build the environment.
    2.  Add the environment path to `.Renviron` as `PYTHON_ENV` (so scripts can access it via `Sys.getenv("PYTHON_ENV")`.
        -   Alternatively, provide the full path directly in the scripts instead of referencing `PYTHON_ENV`.

### Running scanpro with reticulate

-   The script `python/cluster_prop_aim.py` produces several objects in Python's main module.

-   These results are accessed in R through **reticulate** using the syntax:

    ``` r
    py$<object_name>
    ```

    where `<object_name>` corresponds to the Python variable to be retrieved
