
FROM rocker/rstudio:4.4.0

# Install system dependencies ##################################################

RUN apt-get update && apt-get install -y \
    python3 python3-venv python3-pip \
    wget curl bzip2 \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgit2-dev \
    liblapack-dev \
    libblas-dev \
    libpng-dev \
    libjpeg-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libcairo2-dev \
    libglpk-dev \
    libglpk40 \
    make \
    g++ \
    pkg-config \
    libhdf5-dev \
    gawk \
    && rm -rf /var/lib/apt/lists/*

# Create project directory and make it writable ################################

RUN mkdir -p /fisher_garcia_frattari_naasz_sc_analysis \
 && chown -R rstudio:rstudio /fisher_garcia_frattari_naasz_sc_analysis \
 && chmod -R a+rwX /fisher_garcia_frattari_naasz_sc_analysis

WORKDIR /fisher_garcia_frattari_naasz_sc_analysis

# Restore R packages ###########################################################

COPY renv.lock renv.lock
RUN R -e "install.packages('renv'); renv::restore()"

# Add hdf5r to read h5 files
RUN R -e "renv::install('hdf5r')"

# Python environment setup #####################################################

# Copy requirements into the container
COPY python/python_env/pipy_requirements.txt /tmp/pipy_requirements.txt

# Create a venv and install Python dependencies
RUN python3 -m venv /opt/pipy_env \
 && . /opt/pipy_env/bin/activate \
 && pip install --upgrade pip wheel setuptools \
 && pip install -r /tmp/pipy_requirements.txt

# Copy analysis code ###########################################################

COPY . /fisher_garcia_frattari_naasz_sc_analysis

# After COPY, re-assert ownership/permissions (COPY runs as root) ##############
RUN chown -R rstudio:rstudio /fisher_garcia_frattari_naasz_sc_analysis \
 && chmod -R a+rwX /fisher_garcia_frattari_naasz_sc_analysis

# Build .Rprofile and .Renviron files ##########################################

# Add parameters to .Renviron file
RUN echo 'PYTHON_ENV="/opt/pipy_env"' >> /fisher_garcia_frattari_naasz_sc_analysis/.Renviron
RUN echo 'PROCESSED_DATA="data/processed_data/"' >> /fisher_garcia_frattari_naasz_sc_analysis/.Renviron
RUN echo 'RAW_DATA="data/raw_data/single_cell_data/sc_data_by_donor/"' >> /fisher_garcia_frattari_naasz_sc_analysis/.Renviron
RUN echo 'R_LIBS_USER=/home/rstudio/R/x86_64-pc-linux-gnu-library/4.4' >> /fisher_garcia_frattari_naasz_sc_analysis/.Renviron

# Add commands to .Rprofile
RUN echo 'raw_data <- Sys.getenv("RAW_DATA")' >> /fisher_garcia_frattari_naasz_sc_analysis/.Rprofile
RUN echo 'processed_data <- Sys.getenv("PROCESSED_DATA")' >> /fisher_garcia_frattari_naasz_sc_analysis/.Rprofile
RUN echo 'set.seed(123)' >> /fisher_garcia_frattari_naasz_sc_analysis/.Rprofile

# Set default working directory for RStudio Server + remove timeout logout
RUN echo "session-default-working-dir=/fisher_garcia_frattari_naasz_sc_analysis" >> /etc/rstudio/rsession.conf

# Create writeable .Rhistory file (necessary for repl_python) ##################

RUN touch /fisher_garcia_frattari_naasz_sc_analysis/.Rhistory \
 && chown rstudio:rstudio /fisher_garcia_frattari_naasz_sc_analysis/.Rhistory \
 && chmod 664 /fisher_garcia_frattari_naasz_sc_analysis/.Rhistory
