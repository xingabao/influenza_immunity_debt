# Base image: RStudio Server with R 4.5.2
FROM rocker/rstudio:4.5.2

# Maintainer
LABEL maintainer="Abao Xing <albertxn@126.com>"

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages from CRAN
RUN install2.r --error --skipinstalled \
    "glue" "dplyr" "tidyr" "data.table" "stringr" "forcats" \
    "ggplot2" "patchwork" "cowplot" "scales" "ggtext" "RColorBrewer" \
    "ggrepel" "ggh4x" "ggsci" "gstat" "ggspatial" "ggalluvial" \
    "lubridate" "zoo" "ISOweek" \
    "dlnm" "splines" "mgcv" "MASS" "lmtest" "mem" \
    "bootnet" "qgraph" "networktools" "NetworkComparisonTest" "RTransferEntropy" "rEDM" \
    "future" "openxlsx" "readxl" "writexl" "sf" "this.path" \
    && rm -rf /tmp/downloaded_packages

# Install R packages
ARG GITHUB_PAT
ENV GITHUB_PAT=${GITHUB_PAT}
RUN install2.r --error --skipinstalled "remotes" && rm -rf /tmp/downloaded_packages \
	&& Rscript -e 'remotes::install_github("hrbrmstr/ggchicklet")' \
	&& Rscript -e 'remotes::install_git("https://codeberg.org/hrbrmstr/hrbrthemes.git")'
	
# Install Bioconductor packages
RUN Rscript -e "install.packages('BiocManager'); \
    BiocManager::install(c('treeio', 'ggtree'), ask=FALSE, update=FALSE)"

# Copy local directories into the Docker image
COPY Scripts/ /home/rstudio/Scripts/
COPY data/ /home/rstudio/data/
COPY Figures/ /home/rstudio/Figures/
COPY Tables/ /home/rstudio/Tables/

# Fix permissions (Critical for RStudio user)
RUN chown -R rstudio:rstudio /home/rstudio

# Expose port
EXPOSE 8787