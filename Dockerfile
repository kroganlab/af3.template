FROM r-base:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/* \
    && R -q -e 'install.packages(c("bio3d", "data.table", "ggplot2", "ggrepel", "jsonlite", "magrittr", "rjson", "stringr", "usedist"), repos = "https://ftp.osuosl.org/pub/cran/")' \
    && R -q -e 'install.packages("BiocManager", repos = "https://ftp.osuosl.org/pub/cran/")' \
    && R -q -e 'BiocManager::install(c("Biostrings","ComplexHeatmap"))'
