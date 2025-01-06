# Dockerfile for Seurat 4.3.0
FROM rocker/r-ver:4.3.0

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    libclang-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    python3-venv \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    liblzma5 \
    liblzma-dev \
    pkg-config

RUN apt-get install -y llvm-14

# Install system library for rgeos
RUN apt-get install -y libgeos-dev

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-14/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn
RUN pip3 install leidenalg

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

#RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix', version='1.6.0')"
RUN R --no-echo --no-restore --no-save -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz', upgrade = FALSE)"
# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager', version = '3.18', upgrade = FALSE)" 
# RUN R --no-echo --no-restore --no-save -e "BiocManager::install(version = '3.18')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install('Rhtslib', upgrade = FALSE)"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'), upgrade = FALSE)"

# Install CRAN suggests

RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools'), upgrade = FALSE)"

# Install spatstat
RUN R --no-echo --no-restore --no-save -e "install.packages(c('spatstat.explore', 'spatstat.geom'), upgrade = FALSE)"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r', upgrade = FALSE)"

# Install latest Matrix
#RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix', version='1.6.1')"

# Install rgeos
RUN R --no-echo --no-restore --no-save -e "install.packages('rgeos', upgrade = FALSE)"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes', upgrade = FALSE)"
#RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat',version='4.3.0')"
RUN R --no-echo --no-restore --no-save -e "remotes::install_version('SeuratObject', '4.1.3', repos = c('https://satijalab.r-universe.dev', getOption('repos')),upgrade=FALSE)"
RUN R --no-echo --no-restore --no-save -e "remotes::install_version('Seurat', '4.3.0', repos = c('https://satijalab.r-universe.dev', getOption('repos')),upgrade=FALSE)"
RUN R --no-echo --no-restore --no-save -e "install.packages('dplyr',upgrade = FALSE)"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk', upgrade = FALSE)"
RUN /usr/bin/python3.10 -m venv /root/.virtualenvs/r-reticulate
RUN /root/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade pip wheel setuptools
RUN /root/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user numpy
RUN /root/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user umap-learn
RUN /root/.virtualenvs/r-reticulate/bin/python -m pip install --upgrade --no-user leidenalg
# Clone GitHub repository
#RUN git clone https://github.com/GSbioinfo/HSPC_SingleCellAnalysis /scratch/CD34cells_Palantir/
# Set working directory
WORKDIR /scratch/CD34cells_Palantir/

# Default command
CMD ["R"]

#docker build --no-cache -t interactive_seurat430_env:1.0 --file Dockerfile . --progress=plain &> palantir_out.log
#docker save -o /scratch/DockerImages/interactive_seurat430_env.tar interactive_seurat430_env:1.0

#docker run -it -v /scratch/CD34cells_Palantir/:/scratch/CD34cells_Palantir/ -v /Github_repos/HSPC_SingleCellAnalysis/:/scratch/CD34cells_Palantir/scripts interactive_seurat430_env:1.0