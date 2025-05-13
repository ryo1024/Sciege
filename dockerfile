# Start from the official R base image
FROM rocker/r-ver:4.3.1

# Install system libraries and tools
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \ 
    liblzma-dev \
    zlib1g-dev \
    build-essential \
    git \
    wget \
    samtools \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and CRAN packages
RUN R -e "install.packages('BiocManager')"

RUN R -e "BiocManager::install(c('bambu'))"

RUN R -e "install.packages(c('devtools', 'testthat'))"

RUN R -e "install.packages('BisqueRNA')"

RUN R -e "install.packages('plyr')"

RUN R -e "install.packages(c('nnls', 'boot'))"

RUN wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar -xzf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin/ && \
    rm -r gffread-0.12.7.Linux_x86_64*

RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -xzf kallisto_linux-v0.46.1.tar.gz && \
    mv kallisto/kallisto /usr/local/bin/ && \
    rm -r kallisto*

RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm -r seqkit*

RUN wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    tar -xzf gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    mv gffcompare-0.12.6.Linux_x86_64/gffcompare /usr/local/bin && \
    rm -r gffcompare*

# Copy your package source code into the container
COPY . /Sciege

# Install your package (assumes proper DESCRIPTION, NAMESPACE, etc.)
RUN R -e "devtools::install('/Sciege')"

# Set default working directory
WORKDIR /Sciege

# Set default command: launch R console
CMD ["R"]