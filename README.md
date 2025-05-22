# Sciege <img src="https://img.shields.io/badge/R-%23276DC3?style=flat&logo=r&logoColor=white" align="right" height="25"/>

**Sciege** is an R package for cell type level isoform expression analysis in bulk RNA-seq, intergrating with long-read RNA-seq and scRNA-seq. It leverages transcript discovery from long-read RNA-seq (supported by Bambu) and cell type decomposition using scRNA-seq as a reference (supported by Bisque). Using bulk short-read RNA-seq, we provide robust population mean estimate for cell type isoform abundance. Sciege also supports differential testing with stratified sample groups with bootstrapped estimates. Preprint analyzing GTEx and ROSMAP samples coming soon! 

---

## 🔧 Features

- End-to-end or modular cell type isoform expression pipeline
- Compatible with FASTQ, BAM, GTF, and external cell reference data
- Differential isoform expression at the cell-type level using bootstrapped statistics
- Docker support for reproducibility

---
## 🔧 System requirements

- Linux system (Tested on Ubuntu 20.0.4)
- Windows system with Docker enviromnment (Tested on Windows 10)

---
## 🔧 Software requirements

Install following softwares and make sure the binary files are in the environmental path.

- kallisto (Tested on v0.46.1, https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz)
- gffread (Tested on v0.12.7, https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz)
- seqkit (Tested on v2.1.0, https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz)
- gffcompare (Tested on v0.12.6, http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz)

For example,
```bash
wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz && \
    tar -xzf gffread-0.12.7.Linux_x86_64.tar.gz && \
    mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin/ && \
    rm -r gffread-0.12.7.Linux_x86_64*

wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -xzf kallisto_linux-v0.46.1.tar.gz && \
    mv kallisto/kallisto /usr/local/bin/ && \
    rm -r kallisto*

wget https://github.com/shenwei356/seqkit/releases/download/v2.10.0/seqkit_linux_amd64.tar.gz && \
    tar -xzf seqkit_linux_amd64.tar.gz && \
    mv seqkit /usr/local/bin/ && \
    rm -r seqkit*

wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    tar -xzf gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    mv gffcompare-0.12.6.Linux_x86_64/gffcompare /usr/local/bin && \
    rm -r gffcompare*
```
---

## 📦 Installation

Install from GitHub using [devtools](https://cran.r-project.org/package=devtools):

```r
# install.packages("devtools")
devtools::install_github("ryo1024/Sciege")
```

Or build locally:

```bash
git clone https://github.com/ryo1024/Sciege.git
cd Sciege
R CMD build .
R CMD INSTALL Sciege_*.tar.gz
```

If you run into installation issues we strongly recommend using docker environment to avoid any depedency conflicts in R. In normal R environment, installation should take less than 10 minutes.

---

## 🐳 Docker Support


Please refer to official docker website for installation (https://www.docker.com/). 
In your Sciege directory, build the Docker image (typically takes around 15 mins):

```bash
docker build -t sciege .
```

Run the pipeline (adjust file paths as needed):

```bash
docker run -v $PWD:/data sciege Rscript -e 'Sciege::run_Sciege(...)'
```

---

## 🚀 Example Usage

```r
library(Sciege)

# Run the full pipeline
run_pipeline(
  fastq_dir = "fastqs/",
  gtf_file = "annotations.gtf",
  genome_fasta = "genome.fa",
  ref_scrna = "scRNA_reference.rds",
  metadata_file = "sample_metadata.csv",
  out_dir = "sciege_output/"
)

# Or run individual steps
ref <- create_transcript_reference(bam_file = "input.bam", gtf_file = "annotations.gtf", genome_fasta = "genome.fa")
quant <- quantify_isoforms_kallisto(fastq_dir = "fastqs/", index = "transcript_index", output_dir = "quant_output/")
decomp <- decompose_cell_types_bisque(expr_matrix = quant, reference_sce = readRDS("ref_scrna.rds"), out_dir = "decomp_output/")
estimates <- estimate_isoform_expression_by_celltype(quant, decomp)
results <- differential_isoform_expression(quant, decomp, metadata, group_col = "condition")
```

---

## 📂 Package Structure

```
Sciege/
├── R/                      # Function scripts
├── man/                    # Manual (.Rd) files
├── tests/                 # Unit tests (testthat)
├── vignettes/             # Usage documentation (optional)
├── Dockerfile             # Docker support
├── DESCRIPTION            # Package metadata
├── NAMESPACE              # Function exports
└── README.md              # This file
```

---

## 🧪 Running Demos

To see example run of the pipeline with small sub set of files run following commands. Typical run of all tests should take less than 5 mins.
From R:

```r
devtools::test()
```

Or from command line (in Docker):

```bash
docker run -v $PWD:/pkg sciege Rscript -e 'devtools::test("/pkg")'
```

If test runs successfully, following results should appear at the end of stdout.

```
══ Results ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════Duration: 17.4 s

[ FAIL 0 | WARN 3 | SKIP 0 | PASS 25 ]
```

---

## 📃 License

MIT License

---

## 👥 Contributors

- Ryo Yamamoto (@ryo1024) UCLA Bioinformatics IDP

---

## 📫 Contact

If you encounter issues or have feature suggestions, feel free to open an [Issue](https://github.com/ryo1024/Sciege/issues), submit a Pull Request or shoot me an email at ryo10244201@g.ucla.edu.
