[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# Tomato Seed Microbiome Analysis

This repository contains a comprehensive analysis of seed microbiome plasticity across 100 different tomato (*Solanum lycopersicum* L.) genotypes. The study explores how genotype variation influences seed microbiome composition and its relationship with host phenotypic traits.

## Table of Contents
- [Tomato Seed Microbiome Analysis](#tomato-seed-microbiome-analysis)
  - [Table of Contents](#table-of-contents)
  - [Overview](#overview)
  - [Key Features](#key-features)
  - [Data](#data)
  - [Project Structure](#project-structure)
  - [Analysis Pipeline](#analysis-pipeline)
    - [1. Data Preprocessing](#1-data-preprocessing)
    - [2. Alpha Diversity Analysis](#2-alpha-diversity-analysis)
    - [3. Beta Diversity Analysis](#3-beta-diversity-analysis)
    - [4. Machine Learning](#4-machine-learning)
    - [5. Microbiome-Phenotype Association](#5-microbiome-phenotype-association)
  - [Key Findings](#key-findings)
  - [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Installation and Setup](#installation-and-setup)
    - [Running the Analysis](#running-the-analysis)
    - [Reports](#reports)
  - [Dependencies](#dependencies)
  - [Contributing](#contributing)
  - [License](#license)
  - [Citation](#citation)
  - [Authors](#authors)
  - [Contact](#contact)

## Overview
The seed microbiome is crucial for plant health, development, and productivity. This project investigates the extent to which tomato genotype variation affects seed microbiome composition using a large-scale analysis of 100 diverse genotypes.

## Key Features
- **Alpha & Beta Diversity**: Comprehensive analysis of microbial diversity within and between samples.
- **Machine Learning**: Random Forest models to predict host genotype and phenotype from microbiome data.
- **Statistical Analysis**: Robust testing with PERMANOVA, pairwise comparisons, and dispersion analysis.
- **Reproducible Reporting**: A full analysis pipeline documented in a Quarto (`.qmd`) report for complete reproducibility.
- **Data Visualization**: PCoA, NMDS, heatmaps, and various plots to explore data patterns.

## Data
The primary data files used for this analysis are:
- **`ps2.RData`**: A pre-processed `phyloseq` object containing the OTU table, taxonomy, metadata, and phylogenetic tree. This is the main data object for most analyses.
- **`Metadata.csv`**: Comprehensive sample metadata, including genotype, origin, and phenotypic traits.
- **`feature-table.biom`**: The core OTU abundance data in BIOM format.
- **`Taxonomy.txt`**: Taxonomic assignments for each OTU.
- **`tree.nwk`**: The phylogenetic tree used for diversity calculations.

## Project Structure
The repository is organized to separate scripts, data, and results.
```
Tomato-seed-microbiome/
├── R/                                  # Helper R scripts and functions
├── data/                               # Core data files (OTU, taxonomy, metadata)
├── outputs/                            # Main output directory for figures and tables
│   ├── Alpha-Diversity/
│   ├── Beta-Diversity/
│   └── Random_Forest/
├── .gitignore                          # Specifies files for Git to ignore
├── .gitattributes                      # Configures Git LFS for large files
├── README.md                           # This file
├── Tomato-seed-microbiome.Rproj        # R project file
├── Tomato_seed_microbiome-reviews.qmd  # Main analysis (Quarto)
└── Tomato_seed_microbiome-reviews.html # Rendered HTML report
```

## Analysis Pipeline

### 1. Data Preprocessing
- Quality filtering and removal of contaminants.
- Taxonomic classification of microbial sequences.
- Rarefaction analysis to determine an appropriate sequencing depth.
- Prevalence filtering with a 0.25% threshold.

### 2. Alpha Diversity Analysis
- Rarefaction curves across multiple sequencing depths (500-100,000 reads)
- Diversity metrics: Observed OTUs, Shannon, Simpson, Inverse Simpson, Chao1, ACE
- Statistical comparisons with Dunn's test

### 3. Beta Diversity Analysis
- Multiple distance metrics: Bray-Curtis, Jaccard, UniFrac
- Normalization methods: Relative abundance, CSS, Hellinger transformation
- Ordination: PCoA and NMDS
- Statistical testing: PERMANOVA with pairwise comparisons
- Dispersion analysis: `betadisper` for homogeneity testing

### 4. Machine Learning
- Random Forest models for:
  - Genotype prediction (90% accuracy)
  - Production site classification (82% accuracy)
  - Phenotypic trait prediction (e.g., fruit shape: 84%)
- Feature importance analysis and confusion matrix evaluation.

### 5. Microbiome-Phenotype Association
- Genotype-specific microbiome signatures
- Regional variation analysis across Chinese provinces
- Correlation with phenotypic traits

## Key Findings
- **Genotype-Specific Microbiomes**: Tomato genotypes harbor distinct and predictable seed microbiome signatures.
- **Phenotype Correlation**: Seed microbiome composition correlates strongly with key phenotypic traits such as fruit color, shape, and taste.
- **Regional Influence**: Geographical location and environmental factors significantly influence microbiome assembly, though genotype effects remain robust.
- **Predictive Power**: Machine learning models can predict a tomato plant's genotype and key phenotypic traits with high accuracy based solely on its seed microbiome profile.

## Getting Started

### Prerequisites
- Git and Git LFS
- R version 4.0 or higher
- RStudio (Recommended)

### Installation and Setup
1. Clone the repository:
   ```bash
   git clone https://github.com/kaboyo/Tomato-seed-microbiome.git
   cd Tomato-seed-microbiome
   ```
2. Set up and pull large files with Git LFS:
   ```bash
   git lfs install
   git lfs pull
   ```

2. Open the R project file in RStudio:
   ```r
   # Open Tomato-seed-microbiome.Rproj
   ```

3. Install the required packages by running the script in the Dependencies section.

### Running the Analysis
The entire analysis can be reproduced by rendering the main Quarto document:
- Open `Tomato_seed_microbiome-reviews.qmd` in RStudio.
- Click the "Render" button.

### Reports
The final interactive report is available in `Tomato_seed_microbiome-reviews.html`.

## Dependencies

This project requires R (version 4.0 or higher). All required packages can be installed by running the following script in your R console:

```r
install_dependencies <- function() {
  # CRAN packages
  cran_packages <- c(
    "vegan", "psych", "ggplot2", "ggpubr", "lawstat", "lme4", "broom", "rcompanion", 
    "FSA", "explore", "gridExtra", "tidyverse", "phia", "AICcmodavg", "multcomp", 
    "gcookbook", "patchwork", "rmarkdown", "ggh4x", "ggrepel", "pals", "ggforce", 
    "data.table", "dendextend", "knitr", "tidyplots", "egg"
  )
  
  # Bioconductor packages
  bioc_packages <- c(
    "phyloseq", "genefilter", "impute", "DESeq2", "microbiome", "microbiomeutilities", 
    "metagenomeSeq", "ALDEx2", "SummarizedExperiment", "Biostrings", "ComplexHeatmap", 
    "DirichletMultinomial", "HMP", "MicrobiotaProcess", "qiime2R", "speedyseq", 
    "phyloseq.extended", "ggtree", "ggtreeExtra", "phangorn"
  )
  
  # GitHub packages
  github_packages <- c(
    "pmartinezarbizu/pairwiseAdonis",
    "david-barnett/microViz"
  )
  
  # Install CRAN packages
  new_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
  if (length(new_cran)) install.packages(new_cran)
  
  # Install Bioconductor packages
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
  if (length(new_bioc)) BiocManager::install(new_bioc)
  
  # Install GitHub packages
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github(github_packages, quiet = TRUE)
  
  print("All dependencies installed.")
}

install_dependencies()
```


## Contributing
Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this code or data in your research, please cite this repository.

## Authors
- **Expedito Olimi** - Primary researcher and analyst

## Contact
For questions or collaborations, please open an issue in this repository or contact Expedito Olimi.

---