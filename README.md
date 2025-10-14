# Tomato Seed Microbiome Analysis

This repository contains a comprehensive analysis of seed microbiome plasticity across 100 different tomato (*Solanum lycopersicum* L.) genotypes. The study explores how genotype variation influences seed microbiome composition and its relationship with host phenotypic traits.

## Table of Contents
- [Overview](#overview)
- [Data](#data)
- [Analysis Pipeline](#analysis-pipeline)
- [Key Findings](#key-findings)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Project Structure](#project-structure)
- [Contributing](#contributing)
- [License](#license)

## Overview
The seed microbiome is crucial for plant health, development, and productivity. This project investigates the extent to which tomato genotype variation affects seed microbiome composition using a large-scale analysis of 100 diverse genotypes.

The analysis encompasses:
- **Alpha diversity**: Richness, Shannon index, Simpson index, evenness, and phylogenetic diversity across different sequencing depths
- **Beta diversity**: Community composition analysis using Bray-Curtis, Jaccard, and weighted UniFrac distances with multiple normalization methods (relative abundance, CSS, Hellinger)
- **Machine learning**: Random Forest classification for genotype prediction and phenotypic trait association
- **Statistical testing**: PERMANOVA, pairwise comparisons, and dispersion analysis
- **Visualization**: PCoA, NMDS, hierarchical clustering, and correlation analyses
## Data
- **Feature table**: OTU/ASV abundance data (`feature-table.biom`, `OTUs_Table-norm.tab`)
- **Taxonomy**: Taxonomic classification of microbial features (`Taxonomy.txt`)
- **Metadata**: Comprehensive sample metadata including:
  - Genotype information (`English_Name_2`)
  - Geographical origin (`Production_site_2`)
  - Phenotypic traits (fruit color, shape, size, taste, yield, disease resistance)
  - Other traits (seed weight, ovary number, etc.)
- **Phylogenetic tree**: For phylogenetic diversity calculations (`tree.nwk`)
- **DNA sequences**: Raw sequence data (`dna-sequences.txt`)

## Analysis Pipeline

### 1. Data Preprocessing
- Quality filtering and removal of contaminants.
- Taxonomic classification of microbial sequences.
- Normalization and rarefaction analysis to standardize sampling depth.
- Prevalence filtering with a 0.25% threshold.

### 2. Alpha Diversity Analysis
- Rarefaction curves across multiple sequencing depths (500-100,000 reads)
- Diversity metrics: Observed OTUs, Shannon, Simpson, Inverse Simpson, Chao1, ACE
- Statistical comparisons with Dunn's test and significance letters

### 3. Beta-Diversity Analysis
- Multiple distance metrics: Bray-Curtis, Jaccard, UniFrac
- Normalization methods: Relative abundance, CSS, Hellinger transformation
- Ordination: PCoA and NMDS
- Statistical testing: PERMANOVA with pairwise comparisons
- Dispersion analysis: BETADISPER for homogeneity testing

### 4. Machine Learning
- Random Forest models for:
  - Genotype prediction (90% accuracy)
  - Production site classification (82% accuracy)
  - Phenotypic trait prediction (fruit color: 77%, shape: 84%, taste: 84%, etc.)
- Feature importance analysis
- Confusion matrix evaluation

### 5. Microbiome-Phenotype Association
- Genotype-specific microbiome signatures
- Regional variation analysis across Chinese provinces
- Correlation with phenotypic traits
- Differential abundance analysis

## Key Findings
- **Genotype-Specific Microbiomes**: Tomato genotypes harbor distinct and predictable seed microbiome signatures.
- **Phenotype Correlation**: Seed microbiome composition correlates strongly with key phenotypic traits such as fruit color, shape, and taste.
- **Regional Influence**: Geographical location and environmental factors significantly influence microbiome assembly, though genotype effects remain robust.
- **Predictive Power**: Machine learning models can predict a tomato plant's genotype and jeho-phenotypic traits with high accuracy based solely on its seed microbiome profile.

## Usage

### Running the Analysis
1. Clone the repository:
   ```bash
   git clone https://github.com/kaboyo/Tomato-seed-microbiome.git
   cd Tomato-seed-microbiome
   ```

2. Open the R project file in RStudio:
   ```r
   # Open Tomato-seed-microbiome.Rproj
   ```

3. Install required packages (see Dependencies section)

4. Run the main analysis:
```r
source("Tomato100seeds.R")
```

### Reports
- **Main analysis**: Render `Tomato_seed_microbiome-reviews.qmd` for comprehensive results
- **Random Forest analysis**: Available in `Random_for/Random_forest_final.Rmd`
- **HTML output**: View `Tomato_seed_microbiome-reviews.html` for interactive results

### Reproducing Results
All analyses are reproducible with the provided R scripts and data files. The Quarto documents contain embedded code for complete reproducibility.

## Dependencies

This project requires R (version 4.0 or higher). All required packages can be installed by running the following script in your R console:

```r
### Install Dependencies
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
  remotes::install_github(github_packages)
  
  print("All dependencies installed.")
}

install_dependencies()
```

## Project Structure
```
Tomato-seed-microbiome/
├── README.md                           # This file
├── Tomato-seed-microbiome.Rproj        # R project file
├── Tomato_seed_microbiome-reviews.qmd  # Main analysis (Quarto)
├── Tomato_seed_microbiome-reviews.html # Rendered HTML report
├── Tomato100seeds.R                    # Main analysis script
├── Random_forest.Rmd                   # Random forest analysis
├── Random_for/                         # Random forest results
│   ├── ps.RData                        # Processed phyloseq object
│   ├── Random_forest_final.R           # RF implementation
│   └── [various output files]          # Results and figures
├── Alpha-Diversity/                    # Alpha diversity results
├── Beta-Diversity/                     # Beta diversity results
├── outputs_beta/                       # Beta diversity outputs
├── PCoA_individual_plots/              # Ordination plots
├── Original figures/                   # Original figure files
├── Dispersion/                         # Dispersion analysis results
├── Effective_richness/                 # Effective richness data
├── Effective_Shannon/                  # Effective Shannon data
├── Eveness/                           # Evenness analysis
├── PD/                                # Phylogenetic diversity
├── Shannon/                           # Shannon diversity
├── Simpson_index/                     # Simpson diversity
├── richness/                          # Richness data
├── Figures-tomato-seeds.pptx          # PowerPoint figures
├── Metadata.csv                       # Sample metadata
├── feature-table.biom                 # BIOM format feature table
├── OTUs_Table-norm-tax.tab           # Normalized OTU table with taxonomy
├── Taxonomy.txt                       # Taxonomy table
├── tree.nwk                          # Phylogenetic tree
└── [various .tiff and .csv files]     # Figures and data exports
```

## Contributing
Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Citation
If you use this code, data, or analysis in your research, please cite:

[Add appropriate citation when available]

## Authors
- **Expedito Olimi** - Primary researcher and analyst
- Repository maintained by: [GitHub username]

## Contact
For questions, collaborations, or data access requests, please contact the repository maintainer or Expedito Olimi.

---

*This analysis was conducted using R and various bioinformatics tools to explore the complex relationships between tomato genotypes and their seed microbiome communities.*