**Introduction**

#Seed microbiome is crucial for plant health and productivity. It is currently unclear to what extent genotype variation can occur. Here a large-scale approach with 100 different tomato (Solanum lycopersicum L.) genotypes was implemented to explore the plasticity of the seed microbiome.


**Setup and Installation**
**Clone the repository**

```
git clone https://github.com/kaboyo/Tomato-seed-microbiome.git
```
**Navigate to the project directory:**

```
cd Tomato-seed-microbiome
```

**Install the dependencies**

```
Rscript install_dependencies.R
install.packages(c("phyloseq", "vegan", "genefilter", "psych", "DESeq2", "ggplot2", "reshape2", "microbiome", "ranacapa", 
                   "S4Vectors", "Biobase", "BiocGenerics", "parallel", "ggpubr", "decontam", "RColorBrewer", "multcompView", 
                   "metagenomeSeq", "lawstat", "DirichletMultinomial", "lattice", "xtable", "permute", "grid", "doParallel", 
                   "ade4", "car", "gridExtra", "tidyverse", "phia", "lme4", "broom", "AICcmodavg", "multcomp", "rcompanion", 
                   "hrbrthemes", "viridis", "ape", "microbiomeutilities", "microViz", "GUniFrac", "ggtext", "ggraph", "DT", 
                   "corncob", "VennDiagram", "UpSetR", "rstatix", "ALDEx2", "cowplot", "HMP", "dendextend", "rmarkdown", 
                   "dada2", "Biostrings", "ggalluvial", "ggh4x", "ggrepel", "qiime2R", "pals", "ggforce", "vegan", 
                   "data.table", "picante", "ggtreeExtra", "ggpmisc", "DECIPHER", "phangorn", "ggtree", "plotly", "speedyseq", 
                   "ComplexHeatmap", "agricolae", "pairwiseAdonis", "phyloseq.extended"))
```

**Usage**
```
Rscript Tomato100seeds_Final.Rmd
```

**Contributing**    

**Citation**
If you use this project in your work, please cite the original publication.

For more details, view the repository contents.