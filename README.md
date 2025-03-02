**Introduction**

#Seed microbiome is crucial for plant health and productivity. It is currently unclear to what extent plant genotic variations can influence seed microbiome. In a large-scale study (which can be considered the largest seed microbiome dataset), that involve 100 tomato (Solanum lycopersicum L.) genotypes, we examine the plasticity of the tomato seed microbiome.


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