setwd("/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project")
library(phyloseq);packageVersion("phyloseq")#‘1.46.0’
library(vegan);packageVersion("vegan")# ‘2.6.4’
library(genefilter);packageVersion("genefilter")#‘1.84.0’
library(psych);packageVersion("psych")#‘2.4.1’
library(DESeq2);packageVersion("DESeq2")#‘1.42.1’
library(ggplot2);packageVersion("ggplot2")#‘3.5.0’
library(reshape2);packageVersion("reshape2")#‘1.4.4’
library(microbiome);packageVersion("microbiome")#‘1.24.0’
library(ranacapa);packageVersion("ranacapa")#‘0.1.0’
library(S4Vectors);packageVersion("S4Vectors")#‘0.40.2’
library(Biobase);packageVersion("Biobase")#‘2.62.0’
library(BiocGenerics);packageVersion("BiocGenerics")#‘0.48.1’
library(parallel);packageVersion("parallel")#‘4.3.2’
library(ggpubr);packageVersion("ggpubr")#‘0.6.0’
#library(metagMisc)
library(decontam);packageVersion("decontam")#‘1.22.0’
library(RColorBrewer);packageVersion("RColorBrewer")#‘1.1.3’
library(ggpubr);packageVersion("ggpubr")#‘0.6.0’
library(multcompView);packageVersion("multcompView")#‘0.1.10’
library(metagenomeSeq);packageVersion("metagenomeSeq")#‘1.43.0’
#library(devtools)
library(lawstat);packageVersion("lawstat")#‘3.6’
library(DirichletMultinomial);packageVersion("DirichletMultinomial")#‘1.44.0’
library(lattice);packageVersion("lattice")# ‘0.22.5’
library(xtable);packageVersion("xtable")#‘1.8.4’
library(permute);packageVersion("permute")#‘0.9.7’
library(grid);packageVersion("grid")#‘4.3.2’
library("doParallel");packageVersion("doParallel")#‘1.0.17’
#library(igraph)
library(ade4);packageVersion("ade4")#‘1.7.22’
#library(rgl)
library(car);packageVersion("car")#‘3.1.2’
library(gridExtra);packageVersion("gridExtra")# ‘2.3’
library(tidyverse);packageVersion("tidyverse")#‘2.0.0’
library(phia);packageVersion("phia")#‘0.3.1’
library(lme4);packageVersion("lme4")#‘1.1.35.1’
library(broom);packageVersion("broom")#‘1.0.5’
library(AICcmodavg);packageVersion("AICcmodavg")#‘2.3.3’
library(multcomp);packageVersion("multcomp")#‘1.4.25’
library(rcompanion);packageVersion("rcompanion")#‘2.4.35’
library(hrbrthemes);packageVersion("hrbrthemes")#‘0.8.7’
library(viridis);packageVersion("viridis")#‘0.6.5’
library(multcompView);packageVersion("multcompView")#‘0.1.10’#stats
library(ape);packageVersion("ape")#‘5.7.1’#loading the tree
library(microbiomeutilities);packageVersion("microbiomeutilities")#‘1.0.17’#process taxonomy table
library(microViz);packageVersion("microViz")#‘0.12.1’#clean taxonomy table to be cited
library(GUniFrac);packageVersion("GUniFrac")#‘1.8’
library(ggtext) ;packageVersion("ggtext")# ‘0.1.2’for rotated labels on ord_plot() 
library(ggraph);packageVersion("ggraph")# ‘2.2.1’for taxatree_plots()
library(DT);packageVersion("DT")# ‘0.32’for tax_fix_interactive()
library(corncob);packageVersion("corncob")#‘0.4.1’ for beta binomial models in tax_model()
library(VennDiagram);packageVersion("VennDiagram")#‘1.7.3’
library(UpSetR);packageVersion("UpSetR")#‘1.4.0’
library(rstatix);packageVersion("rstatix")#‘0.7.2’
library(ALDEx2);packageVersion("ALDEx2")#[1] ‘1.34.0’
library(cowplot);packageVersion("cowplot")#‘1.1.3’
library(HMP);packageVersion("HMP")#HMP: Hypothesis Testing and Power Calculations for Comparing Metagenomic Samples from HMP #‘2.0.1’
library(dendextend);packageVersion("dendextend")#‘1.17.1’
#citation("rstatix ")
library(dplyr)
library(phyloseq)
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
library(patchwork)
#library(MicrobiotaProcess)
#install.packages("rmarkdown")
library(rmarkdown)
##########################
library(dada2)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(ggalluvial)
library(ggh4x)
library(tidyverse)
library(ggrepel)
#https://github.com/jbisanz/qiime2R
library(qiime2R)
library(dplyr)
library(DESeq2)
library(pals)
library(ggforce)
library(gridExtra)
library(vegan)
library(data.table)
library(ggpubr)
library(picante)
library(ggtreeExtra)
library(ggpmisc)
library(DECIPHER) #BiocManager::install("DECIPHER")
library(phangorn) #remotes::install_github("KlausVigo/phangorn")
library(ggtree) #BiocManager::install("ggtree")
library(plotly)
library(speedyseq) #remotes::install_github("mikemc/speedyseq")
library(ComplexHeatmap)
library(agricolae)
library(scales)
library(ggpmisc)
library(pairwiseAdonis)#library(devtools) install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#remotes::install_github("mahendra-mariadassou/phyloseq-extended", ref = "dev")
library(phyloseq.extended)#remotes::install_github(repo = "mahendra-mariadassou/phyloseq-extended@v0.1.0.9000-beta", ref = "dev")
# install.packages("remotes")
#The following need to be installed but NOT loaded (same naming of different functions):
#library(metagMisc) #devtools::install_github("vmikk/metagMisc")
#library(metagenomeSeq) #BiocManager::install("metagenomeSeq")
#library(ranacapa) #remotes::install_github("gauravsk/ranacapa")
###################################################################################################
###################################################################################################
#Useful website
# https://github.com/joey711/phyloseq
# https://microbiome.github.io/tutorials/
# https://www.microbiomeanalyst.ca/
#using interactive command-line with parameters: sinteractive --partition=highmem --cpus-per-task=30 --mem=400000 --nodes=6 --time=01-02:03:04. 


###################################################################################################
###################################################################################################
source("Rscript_parwise.adonis.txt")
source("aggregate_top_taxa.txt")
source("PD.txt")
source("aggregate_rare.txt")
source("plot_composition.txt")
source("RSCRIPT_parwise.adonis.txt")
source("plot_hierarchy.R")
source("upset_pq.txt")
#installed.packages("rgl")
mycolors=c("grey35","#00AFBB","goldenrod","#0dba2f","#392682","coral3","#8E9CA3","chartreuse4","grey85","#FF2700","bisque4","#0B775E","darkred","#F39B7FB2"
           ,"#FA41CA","steelblue","blue","#D16103", "#C3D7A4", "#52854C","deepskyblue","#293352","#FFDB6D","Yellow","#CC79A7","#392682", "#CAB2D6","maroon"
           ,"orchid1","khaki2","#86486F","#7E6148B2");mycolors

mycola=c("grey35","chartreuse4","#FE7F9D","darkturquoise", "#F0E442", "#009E73", "#FF0000",
         "Orange","blue","#850dba","green1","#FF6AD5", "#6A3D9A","Yellow","Blue",
         "deeppink","#1C0221","#00AFBB","steelblue","#D16103", "#C3D7A4", "#52854C",
         "#4E84C4","#293352","#FFDB6D","#CC79A7", "#C4961A", "#F4EDCA","dodgerblue2", "#E31A1C", # red
         "green4","#0dba2f", # purple"#FF7F00", # orange
         "black", "gold1", "skyblue2", "#FB9A99", # lt pink "palegreen2",
         "#CAB2D6", # lt purple "#FDBF6F", # lt orange "gray70", "khaki2","maroon", "orchid1", "deeppink1", 
         "blue1", "steelblue4", "Sky Blue",  "grey85","#3B9AB2", "goldenrod", "yellow3",
         "darkorange4", "brown","#999999", "#E69F00", "#56B4E9", "#009E73","#F0E442", 
         "#D55E00", "#E69F00", "#ba0d42",
         "yellow4", "#D55E00", "#CC79A7","darkred" ,
         "#56B4E9","#86486F", "#446455", "bisque4","#F7B1AB",
         "#807182","#E3DDDD", "#A45C3D","coral3","#000099","#E2C59F",
         "#B6C5CC", "#8E9CA3", "#556670","#771C19", "#AA3929", "#E25033", 
         "#F27314", "#F8A31B","#2D2D2D","#91D1C2B2","#DC0000B2", "#7E6148B2",
         "#E64B35B2", "#4DBBD5B2", "#00A087B2", "#3C5488B2","#F39B7FB2","#8491B4B2",
         "#FF2700","#62BBA5","#B58900","#80CDC1","#8073AC","#0072B2","#392682",
         "#FA41CA","#0B775E","#02401B","#F4B5BD","#B6854D","deepskyblue");mycola

mycol = c("bisque4", "deeppink", "chartreuse4","blue","coral3","goldenrod","deepskyblue","darkred","#0dba2f", "Orange","Yellow", "#850dba","Blue", "Sky Blue", "#ba0d42");mycol
###################################################################################################
###################################################################################################
###################################################################################################
## Imported otu table
Bacteria_otu <- read.table("Table.txt", header=TRUE,row.names="OTUID");Bacteria_otu
Bacteria_otu <- otu_table(Bacteria_otu, taxa_are_rows = TRUE);Bacteria_otu 

## Imported taxonomy
#Bacteria_tax <- read.csv("Taxonomy.csv",header=TRUE,row.names="OTUID");Bacteria_tax
Bacteria_tax <- read.table("Taxonomy.txt", header=TRUE, sep="\t",row.names="OTUID");Bacteria_tax
Bacteria_tax <- as.matrix(Bacteria_tax);Bacteria_tax
Bacteria_tax <- tax_table(Bacteria_tax);Bacteria_tax

## Imported metadata file
#metadata <- read.table("Metadata.txt", header=TRUE, sep="\t",row.names="SampleID");metadata
metadata <- read.csv("Metadata.csv",header=TRUE,row.names="SampleID");metadata
#rownames(metadata) <- metadata[,1]
metadata <- sample_data(metadata);metadata
#tree
tree <- read.tree(file = "tree.nwk");tree
tree <- phy_tree(tree)

all(rownames(metadata) %in% colnames(Bacteria_otu))
#Convert to phyloseq object
Bacteria <- merge_phyloseq(Bacteria_otu,Bacteria_tax, metadata,tree);Bacteria#[ 16460 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Kingdom=="Bacteria");Bacteria# [ 6448 taxa and 84 samples ]
Bacteria = subset_taxa(Bacteria, Phylum!="Cyanobacteria");Bacteria#  [ 15588 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Order!="Chloroplast");Bacteria#  [ 15588 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Family!="Mitochondria");Bacteria# [ 15588 taxa and 1197 samples ]
Bacteria <- prune_taxa(taxa_sums(Bacteria)>0, Bacteria);Bacteria# [ 15588 taxa and 1197 samples ]
#OTU_table_Bacteria = as(otu_table(Bacteria), "matrix")
#write the phyloseq object
saveRDS(Bacteria, "Bacteria_TomatoSeed_Microbiome.rds")
#this can be read into R as
#ps = readRDS("Bacteria_TomatoSeed_Microbiome.rds")
###################################################################################################
###################################################################################################

#Summarize 
microbiome::summarize_phyloseq(Bacteria)
#Number of singletons = 68"
#Sparsity = 0.995179763625127"
#Average number of reads = 22794.9197994987"
# Total number of reads = 27281011"
#Max. number of reads = 225892"
#Min. number of reads = 68"
taxonomy <- as.vector(phyloseq::tax_table(Bacteria))

sums_Bacteria <- sample_sums(Bacteria);sums_Bacteria
sums_Bacteria <- data.frame(sums_Bacteria);sums_Bacteria
sums_Bacteria <- as.matrix(sums_Bacteria);sums_Bacteria
colSums(sums_Bacteria) #  27281011      
max(sums_Bacteria) #225892
min(sums_Bacteria) # 68
sums_Bacteria

#If you need to convert back the phyloseq data to individual files
Bacteria_table = as(otu_table(Bacteria), "matrix")
# Coerce to data.frame
Bacteria_table = as.data.frame(Bacteria_table)
write.csv(Bacteria_table ,file = "Bacteria_table.csv")

#Bacteria_tax = as(tax_table(Bacteria), "matrix")
# Coerce to data.frame
Bacteria_tax = as.data.frame(Bacteria_tax)
write.csv(Bacteria_tax ,file = "Bacteria_tax.csv")

Bacteria_metadata = as(sample_data(Bacteria), "matrix")
# Coerce to data.frame
Bacteria_metadata = as.data.frame(Bacteria_metadata)
write.csv(Bacteria_metadata ,file = "Bacteria_metadata.csv")
###################################################################################################
###################################################################################################
#Rarefaction
set.seed(1024)
rareres <- get_rarecurve(obj=Bacteria, chunks=400)
prare1 <- ggrarecurve(obj=rareres,
                      factorNames="Production_site_2",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1")
) +
  scale_color_manual(values=mycola)+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey33"),
        strip.text.x = element_text(face="bold"))+
  theme(text=element_text(size=14,  family="sans"))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=6))+
  labs(tag = "A")+
  theme(legend.position = "bottom")+ 
  theme(legend.position="none");prare1
###################################################################################################
###################################################################################################
prare2 <- ggrarecurve(obj=rareres,
                      factorNames="English_Name_2",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1")
) +
  scale_color_manual(values=mycola)+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey33"),
        strip.text.x = element_text(face="bold"))+
  theme(text=element_text(size=14,  family="sans"))+ 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=2))+
  labs(tag = "B")+
  theme(legend.position = "bottom")+
  theme(legend.position="none");prare2

tiff("Rarefaction.tiff", units="in", width=9, height=6, res=300)
Rarefaction.tiff=prare2+plot_layout(nrow =1);Rarefaction.tiff
ggsave("Rarefaction.tiff", plot = Rarefaction.tiff)
dev.off()
###################################################################################################
###################################################################################################
###Coool rarefaction
#https://microsud.github.io/microbiomeutilities/articles/microbiomeutilities.html
library(microbiomeutilities)
p0 <- Bacteria
# set seed
set.seed(1)
subsamples <- seq(0, 50000, by=500)[-1]
#subsamples = c(10, 5000, 10000, 20000, 30000)

p <- plot_alpha_rcurve(p0, index="observed",
                       subsamples = subsamples,
                       lower.conf = 0.025, 
                       upper.conf = 0.975,
                       group="Production_site_2",
                       label.color = "brown3",
                       label.size = 3,
                       label.min = TRUE) 

p <- p + scale_color_manual(values = mycola) + 
  scale_fill_manual(values = mycola)+ 
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey33"),
        strip.text.x = element_text(face="bold"))+
  theme(text=element_text(size=14,  family="sans"))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=6))+
  theme(legend.position="none")+
  labs(tag = "A")
print(p)


tiff("p.tiff", units="in", width=10, height=6, res=300)
ggsave("p.TIFF", plot = p)

###################################################################################################
###################################################################################################
p2<- plot_alpha_rcurve(p0, index="observed",
                     subsamples = subsamples,
                     lower.conf = 0.025, 
                     upper.conf = 0.975,
                     group="English_Name_2",
                     label.color = "brown3",
                     label.size = 3,
                     label.min = TRUE) 

p2 <- p2 + scale_color_manual(values = mycola) + 
  scale_fill_manual(values = mycola)+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey33"),
        strip.text.x = element_text(face="bold"))+
  theme(text=element_text(size=16,  family="sans"))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=6))+
  theme(legend.position="none")+labs(tag = "A")
print(p2)

tiff("p2.tiff", units="in", width=10, height=6, res=300)
ggsave("p2.TIFF", plot = p2)
###################################################################################################
###################################################################################################
p3<- plot_alpha_rcurve(p0, index="observed",
                       subsamples = subsamples,
                       lower.conf = 0.025, 
                       upper.conf = 0.975,
                       group="TMV_resistance",
                       label.color = "brown3",
                       label.size = 3,
                       label.min = TRUE) 

p3 <- p3 + scale_color_manual(values = mycola) + 
  scale_fill_manual(values = mycola)+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey33"),
        strip.text.x = element_text(face="bold"))+
  theme(text=element_text(size=16,  family="sans"))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=6))+
  theme(legend.position="none")+labs(tag = "B")
print(p3)

tiff("TMV_RES.tiff", units="in", width=12, height=6, res=300)
ggsave("TMV_RES.tiff", plot = p3)
###################################################################################################
###################################################################################################
#Plot correlation of shannon index and total reads
#Shannon diversity, and phylogenetic diversity on the subsampled data (since this is common practice).
corr1=ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(Bacteria),
                               "observed" = phyloseq::estimate_richness(Bacteria, measures = "Observed")[, 1]),
             aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Observed Richness\n")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=15),
        axis.title.y = element_text(face="bold", colour="gray33", size=15),
        axis.text.x = element_text(colour="gray33", size=15),
        axis.text.y  = element_text(colour="gray33", size=15))+
  labs(tag = "A");corr1

#############
corr2=ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(Bacteria),
                               "Shannon" = phyloseq::estimate_richness(Bacteria, measures = "Shannon")[, 1]),
             aes(x = total_reads, y = Shannon)) +
  geom_point() +
  geom_smooth(method="lm", se = FALSE) +
  labs(x = "\nTotal Reads", y = "Shannon index\n")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=15),
        axis.title.y = element_text(face="bold", colour="gray33", size=15),
        axis.text.x = element_text(colour="gray33", size=15),
        axis.text.y  = element_text(colour="gray33", size=15))+
  labs(tag = "B");corr2

tiff("correlation.tiff", units="in", width=12, height=6, res=300)
correlation.tiff=corr1+corr2+plot_layout(ncol=2);correlation.tiff
ggsave("correlation.tiff", plot = correlation.tiff)
dev.off()

###################################################################################################
###################################################################################################
#Obtained from: https://f1000research.com/articles/5-1492/v1;#One of the reasons to filter by prevalence is to avoid spending much time analyzing taxa that were only rarely seen. This also turns out to be a useful filter of noise (taxa that are actually just artifacts of the data collection process)
# Define prevalence of each taxa;# (in how many samples did each taxa appear at least once)
# Define prevalence of each taxa; # (in how many samples did each taxa appear at least once)
ps=Bacteria
prev0 = apply(X = otu_table(ps),
              MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})

prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))# the most prevalent ASVs present therein the composition

#only for reference, in case we want to filter Phyla that appear less than X times
keepPhyla = table(prevdf$Phylum)[(table(prevdf$Phylum) > 0)]

prevdf1 = subset(prevdf, Phylum %in% names(keepPhyla))

# Define prevalence threshold as 0.25% of total samples
prevalenceThreshold = 0.0025 * nsamples(ps)
prevalenceThreshold

# Execute prevalence filter, using `prune_taxa()` function
ps1 = prune_taxa((prev0 > prevalenceThreshold), ps)
ps1

# Filter entries with unidentified Phylum.
ps2 = subset_taxa(ps1, Phylum %in% names(keepPhyla))
ps2

ggplot(prevdf1, aes(TotalAbundance, Prevalence, color = Phylum)) +
  geom_hline(yintercept = prevalenceThreshold, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.65) +
  scale_y_log10() + scale_x_log10() +
  ggtitle(paste0("Prevalence filtering, using a threshold of ", prevalenceThreshold))+
  xlab("Total Abundance") +
  scale_color_manual(values=ocean.balance(38))+
  theme_bw()+
  facet_wrap(~Phylum, nrow = 4)
########################################################################################################
print(paste0("Initial #ASVs:", ntaxa(ps), ". #ASVs after prevalence filtering:", ntaxa(ps2), ". #ASVs removed: ", ntaxa(ps)-ntaxa(ps2)))
###################################################################################################
###################################################################################################
#Save the ps2 otu_table and relative abundance otu table 
#write.table(otu_table(ps2), "otu_table_ps2.txt", quote = F, sep = "\t")

#Or transformed to relative abundances
ps2.prop = transform_sample_counts(ps2, function(x) {x/sum(x)*100});ps2.prop
#write.table(otu_table(ps2.prop), "otu_table_ps2_prop.txt", quote = F, sep = "\t")
###################################################################################################
# 6. Bar plots:
##Showing individual replicates. Top 500 taxa
toptax = names(sort(taxa_sums(ps2), decreasing=TRUE))[1:500];toptax
ps.toptax = transform_sample_counts(ps2, function(x) {x/sum(x)*100});ps.toptax
ps.toptax <- prune_taxa(toptax, ps.toptax);ps.toptax

#Class level
psmelt(ps.toptax)
gg1=ggplot(psmelt(ps.toptax), aes(x=reorder(Sample, English_Name), y=Abundance, fill=Class))+
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values=c(ocean.dense(7), ocean.tempo(8)))+
  #scale_fill_manual(values=c("#25725d", "#439b5b", "#7bc369", "#abd69a", "#d8e8b0", "#cccda7", "#bfb77b", "#7e9998", "#92B2B2", "#83bbbf", "#61acc1", "#58a0c4", "#7087c3","#8973b4", "#512972"))+
  #scale_fill_manual(values = mycola)+
  xlab("") + ylab("Relative abundance (%)") +
  facet_wrap(~reorder(Production_site_2, Sample), scales = "free_x", nrow = 1) + 
  ggtitle("Class, top 500 ASVs") +
  theme_pubclean() + theme(legend.position="bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0)) +
theme(legend.text = element_text(face = c(rep("italic", 5), rep("plain", 5))))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow = 2))+
  theme(text=element_text(size=8,  family="Arial"))+ labs(tag = "A");gg1

library(ggpubr)
tiff("gg1_Region_cultivar.tiff", units="in", width=10, height=6, res=300)
gg1_Region_cultivar.tiff=gg1+plot_layout(nrow = 1);gg1_Region_cultivar.tiff
ggsave("gg1_Region_cultivar.tiff", plot = gg1_Region_cultivar.tiff)
dev.off()
###################################################################################################
###################################################################################################

## 6.2 Export Taxa and ASVs table
#Export the tax_table and otu_table at the genus level
ps2.prop.tax = tax_table(ps2.prop) %>% as.data.frame() %>% rownames_to_column("ASVs");ps2.prop.tax
ps2.prop.asv = t(otu_table(ps2.prop)) %>% as.data.frame();ps2.prop.asv
colnames(ps2.prop.asv) = sample_data(ps2.prop)$English_Name_2
ps2.prop.asv =  rownames_to_column(ps2.prop.asv, "ASVs")

#write.table(left_join(ps2.prop.tax, ps2.prop.asv, by="ASVs"), "Tax_ASV_table.txt", quote = F, sep = "\t", row.names = F)

###################################################################################################
## 6.3 CSS normalization
#For certain downstream analyses; for example to compare different abundances of taxa between samples, the counts needs to be normalized to account for different sampling efforts, exposures, baselines, etc. We will normalize the samples using the CSS (Cumulative Sum Scaling) introduced in (Paulson et al. 2013), where the offset of a sample is the cumulative sum of counts in that sample, up to a quantile determined in a data driven way. Calculates scaling factors as the cumulative sum of ASVs abundances up to a data-derived threshold.
#Ussing the CSS
ps2.css = metagMisc::phyloseq_transform_css(ps2, norm = T, log = F);ps2.css#log2 FALSE!
ps2.css = orient_taxa(ps2.css, "columns") #change taxa to columns
#and the the TSS (proportion; relative abundance)
ps2.css.prop = transform_sample_counts(ps2.css, function(x) {x/sum(x)*100});ps2.css.prop


gg2=ggplot(psmelt(prune_taxa(toptax, ps2.css.prop)), aes(x=reorder(English_Name_2, Sample), y=Abundance, fill=Class))+
  geom_bar(stat = "identity", position = "stack", width = 1) +
  scale_fill_manual(values = mycola)+
  #scale_fill_manual(values=c(ocean.dense(7), ocean.tempo(8)))+
  #scale_fill_manual(values=c("#25725d", "#439b5b", "#7bc369", "#abd69a", "#d8e8b0", "#cccda7", "#bfb77b", "#7e9998", "#92B2B2", "#83bbbf", "#61acc1", "#58a0c4", "#7087c3", "#8973b4", "#512972"))+
  xlab("") + ylab("CSS-Normalized relative abundance (%)") +
  facet_wrap(~reorder(Production_site_2, Sample), scales = "free_x", nrow = 1) + 
  ggtitle("Class, top 500 ASVs") + xlab("CSS-Normalized relative abundance (%)")+
theme_pubclean() + theme(legend.position="top", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
theme(legend.text = element_text(face = c(rep("italic", 5), rep("plain", 5))))+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow = 2))+
  theme(text=element_text(size=10,  family="Arial"))+ labs(tag = "A");gg2

# Move to a new page
grid.newpage()
# Create layout : nrow = 2, ncol = 2
pushViewport(viewport(layout = grid.layout(2, 2)))
# A helper function to define a region on the layout
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
# Arrange the plots
print(gg1, vp=define_region(1, 1:2))
print(gg2, vp = define_region(2, 1))
#print(ydensity, vp = define_region(2, 2))
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


#https://microbiome.github.io/tutorials/CoremicrobiotaAmplicon.html
#consider the abundance and occupancy approach
library(microbiome)
#Bacteria
# keep only taxa with positive sums
Bac <- prune_taxa(taxa_sums(ps2) > 0, ps2)
print(Bac)

# Calculate compositional version of the data
# (relative abundances)
Bac_rel <- microbiome::transform(Bac, "compositional");Bac_rel

############################################
#Linegraph visualization of the Core microbiome
# With compositional (relative) abundances
det <- c(0, 0.01, 0.02, 0.1, 0.5, 2, 5, 20)/100;det
prevalences <- seq(.05, 1, .05);prevalences

plot_core(Bac_rel, prevalences = prevalences, 
          detections = det, plot.type = "lineplot") + 
  xlab("Relative Abundance (%)") + 
  theme_bw()
############################################
# 0.0002, prevalence = 30/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 30/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .3);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_30.csv")
######################################################################################################################
######################################################################################################################
# 0.0002, prevalence = 40/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 40/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .4);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_40.csv")
########################################################################################################
######################################################################################################################
# 0.0002, prevalence = 50/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 50/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .5);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_50.csv")
######################################################################################################################
########################################################################################################
# 0.0002, prevalence = 60/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 60/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .6);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_60.csv")
######################################################################################################################
######################################################################################################################
# 0.0002, prevalence = 80/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 80/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .8);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_80.csv")
########################################################################################################
######################################################################################################################
# 0.0002, prevalence = 90/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 90/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .9);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_90.csv")
###################################################################################################
###################################################################################################
# 0.0002, prevalence = 95/100
core.taxa.standard <- core_members(Bac_rel, detection = 0.0002, prevalence = 95/100);core.taxa.standard
######################################################################################################################
#A full phyloseq object of the core microbiota is obtained as follows:
pseq.core <- core(Bac_rel, detection = 0.0002, prevalence = .95);pseq.core

#Retrieving the associated taxa names from the phyloseq object:
core.taxa <- taxa(pseq.core);core.taxa
class(core.taxa)

# get the taxonomy data
tax.mat <- tax_table(pseq.core);tax.mat
tax.df <- as.data.frame(tax.mat);tax.df

# add the OTus to last column
tax.df$OTU <- rownames(tax.df)

# select taxonomy of only 
# those OTUs that are core memebers based on the thresholds that were used.
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
write.csv(core.taxa.class,"core.taxa.class_95.csv")
######################################################################################################################
######################################################################################################################
library(UpSetR)
library(VennDiagram)
vennlist <- get_vennlist(obj=Bacteria, factorNames="English_Name_2")
upsetda <- get_upset(obj=ps, factorNames="Group")
library(VennDiagram)
library(UpSetR)
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "#FD9347"),
                      cat.col=c("#00AED7", "#FD9347"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
grid::grid.draw(vennp)

###################################################################################################
###################################################################################################
## Imported otu table
Bacteria_otu <- read.table("Table.txt", header=TRUE,row.names="OTUID");Bacteria_otu
## Imported taxonomy
#Bacteria_tax <- read.csv("Taxonomy.csv",header=TRUE,row.names="OTUID");Bacteria_tax
Bacteria_tax <- read.table("Taxonomy.txt", header=TRUE, sep="\t",row.names="OTUID");Bacteria_tax
## Imported metadata file
metadata <- read.table("Metadata.txt", header=TRUE, sep="\t",row.names="SampleID");metadata
#metadata <- read.csv("Metadata.csv",header=TRUE,row.names="SampleID");metadata

###################################################################################################
###################################################################################################
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/ 
#Agglomerate to phylum-level and rename
Bacteria_phylum <- phyloseq::tax_glom(ps2, "Phylum")#[ 38 taxa and 1197 samples ]
phyloseq::taxa_names(Bacteria_phylum) <- phyloseq::tax_table(Bacteria_phylum)[, "Phylum"]
phyloseq::otu_table(Bacteria_phylum)[1:5, 1:5]
sample_data(Bacteria_phylum)

#Melt and plot
phyloseq::psmelt(Bacteria_phylum) %>%
  ggplot(data = ., aes(x = Production_site_2, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")+
  scale_color_manual(values=mycola)
########################################################################################################################
#https://bioinfo.ird.fr/index.php/trainings-fr/trainings-2019-metabarcoding/trainings-2019-metabarcoding-practice/
###################################################################################################
###################################################################################################
#Beta-diversity
###################################################################################################
#Hellinger transformation:method="hellinger")
#hellinger transform (centre log ration transformation)
(ps_hellinger <- microbiome::transform(ps2, "hellinger"))  
phyloseq::otu_table(ps2)[1:5, 1:5]
##PCA via phyloseq
ord_hellinger <- phyloseq::ordinate(ps_hellinger, "RDA")
#Plot scree plot
phyloseq::plot_scree(ord_hellinger) + 
  geom_bar(stat="identity", fill = "blue") +
  labs(x = "\nAxis", y = "Proportion of Variance\n")

head(ord_hellinger$CA$eig)                                                  
#PC1       PC2       PC3       PC4       PC5       PC6 
#216.73977  73.60736  59.27932  46.48439  30.78463  30.62497 
sapply(ord_hellinger$CA$eig[1:5], function(x) x / sum(ord_hellinger$CA$eig))     
#PC1        PC2        PC3        PC4        PC5 
#0.09572831 0.03251045 0.02618213 0.02053094 0.01359677 

#Scale axes and plot ordination
hellinger1 <- ord_hellinger$CA$eig[1] / sum(ord_hellinger$CA$eig)
hellinger2 <- ord_hellinger$CA$eig[2] / sum(ord_hellinger$CA$eig)
phyloseq::plot_ordination(Bacteria, ord_hellinger, type="samples", color="Production_site_2_2") + 
  geom_point(size = 2) +
  coord_fixed(hellinger2 / hellinger1) +
  stat_ellipse(aes(group = Production_site_2), linetype = 2)

phyloseq::plot_ordination(Bacteria, ord_hellinger, type="samples", color="English_Name_2_2") + 
  geom_point(size = 2) +
  coord_fixed(hellinger2 / hellinger1) +
  stat_ellipse(aes(group = Production_site_2), linetype = 2)


#Generate distance matrix
hellinger_dist_matrix <- phyloseq::distance(ps_hellinger, method = "euclidean") 
#ADONIS test
vegan::adonis2(hellinger_dist_matrix ~ phyloseq::sample_data(ps_hellinger)$English_Name_2, permutation=999)
#                                                     Df SumOfSqs      R2      F Pr(>F)    
#phyloseq::sample_data(ps_hellinger)$English_Name_2   99   410.92 0.56909 14.634  0.001 ***
#Residual                                           1097   311.14 0.43091                  
#Total                                              1196   722.06 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#Dispersion test and plot
dispr <- vegan::betadisper(hellinger_dist_matrix, phyloseq::sample_data(ps_hellinger)$English_Name_2);dispr #calculate the dispersion
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Aitchison Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
permutest(dispr)
#             Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
#Groups      99  77206  779.86 12.685    999  0.001 ***
#Residuals 1097  67440   61.48                         
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#There is significant difference in the centroid location according to the the seed cultivars
# Load data
#:::::::::::::::::::::::::::::::::::::::
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
data("ToothGrowth")
df <- ToothGrowth

# Kruskal-wallis rank sum test
#:::::::::::::::::::::::::::::::::::::::::
df %>% kruskal_effsize(len ~ dose)

# Grouped data
df %>%
  group_by(supp) %>%
  kruskal_effsize(len ~ dose)
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#################################################################################################
#https://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.htm
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
library(MicrobiotaProcess)
#Will consider the transformation by hellinger approach
distme <- get_dist(ps2, distmethod ="bray", method="hellinger");distme
sampleda <- data.frame(sample_data(ps2), check.names=FALSE);sampleda
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE];sampleda
sampleda$Production_site_2 <- factor(sampleda$Production_site_2)
adores <- adonis2(distme ~ Production_site_2, data=sampleda, permutation=999);adores#change from adonis to adonis2
#                     Df SumOfSqs      R2      F Pr(>F)    
#Production_site_2   11   35.537 0.11415 13.882  0.001 ***
#Residual          1185  275.776 0.88585                  
#Total             1196  311.313 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1                
###################################################################################################
load(file = "alpha-obs-gen-long.RData")
###################################################################################################
#Consider the cultivars
distme <- get_dist(ps2, distmethod ="bray", method="hellinger");distme
sampleda <- data.frame(sample_data(ps2), check.names=FALSE);sampleda
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE];sampleda
sampleda$English_Name_2 <- factor(sampleda$English_Name_2)
adores <- adonis2(distme ~ English_Name_2, data=sampleda, permutation=999);adores#change from adonis to adonis2
#                 Df  SumOfSqs    R2    F    Pr(>F)    
#English_Name_2   99   184.77 0.59353 16.18  0.001 ***
#Residual       1097   126.54 0.40647                 
#Total          1196   311.31 1.00000                 
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1@@@@*** 
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

#***********Convert to relative abundance or the computation of the micribial composition in the samples
GPf = tax_glom(ps, "Class") %>% transform_sample_counts(function(x) {x * 100/sum(x)})

# The taxa orders are the same in tax_table() and otu_table()
# Use 'rowMeans(otu_table(GPf))' to calculate per-row average in the OTU table
df = data.frame(Class = tax_table(GPf)[,"Class"], Mean = rowMeans(otu_table(GPf)), row.names = NULL)
df = df[order(-df$Mean),]
head(df)#***********
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
# "unifrac", "wunifrac", "manhattan", "euclidean", "canberra", "bray", "kulczynski" ...(vegdist, dist)
pcares <- get_pca(obj=Bacteria, method="hellinger");pcares
pcoares <- get_pcoa(obj=Bacteria, distmethod="bray", method="hellinger");pcoares#pcoa

#PCA_pcoa_countries of origign
#Visulizing the result
tiff("PCA_PCOA_Region.tiff", units="in", width=12, height=10, res=300)
pcaplot1 <- ggordpoint(obj=pcares, pc=c(1, 2), biplot=F, speciesannot=F,
                       factorNames=c("Production_site_2"), ellipse=F)+
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "A");pcaplot1
#ggsave("pcaplot1.TIFF", plot = pcaplot1)
# pc = c(1, 3) to show the first and third principal components.

pcaplot2 <- ggordpoint(obj=pcares, pc=c(1, 3), biplot=F, speciesannot=F,
                       factorNames=c("Production_site_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "B");pcaplot2
#ggsave("pcaplot2.TIFF", plot = pcaplot2);pcaplot2


pcoaplot1 <- ggordpoint(obj=pcoares, biplot=F, speciesannot=F,
                        factorNames=c("Production_site_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=14,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "C");pcoaplot1

# pc = c(1, 3) to show the first and third principal components.

pcoaplot2 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=F, speciesannot=F,
                        factorNames=c("Production_site_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "D");pcoaplot2

library(ggpubr)
tiff("PCA_PCOA_Region.tiff", units="in", width=10, height=10, res=300)
PCA_PCOA_Region.tiff=pcaplot1+pcaplot2+pcoaplot1+pcoaplot2+plot_layout(ncol = 2);PCA_PCOA_Region.tiff
ggsave("PCA_PCOA_Region.tiff", plot = PCA_PCOA_Region.tiff)
dev.off()

###################################################################################################
###################################################################################################
##Consider the different cultivars

pcoaplot5 <- ggordpoint(obj=pcares, pc=c(1, 2), biplot=F, speciesannot=F,
                       factorNames=c("English_Name_2"), ellipse=F)+
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "A");pcoaplot5
#ggsave("pcaplot1.TIFF", plot = pcaplot1)
# pc = c(1, 3) to show the first and third principal components.

pcoaplot6 <- ggordpoint(obj=pcares, pc=c(1, 3), biplot=F, speciesannot=F,
                       factorNames=c("English_Name_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "B");pcoaplot6
#ggsave("pcaplot2.TIFF", plot = pcaplot2);pcaplot2


pcoaplot3 <- ggordpoint(obj=pcoares, biplot=F, speciesannot=F,
                        factorNames=c("English_Name_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "C");pcoaplot3
#ggsave("pcoaplot3.TIFF", plot = pcoaplot3);pcoaplot3

# pc = c(1, 3) to show the first and third principal components.
pcoaplot4 <- ggordpoint(obj=pcoares, pc=c(1, 3), biplot=F, speciesannot=F,
                        factorNames=c("English_Name_2"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "D");pcoaplot4
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################

tiff("PCOA_cultitvars.tiff", units="in", width=8, height=9, res=300)
PCOA_cultitvars.tiff=pcoaplot5+pcoaplot6+pcoaplot3+pcoaplot4+plot_layout(ncol = 2);PCOA_cultitvars.tiff
ggsave("PCOA_cultitvars.tiff", plot = PCOA_cultitvars.tiff)
dev.off()

###################################################################################################
###################################################################################################
Bacteria_ovaries <- subset_samples(Bacteria, !is.na(ovaries));Bacteria_ovaries
distme_ovaries <- get_dist(Bacteria_ovaries, distmethod ="bray", method="hellinger");distme_ovaries
sampleda_ovaries <- data.frame(sample_data(Bacteria_ovaries), check.names=FALSE);sampleda_ovaries
sampleda_ovaries <- sampleda_ovaries[match(colnames(as.matrix(distme_ovaries)),rownames(sampleda_ovaries)),,drop=FALSE];sampleda_ovaries
sampleda_ovaries$ovaries <- factor(sampleda_ovaries$ovaries)
adores <- adonis2(distme_ovaries ~ ovaries, data=sampleda_ovaries, permutation=999);adores#change from adonis to adonis2
#           Df SumOfSqs      R2      F Pr(>F)    
#ovaries     2    5.045 0.01694 9.5556  0.001 ***
#Residual 1109  292.753 0.98306                  
#Total    1111  297.798 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_ovaries=c("#00AFBB","#FB9A99","chartreuse4");myco_ovaries
pcoares_ovaries <- get_pcoa(obj=Bacteria_ovaries, distmethod="bray", method="hellinger");pcoares_ovaries#pcoa
pcoaplot_ovaries <- ggordpoint(obj=pcoares_ovaries, biplot=F, speciesannot=F,
                               factorNames=c("ovaries"), ellipse=F) +
  scale_fill_manual(values=myco_ovaries)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "A");pcoaplot_ovaries

###################################################################################################
###################################################################################################
#Fruitshape
Bacteria_Fruit_shape <- subset_samples(Bacteria, !is.na(Fruit_shape));Bacteria_Fruit_shape
distme_Fshape <- get_dist(Bacteria_Fruit_shape, distmethod ="bray", method="hellinger");distme_Fshape
sampleda_Fshape <- data.frame(sample_data(Bacteria_Fruit_shape), check.names=FALSE);sampleda_Fshape
sampleda_Fshape <- sampleda_Fshape[match(colnames(as.matrix(distme_Fshape)),rownames(sampleda_Fshape)),,drop=FALSE];sampleda_Fshape
sampleda_Fshape$Fruit_shape <- factor(sampleda_Fshape$Fruit_shape)
adores <- adonis2(distme_Fshape ~ Fruit_shape, data=sampleda_Fshape, permutation=999);adores#change from adonis to adonis2
#               Df SumOfSqs     R2      F Pr(>F)    
#Fruit_shape    2     4.37 0.0136 8.1494  0.001 ***
#Residual    1182   317.21 0.9864                  
#Total       1184   321.58 1.0000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_shp=c("darkturquoise","chartreuse4","#FB9A99");myco_shp
pcoares_Fshape <- get_pcoa(obj=Bacteria_Fruit_shape, distmethod="bray", method="hellinger");pcoares_Fshape#pcoa
pcoaplot_Fshape <- ggordpoint(obj=pcoares_Fshape, biplot=F, speciesannot=F,
                              factorNames=c("Fruit_shape"), ellipse=F) +
  scale_fill_manual(values=myco_shp)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "B");pcoaplot_Fshape
##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
###############################################################################################################################################################################################################
#Fruit color
tiff("PCOA_Parameters.tiff", units="in", width=20, height=20, res=300)
Bacteria_Fruit_color <- subset_samples(Bacteria, !is.na(Fruit_color));Bacteria_Fruit_color
distme_Fcol <- get_dist(Bacteria_Fruit_color, distme_Fcolthod ="bray", method="hellinger");distme_Fcol
sampleda_Fcol <- data.frame(sample_data(Bacteria_Fruit_color), check.names=FALSE);sampleda_Fcol
sampleda_Fcol <- sampleda_Fcol[match(colnames(as.matrix(distme_Fcol)),rownames(sampleda_Fcol)),,drop=FALSE];sampleda_Fcol
sampleda_Fcol$Fruit_color <- factor(sampleda_Fcol$Fruit_color)
adores <- adonis2(distme_Fcol ~ Fruit_color, data=sampleda_Fcol, permutation=999);adores#change from adonis to adonis2
#             Df SumOfSqs      R2      F Pr(>F)    
#Fruit_color    4    15.05 0.02084 6.2792  0.001 ***
#Residual    1180   707.05 0.97916                  
#Total       1184   722.10 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_col=c("grey35","#00AFBB", "yellow3","chartreuse4","#FB9A99");myco_col
pcoares_Fcol <- get_pcoa(obj=Bacteria_Fruit_color, distmethod="bray", method="hellinger");pcoares_Fcol#pcoa
pcoaplot_Fcol <- ggordpoint(obj=pcoares_Fcol, biplot=F, speciesannot=F,
                        factorNames=c("Fruit_color"), ellipse=F) +
  scale_fill_manual(values=mycola)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "C");pcoaplot_Fcol

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
Bacteria_Fruit_taste <- subset_samples(Bacteria, Fruit_taste %in% c("Sweet", "Sour"));Bacteria_Fruit_taste#[ 15588 taxa and 1102 samples ]
distme_Fruit_taste <- get_dist(Bacteria_Fruit_taste, distmethod ="bray", method="hellinger");distme_Fruit_taste
sampleda_Fruit_taste <- data.frame(sample_data(Bacteria_Fruit_taste), check.names=FALSE);sampleda_Fruit_taste
sampleda_Fruit_taste <- sampleda_Fruit_taste[match(colnames(as.matrix(distme_Fruit_taste)),rownames(sampleda_Fruit_taste)),,drop=FALSE];sampleda_Fruit_taste
sampleda_Fruit_taste$Fruit_taste <- factor(sampleda_Fruit_taste$Fruit_taste)
adores <- adonis2(distme_Fruit_taste ~ Fruit_taste, data=sampleda_Fruit_taste, permutation=999);adores#change from adonis to adonis2
#                       Df SumOfSqs     R2      F Pr(>F)    
#Fruit_taste    2    5.101 0.0174 9.6611  0.001 ***
#Residual            1091  288.020 0.9826                  
#Total               1093  293.121 1.0000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_taste=c("#FB9A99","chartreuse4");myco_taste
pcoares_Fruit_taste <- get_pcoa(obj=Bacteria_Fruit_taste, distmethod="bray", method="hellinger");pcoares_Fruit_taste#pcoa
pcoaplot_Fruit_taste <- ggordpoint(obj=pcoares_Fruit_taste, biplot=F, speciesannot=F,
                                   factorNames=c("Fruit_taste"), ellipse=F) +
  scale_fill_manual(values=myco_taste)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "D");pcoaplot_Fruit_taste

###################################################################################################
###################################################################################################
Bacteria_yield <- subset_samples(Bacteria, !is.na(yield));Bacteria_yield
distme_yield <- get_dist(Bacteria_yield, distmethod ="bray", method="hellinger");distme_yield
sampleda_yield <- data.frame(sample_data(Bacteria_yield), check.names=FALSE);sampleda_yield
sampleda_yield <- sampleda_yield[match(colnames(as.matrix(distme_yield)),rownames(sampleda_yield)),,drop=FALSE];sampleda_yield
sampleda_yield$yield <- factor(sampleda_yield$yield)
adores <- adonis2(distme_yield ~ yield, data=sampleda_yield, permutation=999);adores#change from adonis to adonis2
#           Df SumOfSqs     R2      F Pr(>F)    
#yield       4    8.146 0.0293 7.8241  0.001 ***
#Residual 1037  269.921 0.9707                  
#Total    1041  278.067 1.0000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_yield=c("chartreuse4","yellow3","#0B775E","#00AFBB","grey35");myco_yield
pcoares_yield <- get_pcoa(obj=Bacteria_yield, distmethod="bray", method="hellinger");pcoares_yield#pcoa
pcoaplot_yield <- ggordpoint(obj=pcoares_yield, biplot=F, speciesannot=F,
                             factorNames=c("yield"), ellipse=F) +
  scale_fill_manual(values=myco_yield)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "E");pcoaplot_yield

###################################################################################################
###################################################################################################
Bacteria_seedweight_1000 <- subset_samples(Bacteria, !is.na(seedweight_1000));Bacteria_seedweight_1000#
distme_seedweight_1000 <- get_dist(Bacteria_seedweight_1000, distmethod ="bray", method="hellinger");distme_seedweight_1000
sampleda_seedweight_1000 <- data.frame(sample_data(Bacteria_seedweight_1000), check.names=FALSE);sampleda_seedweight_1000
sampleda_seedweight_1000 <- sampleda_seedweight_1000[match(colnames(as.matrix(distme_seedweight_1000)),rownames(sampleda_seedweight_1000)),,drop=FALSE];sampleda_seedweight_1000
sampleda_seedweight_1000$seedweight_1000 <- factor(sampleda_seedweight_1000$seedweight_1000)
adores <- adonis2(distme_seedweight_1000 ~ seedweight_1000, data=sampleda_seedweight_1000, permutation=999);adores#change from adonis to adonis2
#                 Df SumOfSqs     R2     F Pr(>F)    
#seedweight_1000   3    6.675 0.0296 8.593  0.001 ***
#Residual        845  218.811 0.9704                 
#Total           848  225.487 1.0000                 
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##############################################################################################################################################################################################################
myco_weight=c("#0B775E","#FB9A99","yellow3","grey35");myco_weight
pcoares_seedweight_1000 <- get_pcoa(obj=Bacteria_seedweight_1000, distmethod="bray", method="hellinger");pcoares_seedweight_1000#pcoa
pcoaplot_seedweight_1000 <- ggordpoint(obj=pcoares_seedweight_1000, biplot=F, speciesannot=F,
                                       factorNames=c("seedweight_1000"), ellipse=F) +
  scale_fill_manual(values=myco_weight)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "F");pcoaplot_seedweight_1000#

###################################################################################################
###################################################################################################
Bacteria_Fruit_hori_diameter <- subset_samples(Bacteria, !is.na(Fruit_hori_diameter));Bacteria_Fruit_hori_diameter
distme_Fruit_hori_diameter <- get_dist(Bacteria_Fruit_hori_diameter, distmethod ="bray", method="hellinger");distme_Fruit_hori_diameter
sampleda_Fruit_hori_diameter <- data.frame(sample_data(Bacteria_Fruit_hori_diameter), check.names=FALSE);sampleda_Fruit_hori_diameter
sampleda_Fruit_hori_diameter <- sampleda_Fruit_hori_diameter[match(colnames(as.matrix(distme_Fruit_hori_diameter)),rownames(sampleda_Fruit_hori_diameter)),,drop=FALSE];sampleda_Fruit_hori_diameter
sampleda_Fruit_hori_diameter$Fruit_hori_diameter <- factor(sampleda_Fruit_hori_diameter$Fruit_hori_diameter)
adores <- adonis2(distme_Fruit_hori_diameter ~ Fruit_hori_diameter, data=sampleda_Fruit_hori_diameter, permutation=999);adores#change from adonis to adonis2
#                       Df SumOfSqs     R2      F Pr(>F)    
#Fruit_hori_diameter    2    5.101 0.0174 9.6611  0.001 ***
#Residual            1091  288.020 0.9826                  
#Total               1093  293.121 1.0000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_horiz=c("chartreuse4","#00AFBB","#FB9A99");myco_horiz
pcoares_Fruit_hori_diameter <- get_pcoa(obj=Bacteria_Fruit_hori_diameter, distmethod="bray", method="hellinger");pcoares_Fruit_hori_diameter#pcoa
pcoaplot_Fruit_hori_diameter <- ggordpoint(obj=pcoares_Fruit_hori_diameter, biplot=F, speciesannot=F,
                                           factorNames=c("Fruit_hori_diameter"), ellipse=F) +
  scale_fill_manual(values=myco_horiz)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "A");pcoaplot_Fruit_hori_diameter

###################################################################################################
###################################################################################################
Bacteria_Fruit_long_diameter <- subset_samples(Bacteria, !is.na(Fruit_long_diameter));Bacteria_Fruit_long_diameter
distme_Fruit_long_diameter <- get_dist(Bacteria_Fruit_long_diameter, distmethod ="bray", method="hellinger");distme_Fruit_long_diameter
sampleda_Fruit_long_diameter <- data.frame(sample_data(Bacteria_Fruit_long_diameter), check.names=FALSE);sampleda_Fruit_long_diameter
sampleda_Fruit_long_diameter <- sampleda_Fruit_long_diameter[match(colnames(as.matrix(distme_Fruit_long_diameter)),rownames(sampleda_Fruit_long_diameter)),,drop=FALSE];sampleda_Fruit_long_diameter
sampleda_Fruit_long_diameter$Fruit_long_diameter <- factor(sampleda_Fruit_long_diameter$Fruit_long_diameter)
adores <- adonis2(distme_Fruit_long_diameter ~ Fruit_long_diameter, data=sampleda_Fruit_long_diameter, permutation=999);adores#change from adonis to adonis2
#                       Df SumOfSqs      R2      F Pr(>F)    
#Fruit_long_diameter    2     5.43 0.01882 10.347  0.001 ***
#Residual            1079   283.12 0.98118                  
#Total               1081   288.55 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_longdiam=c("chartreuse4","#00AFBB","#FB9A99" );myco_longdiam
pcoares_Fruit_long_diameter <- get_pcoa(obj=Bacteria_Fruit_long_diameter, distmethod="bray", method="hellinger");pcoares_Fruit_long_diameter#pcoa
pcoaplot_Fruit_long_diameter <- ggordpoint(obj=pcoares_Fruit_long_diameter, biplot=F, speciesannot=F,
                                           factorNames=c("Fruit_long_diameter"), ellipse=F) +
  scale_fill_manual(values=myco_longdiam)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "B");pcoaplot_Fruit_long_diameter

###################################################################################################
###################################################################################################
#TMV resisitance
Bacteria_TMV_resistance <- subset_samples(Bacteria, !is.na(TMV_resistance));Bacteria_TMV_resistance
distme_TmvRes <- get_dist(Bacteria_TMV_resistance, distmethod ="bray", method="hellinger");distme_TmvRes
sampleda_TmvRes <- data.frame(sample_data(Bacteria_TMV_resistance), check.names=FALSE);sampleda_TmvRes
sampleda_TmvRes <- sampleda_TmvRes[match(colnames(as.matrix(distme_TmvRes)),rownames(sampleda_TmvRes)),,drop=FALSE];sampleda_TmvRes
sampleda_TmvRes$TMV_resistance <- factor(sampleda_TmvRes$TMV_resistance)
adores <- adonis2(distme_TmvRes ~ TMV_resistance, data=sampleda_TmvRes, permutation=999);adores#change from adonis to adonis2
#                 Df SumOfSqs      R2      F Pr(>F)    
#TMV_resistance   2    2.327 0.01441 4.5263  0.001 ***
#Residual       619  159.134 0.98559                  
#Total          621  161.461 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_TMV=c("#0B775E","#FB9A99","grey35");myco_TMV
pcoares_TmvRes <- get_pcoa(obj=Bacteria_TMV_resistance, distmethod="bray", method="hellinger");pcoares_TmvRes#pcoa
pcoaplot_TmvRes <- ggordpoint(obj=pcoares_TmvRes, biplot=F, speciesannot=F,
                              factorNames=c("TMV_resistance"), ellipse=F) +
  scale_fill_manual(values=myco_TMV)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=1))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "C");pcoaplot_TmvRes

##############################################################################################################################################################################################################
##############################################################################################################################################################################################################
Bacteria_Insect_resistance <- subset_samples(Bacteria, !is.na(Insect_resistance));Bacteria_Insect_resistance
distme_InsectRes <- get_dist(Bacteria_Insect_resistance, distmethod ="bray", method="hellinger");distme_InsectRes
sampleda_InsectRes <- data.frame(sample_data(Bacteria_Insect_resistance), check.names=FALSE);sampleda_InsectRes
sampleda_InsectRes <- sampleda_InsectRes[match(colnames(as.matrix(distme_InsectRes)),rownames(sampleda_InsectRes)),,drop=FALSE];sampleda_InsectRes
sampleda_InsectRes$Insect_resistance <- factor(sampleda_InsectRes$Insect_resistance)
adores <- adonis2(distme_InsectRes ~ Insect_resistance, data=sampleda_InsectRes, permutation=999);adores#change from adonis to adonis2
#                   Df SumOfSqs      R2      F Pr(>F)    
#Insect_resistance   5   10.617 0.06632 8.4097  0.001 ***
#Residual          592  149.474 0.93368                  
#Total             597  160.091 1.00000                  
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##############################################################################################################################################################################################################
myco_insRes=c("yellow3","grey35","chartreuse4","#0B775E","#00AFBB", "#FB9A99");myco_insRes
pcoares_InsectRes <- get_pcoa(obj=Bacteria_Insect_resistance, distmethod="bray", method="hellinger");pcoares_InsectRes#pcoa
pcoaplot_InsectRes <- ggordpoint(obj=pcoares_InsectRes, biplot=F, speciesannot=F,
                              factorNames=c("Insect_resistance"), ellipse=F) +
  scale_fill_manual(values=myco_insRes)+
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=2))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "bottom")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=12,  family="sans"))+ labs(tag = "D");pcoaplot_InsectRes

###################################################################################################
###################################################################################################
library(patchwork)
tiff("PCOA_Parameters.tiff", units="in", width=9, height=12, res=300)
PCOA_Parameters.tiff=pcoaplot_ovaries+pcoaplot_Fshape+pcoaplot_Fcol+pcoaplot_Fruit_taste+pcoaplot_yield+pcoaplot_seedweight_1000+
  plot_layout(ncol = 2);PCOA_Parameters.tiff
ggsave("PCOA_Parameters.tiff", plot = PCOA_Parameters.tiff)
dev.off()
#####################################################################################################################################################################################################
tiff("PCOA_Parameters2.tiff", units="in", width=11, height=11, res=300)
PCOA_Parameters2.tiff=pcoaplot_Fruit_hori_diameter+pcoaplot_Fruit_long_diameter+pcoaplot_TmvRes+pcoaplot_InsectRes+plot_layout(ncol = 2);PCOA_Parameters2.tiff
ggsave("PCOA_Parameters2.tiff", plot = PCOA_Parameters2.tiff)
dev.off()

###################################################################################################
###################################################################################################
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
#Good representation
library(ggplot2)
library(MicrobiotaProcess)
citation("MicrobiotaProcess")
classtaxa <- get_taxadf(obj=Bacteria, taxlevel=3);classtaxa
phylataxa <- get_taxadf(obj=Bacteria, taxlevel=2);phylataxa
ordertaxa <- get_taxadf(obj=Bacteria, taxlevel=4);ordertaxa
Famtaxa <- get_taxadf(obj=Bacteria, taxlevel=5);Famtaxa
Gentaxa <- get_taxadf(obj=Bacteria, taxlevel=6);Gentaxa
###########
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
tiff("fphyla.tiff", units="in", width=12, height=20, res=300)
fphyla <- ggbartax(obj=phylataxa, facetNames="English_Name_2", plotgroup=TRUE, topn=15) +
  xlab("Cultivars") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 270)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11))+ labs(tag = "a")+
  theme(legend.position = "top")+ 
  labs(fill='Phyla')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"));fphyla

#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
## Show the abundance in different groups considering cultivar summary.
#order
tiff("fclass.tiff", units="in", width=12, height=20, res=300)
fclass <- ggbartax(obj=classtaxa, facetNames="English_Name_2", plotgroup=TRUE, topn=20) +
  xlab("Cultivars") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=3))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11))+ labs(tag = "a")+
  theme(legend.position = "bottom")+ 
  labs(fill='Classes')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();fclass
####################
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
#order
tiff("forder.tiff", units="in", width=12, height=20, res=300)
forder <- ggbartax(obj=ordertaxa, facetNames="English_Name_2", plotgroup=TRUE, topn=30) +
  xlab("Cultivars") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=3))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11))+ labs(tag = "a")+
  theme(legend.position = "bottom")+ 
  labs(fill='Orders')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();forder
dev.off()
########################
#https://yulab-smu.top/MicrobiotaProcessWorkshop/articles/MicrobiotaProcessWorkshop.html
#Family
tiff("ffam.tiff", units="in", width=12, height=20, res=300)
ffam <- ggbartax(obj=Famtaxa, facetNames="English_Name_2", plotgroup=TRUE, topn=30) +
  xlab("Cultivars") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=3))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11))+ labs(tag = "a")+
  theme(legend.position = "bottom")+ 
  labs(fill='Family')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam
dev.off()
########################
#Genus
tiff("fgen.tiff", units="in", width=12, height=20, res=300)
fgen <- ggbartax(obj=Gentaxa, facetNames="English_Name_2", plotgroup=TRUE, topn=30) +
  xlab("Cultivars") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=3))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=11),
        axis.text.y  = element_text(colour="black", size=11))+ labs(tag = "a")+
  theme(legend.position = "bottom")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();fgen
dev.off()

###################################################################################################
###################################################################################################
#SPECIFIC PARAMETERS
#ovaries
ffam_ovaries <- ggbartax(obj=Famtaxa, facetNames="ovaries", plotgroup=TRUE, topn=30) +
  xlab("No. of Ovaries") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "A")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_ovaries
#############################
#Fruit_shape
ffam_Fruit_shape <- ggbartax(obj=Famtaxa, facetNames="Fruit_shape", plotgroup=TRUE, topn=30) +
  xlab("Berry Shape") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "B")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Fruit_shape
#############################
#Fruit_taste

ffam_Fruit_taste <- ggbartax(obj=Famtaxa, facetNames="Fruit_taste", plotgroup=TRUE, topn=30) +
  xlab("Berry Taste") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "C")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Fruit_taste
#############################
#Fruit_color

ffam_Fruit_color <- ggbartax(obj=Famtaxa, facetNames="Fruit_color", plotgroup=TRUE, topn=30) +
  xlab("Berry Color") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "D")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Fruit_color

#############################
#Yield
ffam_yield <- ggbartax(obj=Famtaxa, facetNames="yield", plotgroup=TRUE, topn=30) +
  xlab("Yield (kg/sq.meter)") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "E")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_yield
#############################
#seedweight_1000
ffam_seedweight_1000 <- ggbartax(obj=Famtaxa, facetNames="seedweight_1000", plotgroup=TRUE, topn=30) +
  xlab("weight of 1000 seeds/g") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "F")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_seedweight_1000
###############################
#Fruit_hori_diameter
ffam_Fruit_hori_diameter <- ggbartax(obj=Famtaxa, facetNames="Fruit_hori_diameter", plotgroup=TRUE, topn=30) +
  xlab("Hor. Diam./cm") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "G")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Fruit_hori_diameter
###############################
#Fruit_long_diameter
ffam_Fruit_long_diameter <- ggbartax(obj=Famtaxa, facetNames="Fruit_long_diameter", plotgroup=TRUE, topn=30) +
  xlab("Long. Diam./cm") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "H")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Fruit_long_diameter
###############################
#TMV_resistance
ffam_TMV_resistance <- ggbartax(obj=Famtaxa, facetNames="TMV_resistance", plotgroup=TRUE, topn=30) +
  xlab("TMV Res.") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=5))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "I")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_TMV_resistance
###############################
#Insect_resistance
ffam_Insect_resistance <- ggbartax(obj=Famtaxa, facetNames="Insect_resistance", plotgroup=TRUE, topn=30) +
  xlab("Insect Res.") +
  ylab("relative abundance (%)") +scale_fill_manual(values=mycolors)+theme(axis.text.x = element_text(angle = 90)) +
  #scale_color_manual(values=c("chartreuse4","bisque4","deepskyblue","Sky Blue","goldenrod", "#ba0d42", "grey85","Orange","coral3","darkred","#0dba2f","#850dba","grey35","Blue","deeppink","blue","Yellow")) + 
  guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, ncol=4))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=0.4))+ 
  theme(axis.title.x = element_text(face="bold", colour="black", size=14),
        axis.title.y = element_text(face="bold", colour="black", size=14),
        axis.text.x = element_text(colour="black", size=12),
        axis.text.y  = element_text(colour="black", size=12))+ labs(tag = "J")+
  theme(legend.position = "")+ 
  labs(fill='Genera')+
  theme(legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12))+
  theme(text=element_text(family="sans"))+ coord_flip();ffam_Insect_resistance
########################################################################################################################################################################################
tiff("Barplot_Parameters.tiff", units="in", width=15, height=10, res=300)
Barplot_Parameters.tiff=ffam_ovaries+ffam_Fruit_shape+ffam_Fruit_taste+ffam_Fruit_color+ffam_yield+ffam_seedweight_1000+ffam_Fruit_hori_diameter+ffam_Fruit_long_diameter+ ffam_TMV_resistance+
  ffam_Insect_resistance+plot_layout(ncol = 3);Barplot_Parameters.tiff
ggsave("Barplot_Parameters.tiff", plot = Barplot_Parameters.tiff)
dev.off()
########################################################################################################################################################################################
ffam_Insect_resistance+
###################################################################################################
###################################################################################################
#https://zachcp.github.io/phylogeo/projections.html
#Spatial reprwsentation of the data (phylogeo)
#devtools::install_github("zachcp/phylogeo")
#install.packages("rmarkdown")
library(phylogeo)
library(ggplot2)
library(gridExtra)
library(phylogeo)
library(rmarkdown)
#render("Tomato100seeds.R")#make markdown
map=map_phyloseq(Bacteria);map
sample_data(Bacteria)
#zoomed in map, colored by pH and jittered in order to see the points better
map_phyloseq(Bacteria, region="China", jitter=TRUE, jitter.x=10,jitter.y=10, color="Production_site_2_2" )

htmlmap_phyloseq(Bacteria)

# customize maps with leaflet
m <- htmlmap_phyloseq(Bacteria);m
#simple network map without lines 
map_network(Bacteria)
htmlmap_network(Bacteria)
map_tree(Bacteria)
map_clusters(Bacteria, clusternum=4)

################################
east_asia <- map_data("world", region = c("China"));east_asia
write_csv(east_asia,"east_asia.csv")
East_asia <- read.table ("east_asia.txt", check.names = FALSE, header = TRUE, sep = "\t");East_asia

# Map region to fill color
ggplot(east_asia, aes(x = long, y = lat, group = subregion, fill = subregion)) +
  geom_polygon(colour = "black") +
  scale_colour_manual(values = mycola)+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())
###################################################################################################
###################################################################################################

#adopted from https://microbiome.github.io/tutorials/Composition.html
Bacteria.com <- Bacteria
source("aggregate_top_taxa.txt")
###########
Bacteria.com.phy <- aggregate_top_taxa(Bacteria.com, "Phylum", top = 5)
# Use traqnsform function of microbiome to convert it to rel abun.
Bacteria.com.phy.rel <- microbiome::transform(Bacteria.com.phy, "compositional")
plot.composition.relAbun <- plot_composition(Bacteria.com.phy.rel,
                                             sample.sort = "Firmicutes",
                                             x.label = "Samples");plot.composition.relAbun
plot.composition.relAbun <- plot.composition.relAbun + theme(legend.position = "bottom");plot.composition.relAbun 
plot.composition.relAbun <- plot.composition.relAbun + scale_fill_brewer("Phylum", palette = "Paired") + theme_bw();plot.composition.relAbun 
plot.composition.relAbun <- plot.composition.relAbun + theme(axis.text.x = element_text(angle = 90));plot.composition.relAbun
plot.composition.relAbun <- plot.composition.relAbun + ggtitle("Relative abundance") + theme(legend.title = element_text(size = 18))+
  scale_colour_manual(values = mycola );plot.composition.relAbun

print(plot.composition.relAbun)
#################################################################################################################################################################################
#################################################################################################################################################################################
#good representation of taxa
#Following is an example of customizing the plot using ggpubr.
# we use count data at family level from the barplot for counts
ps_df <- microbiomeutilities::phy_to_ldf(Bacteria.com.phy, 
                                         transform.counts = "compositional")
colnames(ps_df)
# this data.frame can be used to customize several plots.  
# example boxplot at family level
p.box.tiff <- ggstripchart(ps_df, "Cultivar_summary", "Abundance", 
                      facet.by = "Phylum", color = "Phylum",
                      palette = "jco")+rremove("x.text")+scale_colour_manual(values = mycola)+theme(legend.position = "top")+
  theme(axis.title.x = element_text(face="bold", colour="black", size=14,family="sans"),
        axis.title.y = element_text(face="bold", colour="black", size=14,family="sans"),
        axis.text.y  = element_text(face="bold", colour="black", size=12,family="sans"))+ labs(tag = "A")+ 
  xlab("Cultivar") +
  ylab("Abundance");p.box.tiff

tiff("p.box.tiff", units="in", width=10, height=5, res=300)
p.box.tiff=p.box.tiff+plot_layout(ncol = 1);p.box.tiff
ggsave("p.box.tiff", plot = p.box.tiff)
dev.off()
#################################################################################################################################################################################
#################################################################################################################################################################################
p.box1.tiff <- ggstripchart(ps_df, "Cultivar_summary", "Abundance", 
                      facet.by = "Phylum", color = "Production_site_2",
                      palette = "jco")+ rremove("x.text")+scale_colour_manual(values = mycola)+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.title.y = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.text.y  = element_text(face="bold", colour="gray33", size=12,family="sans"))+ labs(tag = "E")+ 
  xlab("Cultivar") +
  ylab("Abundance")+ guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=2))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "top")


tiff("p.box1.tiff", units="in", width=10, height=6, res=300)
p.box1.tiff=p.box1.tiff+plot_layout(ncol = 1);p.box1.tiff
ggsave("p.box1.tiff", plot = p.box1.tiff)
dev.off()

write_csv(ps_df,"ps_df_relative abundances_phylum.csv")
#################################################################################################################################################################################
#################################################################################################################################################################################
#Consider family level
Bacteria.com <- Bacteria
Bacteria.com.fam <- aggregate_top_taxa(Bacteria.com, "Family", top = 11)
Bacteria.com.fam.rel <- microbiome::transform(Bacteria.com.fam, "compositional")
#################################
ps_df <- microbiomeutilities::phy_to_ldf(Bacteria.com.fam.rel, 
                                         transform.counts = "compositional")
colnames(ps_df)
# this data.frame can be used to customize several plots.  
# example boxplot at family level
p.box2.tiff <- ggstripchart(ps_df, "Cultivar_summary", "Abundance", 
                      facet.by = "Family", color = "Family",
                      palette = "jco")+ rremove("x.text")+scale_colour_manual(values = mycola)+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.title.y = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.text.y  = element_text(face="bold", colour="gray33", size=12,family="sans"))+ labs(tag = "E")+ 
  xlab("Cultivar") +
  ylab("Abundance")+ guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=2))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "top")

tiff("p.box2.tiff", units="in", width=10, height=6, res=300)
p.box2.tiff=p.box2.tiff+plot_layout(ncol = 1);p.box2.tiff
ggsave("p.box2.tiff", plot = p.box2.tiff)
dev.off()

#################################################################################################################################################################################
#################################################################################################################################################################################
ps_df <- microbiomeutilities::phy_to_ldf(Bacteria.com.fam.rel, transform.counts = "compositional")
colnames(ps_df)
# this data.frame can be used to customize several plots.  
# example boxplot at family level
p.box3.tiff <- ggstripchart(ps_df, "Cultivar_summary", "Abundance", 
                       facet.by = "Family", color = "Production_site_2",
                       palette = "jco") + rremove("x.text")+scale_colour_manual(values = mycola)+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.title.y = element_text(face="bold", colour="gray33", size=14,family="sans"),
        axis.text.y  = element_text(face="bold", colour="gray33", size=12,family="sans"))+ labs(tag = "E")+ 
  xlab("Cultivar") +
  ylab("Abundance")+ guides(fill= guide_legend(keywidth = 0.6, keyheight = 0.6, nrow=2))+
  theme(axis.text.x = element_text(size=12, face="bold", color = "black"),
        axis.text.y = element_text(size=12, face="bold", color = "black"))+
  theme(text=element_text(size=12,  family="sans"))+theme(legend.position = "")

tiff("p.box3.tiff", units="in", width=10, height=6, res=300)
p.box3.tiff=p.box3.tiff+plot_layout(ncol = 1);p.box3.tiff
ggsave("p.box3.tiff", plot = p.box3.tiff)
dev.off()

write_csv(ps_df,"ps_df_relative abundances_Family.csv")

###################################################################################################
###################################################################################################
#Plotting a ggstripchart
#Following is an example of customizing the plot using ggpubr.
# we use count data at family level from the barplot for counts
ps_df <- microbiomeutilities::phy_to_ldf(Bacteria.com.fam, 
                                         transform.counts = "compositional")
colnames(ps_df)
# this data.frame can be used to customize several plots.  
# example boxplot at family level
p.box <- ggstripchart(ps_df, "Cultivar_summary", "Abundance", 
                      facet.by = "Family", color = "Family",
                      palette = "jco") 
p.box + rremove("x.text")+scale_colour_manual(values =mycol)
###################################################################################################
###################################################################################################
#https://uw-madison-microbiome-hub.github.io/Microbiome_analysis_in-_R/ (good)
hyseq_fam <- microbiome::aggregate_rare(Bacteria, level = "Family", detection = 50/100, prevalence = 70/100)
# Alternative
Bacteria_fam <- aggregate_top_taxa2(Bacteria, level = "Family", top = 25)
Bacteria.fam.rel <- microbiome::transform(Bacteria_fam, "compositional")
Bacteria.fam.rel <- Bacteria %>%
  aggregate_rare(level = "Family", detection = 50/100, prevalence = 70/100) %>%
  microbiome::transform(transform = "compositional")

########################################################################################################################################################################
########################################################################################################################################################################
#Alpha diversity measures
########################################################################################################################################################################
library(microViz)
library(ggtext)
library(microViz)
library(ggraph)
library(DT)
library(corncob)
library(ggplot2)
library(patchwork)
library(ggplot2) # 3.2.1
library(dplyr)
library(lubridate)
library(cowplot) # 1.0.0
library(egg) # 0.4.5
#################################################################################################################################################################################
#################################################################################################################################################################################
#Rarefying the dataset
library(microbiome)
library(knitr)
library(tidyplots)
summary(sample_sums(Bacteria))#[ 15421 taxa and 1197 samples ]
Bacteria_subsamp <- rarefy_even_depth(Bacteria, sample.size=500 , rngseed=5163, replace=FALSE, trimOTUs=TRUE);Bacteria_subsamp#[ 1128 samples by 20 sample variables ]

#lost Samples after rarefaction: 69 samples
#Diversity was computed for using the functions mentioned in Rhea package considering ASVs with abundance >0.25%
#Load the alpha diversity data (rarefied at 500 depth)#alpha-diversity_500.tab
#condider the phylogenetic diversity plot
#pd=estimate_pd(Bacteria_subsamp);pd
#write.csv(pd, file="pd.csv")
#Combime the pg onto the data frame for alha at 500 depth (alpha_500_PD)
alpha_500_PD <- read.table ("alpha-diversity_500.tab.txt", check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "");alpha_500_PD
alpha_500_PD$yield <- cut(alpha_500_PD$`Yield_(kg/sq_meter)`, breaks = seq(0,150,20), include.lowest = TRUE)
alpha_500_PD$yield
alpha_500_PD$Fruit_long_diameter=cut(alpha_500_PD$Fruit_longitudinal_diamete, breaks = seq(0,10,3), include.lowest = TRUE)
alpha_500_PD$Fruit_long_diameter
alpha_500_PD$Fruit_hori_diameter=cut(alpha_500_PD$Fruit_horizontal_diameter, breaks = seq(0,10,3), include.lowest = TRUE)
alpha_500_PD$Fruit_hori_diameter
alpha_500_PD$seedweight_1000=cut(alpha_500_PD$`1000_seeds_weight_(units)`, breaks = seq(0,5,1), include.lowest = TRUE)
alpha_500_PD$seedweight_1000
alpha_500_PD$ovaries =cut(alpha_500_PD$Number_of_ovary, breaks = seq(0,12,4), include.lowest = TRUE)
alpha_500_PD$ovaries
#write_csv(alpha_500_PD,"alpha_500_PD.csv")
#subset for the different parameters
alpha_500_PD_Fruit_taste=alpha_500_PD[!is.na(alpha_500_PD$Fruit_taste),];alpha_500_PD_Fruit_taste
alpha_500_PD_Fruit_color=alpha_500_PD[!is.na(alpha_500_PD$Fruit_color),];alpha_500_PD_Fruit_color
alpha_500_PD_Fruit_shape=alpha_500_PD[!is.na(alpha_500_PD$Fruit_shape),];alpha_500_PD_Fruit_shape
alpha_500_PD_TMV_resistance=alpha_500_PD[!is.na(alpha_500_PD$TMV_resistance),];alpha_500_PD_TMV_resistance
alpha_500_PD_Insect_resistance=alpha_500_PD[!is.na(alpha_500_PD$Insect_resistance),];alpha_500_PD_Insect_resistance
alpha_500_PD_yield=alpha_500_PD[!is.na(alpha_500_PD$yield),];alpha_500_PD_yield
alpha_500_PD_seedweight_1000=alpha_500_PD[!is.na(alpha_500_PD$seedweight_1000),];alpha_500_PD_seedweight_1000
attach(alpha_500_PD)
alpha_500_PD$Insect_resistance_2
#RICHNESS considering the rarefied object
variable.names(alpha_500_PD)
#[1] "PD"                         "Richness"                   "Normalized_Richness"       
#[4] "Effective_Richness"         "Shannon_Index"              "Shannon_Effective"         
#[7] "Simpson_Index"              "Simpson_Effective"          "Evenness"                  
#[10] "rarefaction_depth"          "Cultivar_summary"           "English_Name_2"              
#[13] "Country_of_Origin"          "Production_site_2"            "latitude"                  
#[16] "longitude"                  "Production_year"            "1000_seeds_weight_(units)" 
#[19] "TMV_resistance"             "Insect_resistance"          "Fruit_longitudinal_diamete"
#[22] "Fruit_horizontal_diameter"  "Fruit_shape"                "Fruit_color"               
#[25] "Fruit_taste"                "Number_of_ovary"            "Sowing_date"               
#[28] "Begining_of_Harvest"        "End_of_Harvest"             "Yield_(g/sq_meter)" 
########################################################################################################################library(rstatix)
#citation("rstatix ")
library(rcompanion)
library(multcompView)
library(tidyverse)
library(ggpubr)
library(dunn.test)
library(FSA)
#install.packages("explore")
library(explore)
explore(alpha_500_PD)
alpha_500_PD|> explore()
alpha_500_PD |> describe()
alpha_500_PD |> explore(type)
alpha_500_PD |> explore(Production_year)
alpha_500_PD |> explore(energy_kcal_100ml, target = type)
alpha_500_PD |> explore(alcohol_vol_pct, energy_kcal_100ml, target = type)
# report of all variables
alpha_500_PD %>% report(output_file = "report.html", output_dir = tempdir())

alpha_500_PD %>% get_summary_stats(type = "common")
#variable                       n      min       max   median    iqr     mean       sd      se      ci
#<fct>                      <dbl>    <dbl>     <dbl>    <dbl>  <dbl>    <dbl>    <dbl>   <dbl>   <dbl>
#  1 PD                          1128    1.08     13.4      3.90   2.58     4.21     1.87    0.056   0.109
#2 Richness                    1128    6       164       36     28       41.0     21.5     0.641   1.26 
#3 Normalized_Richness         1128    6       164       36     28       41.0     21.5     0.641   1.26 
#4 Effective_Richness          1128    3        85       24     20       27.6     14.5     0.432   0.848
#5 Shannon_Index               1128    0.11      4.58     2.33   1.19     2.27     0.892   0.027   0.052
#6 Shannon_Effective           1128    1.12     97.2     10.2   12.8     13.9     11.9     0.354   0.694
#7 Simpson_Index               1128    0.017     0.968    0.175  0.236    0.26     0.226   0.007   0.013
#8 Simpson_Effective           1128    1.03     58.0      5.70   6.98     7.89     7.24    0.216   0.423
#9 Evenness                    1128    0.039     0.638    0.457  0.159    0.427    0.132   0.004   0.008
#10 rarefaction_depth           1128  500       500      500      0      500        0       0       0    
#11 latitude                    1128   23.6      48.5     32.9    6.41    35.1      6.55    0.195   0.383
#12 longitude                   1128  102.      128.     119.     2.74   117.       5.99    0.178   0.35 
#13 Production_year             1128 2019      2019     2019      0     2019        0       0       0    
#14 1000_seeds_weight_(units)    849    1.6       5        3.2    0.4      3.16     0.489   0.017   0.033
#15 Fruit_longitudinal_diamete  1104    1.8   45448        4.8    1.3    417.    4307.    130.    254.   
#16 Fruit_horizontal_diameter   1104    2     45479        5.3    2.3    417.    4310.    130.    255.   
#17 Number_of_ovary             1112    2        10        4      3        4.23     2.04    0.061   0.12 
#18 Yield_(kg/sq_meter)          996    8.78    100       45.4   28.5     44.9     19.0     0.6     1.18 
library(tidyverse)
data <- readr::read_delim('https://www.nxn.se/single-cell-studies/data.tsv', delim = '\t')
#STATISTICS FOR ALPHA DIVERSITY
#https://www.rdocumentation.org/packages/rstatix/versions/0.7.2

# RICHNESS
kruskal.test(Richness ~ English_Name_2, data = alpha_500_PD)
#1 Kruskal-Wallis chi-squared = 761.57, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Richness ~ Production_site_2, data = alpha_500_PD)
#data:  Richness by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_richness_Production_site_2=dunn_test(alpha_500_PD, Richness ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_richness_Production_site_2
write_csv(Dunn_tomato_richness_Production_site_2, "Dunn_tomato_richness_Production_site_2.csv")
Dunn_tomato_richness_Production_site_2 <- dunnTest(alpha_500_PD$Richness, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_richness_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_richness_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_richness_Production_site_2$res);CLD_Dunn_tomato_richness_Production_site_2 

#Dunn_tomato_richness_bonferon=dunn_test(alpha_500_PD, Richness ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_richness_bonferon
#write_csv(Dunn_tomato_richness_bonferon, "Dunn_tomato_richness_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_richness_English_Name_2=dunn_test(alpha_500_PD, Richness ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_richness_English_Name_2
write_csv(Dunn_tomato_richness_English_Name_2, "Dunn_tomato_richness_English_Name_2.csv")
Dunn_tomato_richness_English_Name_2 <- dunnTest(alpha_500_PD$Richness, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_richness_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_richness_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_richness_English_Name_2$res);CLD_Dunn_tomato_richness_English_Name_2 


kruskal.test(Richness ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Richness  1036      6.57     1 0.0104 Kruskal-Wallis

kruskal.test(Richness ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Richness  1116      19.7     4 0.000583 Kruskal-Walliss

kruskal.test(Richness ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Richness  1116      38.6     2 0.00000000408 Kruskal-Wallis

kruskal.test(Richness ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Richness   610      12.4     2 0.00207 Kruskal-Wallis

kruskal.test(Richness ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Richness   566      24.5     5 0.000173 Kruskal-Wallis
kruskal.test(Richness ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Richness   849      3.89     3 0.273 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Richness_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Richness ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Richness_Fruit_color
write_csv(Dunn_tomato_Richness_Fruit_color, "Dunn_tomato_Richness_Fruit_color.csv")
Dunn_tomato_Richness_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Richness, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Richness_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Richness_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Richness_Fruit_color$res);CLD_Dunn_tomato_Richness_Fruit_color 

#insect resistance
Dunn_tomato_Richness_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Richness ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Richness_Insect_resistance
write_csv(Dunn_tomato_Richness_Insect_resistance, "Dunn_tomato_Richness_Insect_resistance.csv")
Dunn_tomato_Richness_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Richness, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Richness_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Richness_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Richness_Insect_resistance$res);CLD_Dunn_tomato_Richness_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Richness_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Richness ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Richness_TMV_resistance
write_csv(Dunn_tomato_Richness_TMV_resistance, "Dunn_tomato_Richness_TMV_resistance.csv")
Dunn_tomato_Richness_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Richness, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Richness_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Richness_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Richness_TMV_resistance$res);CLD_Dunn_tomato_Richness_TMV_resistance 

#seed weoght
Dunn_tomato_Richness_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Richness ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Richness_seedweight_1000
write_csv(Dunn_tomato_Richness_seedweight_1000, "Dunn_tomato_Richness_seedweight_1000.csv")
Dunn_tomato_Richness_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Richness, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Richness_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Richness_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Richness_seedweight_1000$res);CLD_Dunn_tomato_Richness_seedweight_1000 


#seed weoght
Dunn_tomato_Richness_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Richness ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Richness_Fruit_shape
write_csv(Dunn_tomato_Richness_Fruit_shape, "Dunn_tomato_Richness_Fruit_shape.csv")
Dunn_tomato_Richness_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Richness, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Richness_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Richness_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Richness_Fruit_shape$res);CLD_Dunn_tomato_Richness_Fruit_shape 

###########################################################################################
#EFFECTIVE RICHNESS
kruskal.test(Effective_Richness ~ English_Name_2, data = alpha_500_PD)
#1 Kruskal-Wallis chi-squared = 761.57, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Effective_Richness ~ Production_site_2, data = alpha_500_PD)
#data:  Effective_Richness by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Effective_Richness_Production_site_2=dunn_test(alpha_500_PD, Effective_Richness ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_Production_site_2
write_csv(Dunn_tomato_Effective_Richness_Production_site_2, "Dunn_tomato_Effective_Richness_Production_site_2.csv")
Dunn_tomato_Effective_Richness_Production_site_2 <- dunnTest(alpha_500_PD$Effective_Richness, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Effective_Richness_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_Production_site_2$res);CLD_Dunn_tomato_Effective_Richness_Production_site_2 

#Dunn_tomato_Effective_Richness_bonferon=dunn_test(alpha_500_PD, Effective_Richness ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Effective_Richness_bonferon
#write_csv(Dunn_tomato_Effective_Richness_bonferon, "Dunn_tomato_Effective_Richness_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Effective_Richness_English_Name_2=dunn_test(alpha_500_PD, Effective_Richness ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_English_Name_2
write_csv(Dunn_tomato_Effective_Richness_English_Name_2, "Dunn_tomato_Effective_Richness_English_Name_2.csv")
Dunn_tomato_Effective_Richness_English_Name_2 <- dunnTest(alpha_500_PD$Effective_Richness, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Effective_Richness_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_English_Name_2$res);CLD_Dunn_tomato_Effective_Richness_English_Name_2 


kruskal.test(Effective_Richness ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 7.5246, df = 1, p-value = 0.006086

kruskal.test(Effective_Richness ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 29.5, df = 4, p-value = 6.186e-06

kruskal.test(Effective_Richness ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Effective_Richness ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Effective_Richness ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Effective_Richness ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283Kruskal-Wallis (no differences were observed)

Dunn_tomato_Effective_Richness_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Effective_Richness ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_Fruit_color
write_csv(Dunn_tomato_Effective_Richness_Fruit_color, "Dunn_tomato_Effective_Richness_Fruit_color.csv")
Dunn_tomato_Effective_Richness_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Effective_Richness, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Effective_Richness_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_Fruit_color$res);CLD_Dunn_tomato_Effective_Richness_Fruit_color 

#insect resistance
Dunn_tomato_Effective_Richness_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Effective_Richness ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_Insect_resistance
write_csv(Dunn_tomato_Effective_Richness_Insect_resistance, "Dunn_tomato_Effective_Richness_Insect_resistance.csv")
Dunn_tomato_Effective_Richness_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Effective_Richness, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Effective_Richness_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_Insect_resistance$res);CLD_Dunn_tomato_Effective_Richness_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Effective_Richness_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Effective_Richness ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_TMV_resistance
write_csv(Dunn_tomato_Effective_Richness_TMV_resistance, "Dunn_tomato_Effective_Richness_TMV_resistance.csv")
Dunn_tomato_Effective_Richness_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Effective_Richness, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Effective_Richness_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_TMV_resistance$res);CLD_Dunn_tomato_Effective_Richness_TMV_resistance 

#seed weoght
Dunn_tomato_Effective_Richness_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Effective_Richness ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_seedweight_1000
write_csv(Dunn_tomato_Effective_Richness_seedweight_1000, "Dunn_tomato_Effective_Richness_seedweight_1000.csv")
Dunn_tomato_Effective_Richness_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Effective_Richness, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Effective_Richness_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_seedweight_1000$res);CLD_Dunn_tomato_Effective_Richness_seedweight_1000 


#seed SHAPE
Dunn_tomato_Effective_Richness_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Effective_Richness ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Effective_Richness_Fruit_shape
write_csv(Dunn_tomato_Effective_Richness_Fruit_shape, "Dunn_tomato_Effective_Richness_Fruit_shape.csv")
Dunn_tomato_Effective_Richness_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Effective_Richness, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Effective_Richness_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Effective_Richness_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Effective_Richness_Fruit_shape$res);CLD_Dunn_tomato_Effective_Richness_Fruit_shape 
###########################################################################################
###########################################################################################
###########################################################################################
#PHYLOGENETIC DIVERSITY
kruskal.test(PD ~ English_Name_2, data = alpha_500_PD)
#1 Kruskal-Wallis chi-squared = 761.57, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(PD ~ Production_site_2, data = alpha_500_PD)
#data:  PD by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_PD_Production_site_2=dunn_test(alpha_500_PD, PD ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_Production_site_2
write_csv(Dunn_tomato_PD_Production_site_2, "Dunn_tomato_PD_Production_site_2.csv")
Dunn_tomato_PD_Production_site_2 <- dunnTest(alpha_500_PD$PD, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_PD_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_Production_site_2$res);CLD_Dunn_tomato_PD_Production_site_2 

#Dunn_tomato_PD_bonferon=dunn_test(alpha_500_PD, PD ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_PD_bonferon
#write_csv(Dunn_tomato_PD_bonferon, "Dunn_tomato_PD_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_PD_English_Name_2=dunn_test(alpha_500_PD, PD ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_English_Name_2
write_csv(Dunn_tomato_PD_English_Name_2, "Dunn_tomato_PD_English_Name_2.csv")
Dunn_tomato_PD_English_Name_2 <- dunnTest(alpha_500_PD$PD, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_PD_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_English_Name_2$res);CLD_Dunn_tomato_PD_English_Name_2 

kruskal.test(PD ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 5.0993, df = 1, p-value = 0.02394

kruskal.test(PD ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(PD ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(PD ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(PD ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(PD ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_PD_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, PD ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_Fruit_color
write_csv(Dunn_tomato_PD_Fruit_color, "Dunn_tomato_PD_Fruit_color.csv")
Dunn_tomato_PD_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$PD, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_PD_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_Fruit_color$res);CLD_Dunn_tomato_PD_Fruit_color 

#insect resistance
Dunn_tomato_PD_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, PD ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_Insect_resistance
write_csv(Dunn_tomato_PD_Insect_resistance, "Dunn_tomato_PD_Insect_resistance.csv")
Dunn_tomato_PD_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$PD, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_PD_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_Insect_resistance$res);CLD_Dunn_tomato_PD_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_PD_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, PD ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_TMV_resistance
write_csv(Dunn_tomato_PD_TMV_resistance, "Dunn_tomato_PD_TMV_resistance.csv")
Dunn_tomato_PD_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$PD, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_PD_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_TMV_resistance$res);CLD_Dunn_tomato_PD_TMV_resistance 

#seed weoght
Dunn_tomato_PD_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, PD ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_seedweight_1000
write_csv(Dunn_tomato_PD_seedweight_1000, "Dunn_tomato_PD_seedweight_1000.csv")
Dunn_tomato_PD_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$PD, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_PD_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_seedweight_1000$res);CLD_Dunn_tomato_PD_seedweight_1000 


#seed SHAPE
Dunn_tomato_PD_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, PD ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_PD_Fruit_shape
write_csv(Dunn_tomato_PD_Fruit_shape, "Dunn_tomato_PD_Fruit_shape.csv")
Dunn_tomato_PD_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$PD, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_PD_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_PD_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_PD_Fruit_shape$res);CLD_Dunn_tomato_PD_Fruit_shape 
###########################################################################################
######################################################################################################################################################################################
#SHANNON INDEX

kruskal.test(Shannon_Index ~ English_Name_2, data = alpha_500_PD)
#Kruskal-Wallis chi-squared = 643.28, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Shannon_Index ~ Production_site_2, data = alpha_500_PD)
#data:  Shannon_Index by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Shannon_Index_Production_site_2=dunn_test(alpha_500_PD, Shannon_Index ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_Production_site_2
write_csv(Dunn_tomato_Shannon_Index_Production_site_2, "Dunn_tomato_Shannon_Index_Production_site_2.csv")
Dunn_tomato_Shannon_Index_Production_site_2 <- dunnTest(alpha_500_PD$Shannon_Index, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Shannon_Index_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_Production_site_2$res);CLD_Dunn_tomato_Shannon_Index_Production_site_2 

#Dunn_tomato_Shannon_Index_bonferon=dunn_test(alpha_500_PD, Shannon_Index ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Shannon_Index_bonferon
#write_csv(Dunn_tomato_Shannon_Index_bonferon, "Dunn_tomato_Shannon_Index_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Shannon_Index_English_Name_2=dunn_test(alpha_500_PD, Shannon_Index ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_English_Name_2
write_csv(Dunn_tomato_Shannon_Index_English_Name_2, "Dunn_tomato_Shannon_Index_English_Name_2.csv")
Dunn_tomato_Shannon_Index_English_Name_2 <- dunnTest(alpha_500_PD$Shannon_Index, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Shannon_Index_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_English_Name_2$res);CLD_Dunn_tomato_Shannon_Index_English_Name_2 

kruskal.test(Shannon_Index ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 12.938, df = 1, p-value = 0.000322

kruskal.test(Shannon_Index ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(Shannon_Index ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Shannon_Index ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Shannon_Index ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Shannon_Index ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Shannon_Index_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Shannon_Index ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_Fruit_color
write_csv(Dunn_tomato_Shannon_Index_Fruit_color, "Dunn_tomato_Shannon_Index_Fruit_color.csv")
Dunn_tomato_Shannon_Index_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Shannon_Index, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Shannon_Index_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_Fruit_color$res);CLD_Dunn_tomato_Shannon_Index_Fruit_color 

#insect resistance
Dunn_tomato_Shannon_Index_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Shannon_Index ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_Insect_resistance
write_csv(Dunn_tomato_Shannon_Index_Insect_resistance, "Dunn_tomato_Shannon_Index_Insect_resistance.csv")
Dunn_tomato_Shannon_Index_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Shannon_Index, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Shannon_Index_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_Insect_resistance$res);CLD_Dunn_tomato_Shannon_Index_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Shannon_Index_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Shannon_Index ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_TMV_resistance
write_csv(Dunn_tomato_Shannon_Index_TMV_resistance, "Dunn_tomato_Shannon_Index_TMV_resistance.csv")
Dunn_tomato_Shannon_Index_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Shannon_Index, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Shannon_Index_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_TMV_resistance$res);CLD_Dunn_tomato_Shannon_Index_TMV_resistance 

#seed weoght
Dunn_tomato_Shannon_Index_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Shannon_Index ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_seedweight_1000
write_csv(Dunn_tomato_Shannon_Index_seedweight_1000, "Dunn_tomato_Shannon_Index_seedweight_1000.csv")
Dunn_tomato_Shannon_Index_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Shannon_Index, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Shannon_Index_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_seedweight_1000$res);CLD_Dunn_tomato_Shannon_Index_seedweight_1000 


#seed SHAPE
Dunn_tomato_Shannon_Index_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Shannon_Index ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Index_Fruit_shape
write_csv(Dunn_tomato_Shannon_Index_Fruit_shape, "Dunn_tomato_Shannon_Index_Fruit_shape.csv")
Dunn_tomato_Shannon_Index_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Shannon_Index, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Shannon_Index_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Index_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Index_Fruit_shape$res);CLD_Dunn_tomato_Shannon_Index_Fruit_shape 
####################################################################################################################################################################################################
###########################################################################################

kruskal.test(Shannon_Effective ~ English_Name_2, data = alpha_500_PD)
#Kruskal-Wallis chi-squared = 643.28, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Shannon_Effective ~ Production_site_2, data = alpha_500_PD)
#data:  Shannon_Effective by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Shannon_Effective_Production_site_2=dunn_test(alpha_500_PD, Shannon_Effective ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_Production_site_2
write_csv(Dunn_tomato_Shannon_Effective_Production_site_2, "Dunn_tomato_Shannon_Effective_Production_site_2.csv")
Dunn_tomato_Shannon_Effective_Production_site_2 <- dunnTest(alpha_500_PD$Shannon_Effective, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Shannon_Effective_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_Production_site_2$res);CLD_Dunn_tomato_Shannon_Effective_Production_site_2 

#Dunn_tomato_Shannon_Effective_bonferon=dunn_test(alpha_500_PD, Shannon_Effective ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Shannon_Effective_bonferon
#write_csv(Dunn_tomato_Shannon_Effective_bonferon, "Dunn_tomato_Shannon_Effective_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Shannon_Effective_English_Name_2=dunn_test(alpha_500_PD, Shannon_Effective ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_English_Name_2
write_csv(Dunn_tomato_Shannon_Effective_English_Name_2, "Dunn_tomato_Shannon_Effective_English_Name_2.csv")
Dunn_tomato_Shannon_Effective_English_Name_2 <- dunnTest(alpha_500_PD$Shannon_Effective, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Shannon_Effective_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_English_Name_2$res);CLD_Dunn_tomato_Shannon_Effective_English_Name_2 

kruskal.test(Shannon_Effective ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 12.944, df = 1, p-value = 0.0003209

kruskal.test(Shannon_Effective ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(Shannon_Effective ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Shannon_Effective ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Shannon_Effective ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Shannon_Effective ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Shannon_Effective_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Shannon_Effective ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_Fruit_color
write_csv(Dunn_tomato_Shannon_Effective_Fruit_color, "Dunn_tomato_Shannon_Effective_Fruit_color.csv")
Dunn_tomato_Shannon_Effective_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Shannon_Effective, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Shannon_Effective_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_Fruit_color$res);CLD_Dunn_tomato_Shannon_Effective_Fruit_color 

#insect resistance
Dunn_tomato_Shannon_Effective_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Shannon_Effective ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_Insect_resistance
write_csv(Dunn_tomato_Shannon_Effective_Insect_resistance, "Dunn_tomato_Shannon_Effective_Insect_resistance.csv")
Dunn_tomato_Shannon_Effective_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Shannon_Effective, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Shannon_Effective_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_Insect_resistance$res);CLD_Dunn_tomato_Shannon_Effective_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Shannon_Effective_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Shannon_Effective ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_TMV_resistance
write_csv(Dunn_tomato_Shannon_Effective_TMV_resistance, "Dunn_tomato_Shannon_Effective_TMV_resistance.csv")
Dunn_tomato_Shannon_Effective_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Shannon_Effective, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Shannon_Effective_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_TMV_resistance$res);CLD_Dunn_tomato_Shannon_Effective_TMV_resistance 

#seed weoght
Dunn_tomato_Shannon_Effective_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Shannon_Effective ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_seedweight_1000
write_csv(Dunn_tomato_Shannon_Effective_seedweight_1000, "Dunn_tomato_Shannon_Effective_seedweight_1000.csv")
Dunn_tomato_Shannon_Effective_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Shannon_Effective, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Shannon_Effective_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_seedweight_1000$res);CLD_Dunn_tomato_Shannon_Effective_seedweight_1000 


#seed SHAPE
Dunn_tomato_Shannon_Effective_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Shannon_Effective ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Shannon_Effective_Fruit_shape
write_csv(Dunn_tomato_Shannon_Effective_Fruit_shape, "Dunn_tomato_Shannon_Effective_Fruit_shape.csv")
Dunn_tomato_Shannon_Effective_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Shannon_Effective, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Shannon_Effective_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Shannon_Effective_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Shannon_Effective_Fruit_shape$res);CLD_Dunn_tomato_Shannon_Effective_Fruit_shape 
############################################################################################################################################
############################################################################################################################################
#EVENNESS
kruskal.test(Evenness ~ English_Name_2, data = alpha_500_PD)
#Kruskal-Wallis chi-squared = 643.28, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Evenness ~ Production_site_2, data = alpha_500_PD)
#data:  Evenness by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Evenness_Production_site_2=dunn_test(alpha_500_PD, Evenness ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_Production_site_2
write_csv(Dunn_tomato_Evenness_Production_site_2, "Dunn_tomato_Evenness_Production_site_2.csv")
Dunn_tomato_Evenness_Production_site_2 <- dunnTest(alpha_500_PD$Evenness, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Evenness_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_Production_site_2$res);CLD_Dunn_tomato_Evenness_Production_site_2 

#Dunn_tomato_Evenness_bonferon=dunn_test(alpha_500_PD, Evenness ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Evenness_bonferon
#write_csv(Dunn_tomato_Evenness_bonferon, "Dunn_tomato_Evenness_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Evenness_English_Name_2=dunn_test(alpha_500_PD, Evenness ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_English_Name_2
write_csv(Dunn_tomato_Evenness_English_Name_2, "Dunn_tomato_Evenness_English_Name_2.csv")
Dunn_tomato_Evenness_English_Name_2 <- dunnTest(alpha_500_PD$Evenness, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Evenness_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_English_Name_2$res);CLD_Dunn_tomato_Evenness_English_Name_2 

kruskal.test(Evenness ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 12.026, df = 1, p-value = 0.0005246

kruskal.test(Evenness ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(Evenness ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Evenness ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Evenness ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Evenness ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Evenness_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Evenness ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_Fruit_color
write_csv(Dunn_tomato_Evenness_Fruit_color, "Dunn_tomato_Evenness_Fruit_color.csv")
Dunn_tomato_Evenness_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Evenness, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Evenness_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_Fruit_color$res);CLD_Dunn_tomato_Evenness_Fruit_color 

#insect resistance
Dunn_tomato_Evenness_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Evenness ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_Insect_resistance
write_csv(Dunn_tomato_Evenness_Insect_resistance, "Dunn_tomato_Evenness_Insect_resistance.csv")
Dunn_tomato_Evenness_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Evenness, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Evenness_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_Insect_resistance$res);CLD_Dunn_tomato_Evenness_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Evenness_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Evenness ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_TMV_resistance
write_csv(Dunn_tomato_Evenness_TMV_resistance, "Dunn_tomato_Evenness_TMV_resistance.csv")
Dunn_tomato_Evenness_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Evenness, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Evenness_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_TMV_resistance$res);CLD_Dunn_tomato_Evenness_TMV_resistance 

#seed weoght
Dunn_tomato_Evenness_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Evenness ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_seedweight_1000
write_csv(Dunn_tomato_Evenness_seedweight_1000, "Dunn_tomato_Evenness_seedweight_1000.csv")
Dunn_tomato_Evenness_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Evenness, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Evenness_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_seedweight_1000$res);CLD_Dunn_tomato_Evenness_seedweight_1000 


#seed SHAPE
Dunn_tomato_Evenness_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Evenness ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Evenness_Fruit_shape
write_csv(Dunn_tomato_Evenness_Fruit_shape, "Dunn_tomato_Evenness_Fruit_shape.csv")
Dunn_tomato_Evenness_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Evenness, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Evenness_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Evenness_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Evenness_Fruit_shape$res);CLD_Dunn_tomato_Evenness_Fruit_shape 

##############################################################################################################################################################################################################################################
#Simpson_Index

kruskal.test(Simpson_Index ~ English_Name_2, data = alpha_500_PD)
#Kruskal-Wallis chi-squared = 643.28, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Simpson_Index ~ Production_site_2, data = alpha_500_PD)
#data:  Simpson_Index by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Simpson_Index_Production_site_2=dunn_test(alpha_500_PD, Simpson_Index ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_Production_site_2
write_csv(Dunn_tomato_Simpson_Index_Production_site_2, "Dunn_tomato_Simpson_Index_Production_site_2.csv")
Dunn_tomato_Simpson_Index_Production_site_2 <- dunnTest(alpha_500_PD$Simpson_Index, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Simpson_Index_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_Production_site_2$res);CLD_Dunn_tomato_Simpson_Index_Production_site_2 

#Dunn_tomato_Simpson_Index_bonferon=dunn_test(alpha_500_PD, Simpson_Index ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Simpson_Index_bonferon
#write_csv(Dunn_tomato_Simpson_Index_bonferon, "Dunn_tomato_Simpson_Index_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Simpson_Index_English_Name_2=dunn_test(alpha_500_PD, Simpson_Index ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_English_Name_2
write_csv(Dunn_tomato_Simpson_Index_English_Name_2, "Dunn_tomato_Simpson_Index_English_Name_2.csv")
Dunn_tomato_Simpson_Index_English_Name_2 <- dunnTest(alpha_500_PD$Simpson_Index, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Simpson_Index_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_English_Name_2$res);CLD_Dunn_tomato_Simpson_Index_English_Name_2 

kruskal.test(Simpson_Index ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 11.059, df = 1, p-value = 0.0008826

kruskal.test(Simpson_Index ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(Simpson_Index ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Simpson_Index ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Simpson_Index ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Simpson_Index ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Simpson_Index_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Simpson_Index ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_Fruit_color
write_csv(Dunn_tomato_Simpson_Index_Fruit_color, "Dunn_tomato_Simpson_Index_Fruit_color.csv")
Dunn_tomato_Simpson_Index_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Simpson_Index, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Simpson_Index_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_Fruit_color$res);CLD_Dunn_tomato_Simpson_Index_Fruit_color 

#insect resistance
Dunn_tomato_Simpson_Index_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Simpson_Index ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_Insect_resistance
write_csv(Dunn_tomato_Simpson_Index_Insect_resistance, "Dunn_tomato_Simpson_Index_Insect_resistance.csv")
Dunn_tomato_Simpson_Index_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Simpson_Index, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Simpson_Index_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_Insect_resistance$res);CLD_Dunn_tomato_Simpson_Index_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Simpson_Index_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Simpson_Index ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_TMV_resistance
write_csv(Dunn_tomato_Simpson_Index_TMV_resistance, "Dunn_tomato_Simpson_Index_TMV_resistance.csv")
Dunn_tomato_Simpson_Index_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Simpson_Index, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Simpson_Index_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_TMV_resistance$res);CLD_Dunn_tomato_Simpson_Index_TMV_resistance 

#seed weoght
Dunn_tomato_Simpson_Index_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Simpson_Index ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_seedweight_1000
write_csv(Dunn_tomato_Simpson_Index_seedweight_1000, "Dunn_tomato_Simpson_Index_seedweight_1000.csv")
Dunn_tomato_Simpson_Index_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Simpson_Index, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Simpson_Index_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_seedweight_1000$res);CLD_Dunn_tomato_Simpson_Index_seedweight_1000 


#seed SHAPE
Dunn_tomato_Simpson_Index_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Simpson_Index ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Index_Fruit_shape
write_csv(Dunn_tomato_Simpson_Index_Fruit_shape, "Dunn_tomato_Simpson_Index_Fruit_shape.csv")
Dunn_tomato_Simpson_Index_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Simpson_Index, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Simpson_Index_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Index_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Index_Fruit_shape$res);CLD_Dunn_tomato_Simpson_Index_Fruit_shape 

##############################################################################################################################################################################################################################################
#Effective simpson index

kruskal.test(Simpson_Effective ~English_Name_2, data = alpha_500_PD)
#Kruskal-Wallis chi-squared = 643.28, df = 99, p-value < 2.2e-16
#There was overall significant differences between the different cultivars

kruskal.test(Simpson_Effective ~ Production_site_2, data = alpha_500_PD)
#data:  Simpson_Effective by Production_site_2
#Kruskal-Wallis chi-squared = 177.59, df = 11, p-value < 2.2e-16

# Dunn's Test 
Dunn_tomato_Simpson_Effective_Production_site_2=dunn_test(alpha_500_PD, Simpson_Effective ~ Production_site_2, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_Production_site_2
write_csv(Dunn_tomato_Simpson_Effective_Production_site_2, "Dunn_tomato_Simpson_Effective_Production_site_2.csv")
Dunn_tomato_Simpson_Effective_Production_site_2 <- dunnTest(alpha_500_PD$Simpson_Effective, alpha_500_PD$Production_site_2, method = "bonferroni");Dunn_tomato_Simpson_Effective_Production_site_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_Production_site_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_Production_site_2$res);CLD_Dunn_tomato_Simpson_Effective_Production_site_2 

#Dunn_tomato_Simpson_Effective_bonferon=dunn_test(alpha_500_PD, Simpson_Effective ~ Production_site_2, p.adjust.method = "BH", detailed = FALSE);Dunn_tomato_Simpson_Effective_bonferon
#write_csv(Dunn_tomato_Simpson_Effective_bonferon, "Dunn_tomato_Simpson_Effective_bonferon.csv")
#*****The use of BH is rather less stringent; however an alternative to bonferroni is holm
Dunn_tomato_Simpson_Effective_English_Name_2=dunn_test(alpha_500_PD, Simpson_Effective ~ English_Name, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_English_Name_2
write_csv(Dunn_tomato_Simpson_Effective_English_Name_2, "Dunn_tomato_Simpson_Effective_English_Name_2.csv")
Dunn_tomato_Simpson_Effective_English_Name_2 <- dunnTest(alpha_500_PD$Simpson_Effective, alpha_500_PD$English_Name, method = "bonferroni");Dunn_tomato_Simpson_Effective_English_Name_2
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_English_Name_2 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_English_Name_2$res);CLD_Dunn_tomato_Simpson_Effective_English_Name_2 

kruskal.test(Simpson_Effective ~ Fruit_taste, data = alpha_500_PD_Fruit_taste)
#Kruskal-Wallis chi-squared = 11.062, df = 1, p-value = 0.000881

kruskal.test(Simpson_Effective ~ Fruit_color, data = alpha_500_PD_Fruit_color)
#Kruskal-Wallis chi-squared = 35.11, df = 4, p-value = 4.41e-07

kruskal.test(Simpson_Effective ~ Fruit_shape, data = alpha_500_PD_Fruit_shape)
#Kruskal-Wallis chi-squared = 26.851, df = 2, p-value = 1.477e-06

kruskal.test(Simpson_Effective ~ TMV_resistance, data = alpha_500_PD_TMV_resistance)
#Kruskal-Wallis chi-squared = 11.961, df = 2, p-value = 0.002527

kruskal.test(Simpson_Effective ~ Insect_resistance, data = alpha_500_PD_Insect_resistance)
#Kruskal-Wallis chi-squared = 46.262, df = 5, p-value = 8.031e-09 Kruskal-Wallis
kruskal.test(Simpson_Effective ~ seedweight_1000, data = alpha_500_PD_seedweight_1000)
#Kruskal-Wallis chi-squared = 12.03, df = 3, p-value = 0.007283 Kruskal-Wallis (no differences were observed)

Dunn_tomato_Simpson_Effective_Fruit_color=dunn_test(alpha_500_PD_Fruit_color, Simpson_Effective ~ Fruit_color, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_Fruit_color
write_csv(Dunn_tomato_Simpson_Effective_Fruit_color, "Dunn_tomato_Simpson_Effective_Fruit_color.csv")
Dunn_tomato_Simpson_Effective_Fruit_color <- dunnTest(alpha_500_PD_Fruit_color$Simpson_Effective, alpha_500_PD_Fruit_color$Fruit_color, method = "bonferroni");Dunn_tomato_Simpson_Effective_Fruit_color
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_Fruit_color = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_Fruit_color$res);CLD_Dunn_tomato_Simpson_Effective_Fruit_color 

#insect resistance
Dunn_tomato_Simpson_Effective_Insect_resistance=dunn_test(alpha_500_PD_Insect_resistance, Simpson_Effective ~ Insect_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_Insect_resistance
write_csv(Dunn_tomato_Simpson_Effective_Insect_resistance, "Dunn_tomato_Simpson_Effective_Insect_resistance.csv")
Dunn_tomato_Simpson_Effective_Insect_resistance <- dunnTest(alpha_500_PD_Insect_resistance$Simpson_Effective, alpha_500_PD_Insect_resistance$Insect_resistance, method = "bonferroni");Dunn_tomato_Simpson_Effective_Insect_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_Insect_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_Insect_resistance$res);CLD_Dunn_tomato_Simpson_Effective_Insect_resistance 

#TMV RESISTANCE
Dunn_tomato_Simpson_Effective_TMV_resistance=dunn_test(alpha_500_PD_TMV_resistance, Simpson_Effective ~ TMV_resistance, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_TMV_resistance
write_csv(Dunn_tomato_Simpson_Effective_TMV_resistance, "Dunn_tomato_Simpson_Effective_TMV_resistance.csv")
Dunn_tomato_Simpson_Effective_TMV_resistance <- dunnTest(alpha_500_PD_TMV_resistance$Simpson_Effective, alpha_500_PD_TMV_resistance$TMV_resistance, method = "bonferroni");Dunn_tomato_Simpson_Effective_TMV_resistance
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_TMV_resistance = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_TMV_resistance$res);CLD_Dunn_tomato_Simpson_Effective_TMV_resistance 

#seed weoght
Dunn_tomato_Simpson_Effective_seedweight_1000=dunn_test(alpha_500_PD_seedweight_1000, Simpson_Effective ~ seedweight_1000, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_seedweight_1000
write_csv(Dunn_tomato_Simpson_Effective_seedweight_1000, "Dunn_tomato_Simpson_Effective_seedweight_1000.csv")
Dunn_tomato_Simpson_Effective_seedweight_1000 <- dunnTest(alpha_500_PD_seedweight_1000$Simpson_Effective, alpha_500_PD_seedweight_1000$seedweight_1000, method = "bonferroni");Dunn_tomato_Simpson_Effective_seedweight_1000
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_seedweight_1000 = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_seedweight_1000$res);CLD_Dunn_tomato_Simpson_Effective_seedweight_1000 


#seed SHAPE
Dunn_tomato_Simpson_Effective_Fruit_shape=dunn_test(alpha_500_PD_Fruit_shape, Simpson_Effective ~ Fruit_shape, p.adjust.method = "bonferroni", detailed = FALSE);Dunn_tomato_Simpson_Effective_Fruit_shape
write_csv(Dunn_tomato_Simpson_Effective_Fruit_shape, "Dunn_tomato_Simpson_Effective_Fruit_shape.csv")
Dunn_tomato_Simpson_Effective_Fruit_shape <- dunnTest(alpha_500_PD_Fruit_shape$Simpson_Effective, alpha_500_PD_Fruit_shape$Fruit_shape, method = "bonferroni");Dunn_tomato_Simpson_Effective_Fruit_shape
# Obtaining significance letters
library(rcompanion)
CLD_Dunn_tomato_Simpson_Effective_Fruit_shape = cldList(P.adj ~ Comparison, data=Dunn_tomato_Simpson_Effective_Fruit_shape$res);CLD_Dunn_tomato_Simpson_Effective_Fruit_shape 

##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
#GGPLOTS FOR ALL DIVERSITY INDICIES 
#1:richness_shannon_pd_Simpson_eveness.tiff
plot_richness=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Richness", fill= "English_Name_2",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=18),
        axis.title.y = element_text(face="bold", colour="gray33", size=18),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Production_site_2);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#SHANON DIVERSITY
plot_Shannon=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Shannon_Index", fill= "English_Name_2", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=18),
        axis.title.y = element_text(face="bold", colour="gray33", size=18),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


#Phylogenetic diversity
plot_PD=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "PD", fill= "English_Name_2", ylab = "PD", xlab = "Tomato cultivars")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=18),
        axis.title.y = element_text(face="bold", colour="gray33", size=18),
        axis.text.x = element_text( colour="gray33", size=16),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD
#ggsave("plot_PD.TIFF", plot = plot_PD)

###########################################################################################################################################################
#tiff("richness_shannon_pd_Simpson_eveness.tiff", units="in", width=20, height=18, res=300)
tiff("richness_shannon_pd_Simpson_eveness.tiff", units="in", width=20, height=18, res=300)
richness_shannon_pd_Simpson_eveness.tiff=(p2/plot_richness/plot_Shannon/plot_PD);richness_shannon_pd_Simpson_eveness.tiff
ggsave("richness_shannon_pd_Simpson_eveness.tiff", plot = richness_shannon_pd_Simpson_eveness.tiff)
dev.off()
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
#ggplots 
plot_Simpson=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Simpson_Index", fill= "English_Name_2", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=18),
        axis.title.y = element_text(face="bold", colour="gray33", size=18),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

#EVENNESS
#ggplots 
plot_evenness=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Evenness", fill= "English_Name_2", ylab = "Evenness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=18),
        axis.title.y = element_text(face="bold", colour="gray33", size=18),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

#Effective_alpha_diversity
tiff("Effective_alpha_diversity.tiff", units="in", width=20, height=15, res=300)
plot_effective_richness=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Effective_Richness", fill= "English_Name_2",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=16),
        axis.title.y = element_text(face="bold", colour="gray33", size=16),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)


plot_effective_Shannon=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Shannon_Effective", fill= "English_Name_2",  ylab = "Eff. Shannon diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('gray', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=16),
        axis.title.y = element_text(face="bold", colour="gray33", size=16),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)


plot_effective_Simpson=ggboxplot(alpha_500_PD, x = "English_Name_2", y = "Simpson_Effective", fill= "English_Name_2",  ylab = "Eff. Simpson diversity", xlab = "Cultivar")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('gray', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=16),
        axis.title.y = element_text(face="bold", colour="gray33", size=16),
        axis.text.x = element_text( colour="gray33", size=16),
        axis.text.y  = element_text( colour="gray33", size=16))+
  theme(text=element_text(size=16,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)


tiff("pd_Simpson_eveness_Effective_alpha_diversity.tiff", units="in", width=21, height=20, res=300)
pd_Simpson_eveness_Effective_alpha_diversity.tiff=(plot_Simpson/plot_evenness/plot_effective_richness/plot_effective_Shannon/plot_effective_Simpson);pd_Simpson_eveness_Effective_alpha_diversity.tiff
ggsave("pd_Simpson_eveness_Effective_alpha_diversity.tiff", plot = pd_Simpson_eveness_Effective_alpha_diversity.tiff)
dev.off()
###########################################################################################################################################################
#########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
###########################################################################################################################################################
#########################################################################################################################################################
###########################################################################################################################################################
#Production site
plot_richness_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Richness", fill= "Production_site_2",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_prod_site
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Production_site_2);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Effective_Richness", fill= "Production_site_2",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_prod_site
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Shannon_Index", fill= "Production_site_2", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_prod_site
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Shannon_Effective", fill= "Production_site_2",  ylab = "Eff. Shannon diversity", xlab = "Region of seed of production")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('gray', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_text( colour="gray33", size=17),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_prod_site
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "PD", fill= "Production_site_2", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_prod_site
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Simpson_Index", fill= "Production_site_2", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_prod_site
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Simpson_Effective", fill= "Production_site_2",  ylab = "Effe. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('gray', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_prod_site
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_prod_site=ggboxplot(alpha_500_PD, x = "Production_site_2", y = "Evenness", fill= "Production_site_2", ylab = "Evenness", xlab = "Region of seed of production")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_color_gradientn(colours = c('grey', 'pink'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=17),
        axis.title.y = element_text(face="bold", colour="gray33", size=17),
        axis.text.x = element_text( colour="gray33", size=17),
        axis.text.y  = element_text( colour="gray33", size=17))+
  theme(text=element_text(size=17,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_prod_site
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_priduction_site.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_priduction_site.tiff=(plot_richness_prod_site+plot_effective_richness_prod_site+plot_Shannon_prod_site+plot_effective_Shannon_prod_site+plot_PD_prod_site+plot_Simpson_prod_site+plot_effective_Simpson_prod_site+plot_evenness_prod_site+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_priduction_site.tiff
ggsave("richness_shannon_pd_Simpson_eveness_priduction_site.tiff", plot = richness_shannon_pd_Simpson_eveness_priduction_site.tiff)
dev.off()
######################################################################################################################################################
######################################################################################################################################################
#Dataframe
alpha_500_PD_Fruit_color=alpha_500_PD[!is.na(alpha_500_PD$Fruit_color),];alpha_500_PD_Fruit_color
#####################################################################################################################################################
myco_col=c("grey35","#00AFBB", "yellow3","chartreuse4","#FB9A99");myco_col
#FRUIT COLOR FRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLORFRUIT COLOR

#alpha_500_PD_Fruit_color
plot_richness_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Richness", fill= "Fruit_color",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Fruit_color
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Fruit_color);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Effective_Richness", fill= "Fruit_color",  ylab = "Eff.Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Fruit_color
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Shannon_Index", fill= "Fruit_color", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Fruit_color
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Shannon_Effective", fill= "Fruit_color",  ylab = "Eff. Shannon diversity", xlab = "Berry Color")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Fruit_color
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "PD", fill= "Fruit_color", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Fruit_color
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Simpson_Index", fill= "Fruit_color", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Fruit_color
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Simpson_Effective", fill= "Fruit_color",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Fruit_color
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Fruit_color=ggboxplot(alpha_500_PD_Fruit_color, x = "Fruit_color", y = "Evenness", fill= "Fruit_color", ylab = "Evenness", xlab = "Berry Color")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_col)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Fruit_color
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_Fruit_color.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_Fruit_color.tiff=(plot_richness_Fruit_color+plot_effective_richness_Fruit_color+plot_Shannon_Fruit_color+plot_effective_Shannon_Fruit_color+plot_PD_Fruit_color+plot_Simpson_Fruit_color+plot_effective_Simpson_Fruit_color+plot_evenness_Fruit_color+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Fruit_color.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Fruit_color.tiff", plot = richness_shannon_pd_Simpson_eveness_Fruit_color.tiff)
dev.off()


##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
alpha_500_PD_Fruit_taste=alpha_500_PD[!is.na(alpha_500_PD$Fruit_taste),];alpha_500_PD_Fruit_taste
myco_taste=c("#FB9A99","chartreuse4");myco_taste
##############################################################################################################################################################################################################################################
plot_richness_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Richness", fill= "Fruit_taste",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Fruit_taste
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Fruit_taste);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Effective_Richness", fill= "Fruit_taste",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Fruit_taste
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Shannon_Index", fill= "Fruit_taste", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Fruit_taste
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Shannon_Effective", fill= "Fruit_taste",  ylab = "Eff. Shannon diversity", xlab = "Berry taste")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Fruit_taste
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "PD", fill= "Fruit_taste", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Fruit_taste
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Simpson_Index", fill= "Fruit_taste", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Fruit_taste
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Simpson_Effective", fill= "Fruit_taste",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Fruit_taste
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Fruit_taste=ggboxplot(alpha_500_PD_Fruit_taste, x = "Fruit_taste", y = "Evenness", fill= "Fruit_taste", ylab = "Evenness", xlab = "Berry taste")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_taste)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Fruit_taste
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_Fruit_taste.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_Fruit_taste.tiff=(plot_richness_Fruit_taste+plot_effective_richness_Fruit_taste+plot_Shannon_Fruit_taste+plot_effective_Shannon_Fruit_taste+plot_PD_Fruit_taste+plot_Simpson_Fruit_taste+plot_effective_Simpson_Fruit_taste+plot_evenness_Fruit_taste+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Fruit_taste.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Fruit_taste.tiff", plot = richness_shannon_pd_Simpson_eveness_Fruit_taste.tiff)
dev.off()
##############################################################################################################################################################################################################################################

alpha_500_PD_Fruit_shape=alpha_500_PD[!is.na(alpha_500_PD$Fruit_shape),];alpha_500_PD_Fruit_shape
myco_shp=c("darkturquoise","chartreuse4","#FB9A99");myco_shp
##############################################################################################################################################################################################################################################
alpha_500_PD_Fruit_shape=alpha_500_PD[!is.na(alpha_500_PD$Fruit_shape),];alpha_500_PD_Fruit_shape
tiff("richness_shannon_pd_Simpson_eveness_Fruit_shape.tiff", units="in", width=20, height=20, res=300)
plot_richness_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Richness", fill= "Fruit_shape",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Fruit_shape
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Fruit_shape);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Effective_Richness", fill= "Fruit_shape",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Fruit_shape
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Shannon_Index", fill= "Fruit_shape", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Fruit_shape
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Shannon_Effective", fill= "Fruit_shape",  ylab = "Eff. Shannon diversity", xlab = "Berry shape")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Fruit_shape
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "PD", fill= "Fruit_shape", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Fruit_shape
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Simpson_Index", fill= "Fruit_shape", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Fruit_shape
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Simpson_Effective", fill= "Fruit_shape",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Fruit_shape
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Fruit_shape=ggboxplot(alpha_500_PD_Fruit_shape, x = "Fruit_shape", y = "Evenness", fill= "Fruit_shape", ylab = "Evenness", xlab = "Berry shape")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_shp)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Fruit_shape
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

richness_shannon_pd_Simpson_eveness_Fruit_shape.tiff=(plot_richness_Fruit_shape+plot_effective_richness_Fruit_shape+plot_Shannon_Fruit_shape+plot_effective_Shannon_Fruit_shape+plot_PD_Fruit_shape+plot_Simpson_Fruit_shape+plot_effective_Simpson_Fruit_shape+plot_evenness_Fruit_shape+ plot_layout(ncol= 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Fruit_shape.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Fruit_shape.tiff", plot = richness_shannon_pd_Simpson_eveness_Fruit_shape.tiff)
dev.off()
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
alpha_500_PD_TMV_resistance=alpha_500_PD[!is.na(alpha_500_PD$TMV_resistance),];alpha_500_PD_TMV_resistance
myco_TMV=c("#0B775E","#FB9A99","grey35");myco_TMV
##############################################################################################################################################################################################################################################
plot_richness_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Richness", fill= "TMV_resistance",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_TMV_resistance
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~TMV_resistance);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Effective_Richness", fill= "TMV_resistance",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_TMV_resistance
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Shannon_Index", fill= "TMV_resistance", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_TMV_resistance
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Shannon_Effective", fill= "TMV_resistance",  ylab = "Eff. Shannon diversity", xlab = "TMV resistance")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_TMV_resistance
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "PD", fill= "TMV_resistance", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_TMV_resistance
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Simpson_Index", fill= "TMV_resistance", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_TMV_resistance
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Simpson_Effective", fill= "TMV_resistance",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x =element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_TMV_resistance
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_TMV_resistance=ggboxplot(alpha_500_PD_TMV_resistance, x = "TMV_resistance", y = "Evenness", fill= "TMV_resistance", ylab = "Evenness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_TMV)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_TMV_resistance
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_TMV_resistance.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_TMV_resistance.tiff=(plot_richness_TMV_resistance+plot_effective_richness_TMV_resistance+plot_Shannon_TMV_resistance+plot_effective_Shannon_TMV_resistance+plot_PD_TMV_resistance+plot_Simpson_TMV_resistance+plot_effective_Simpson_TMV_resistance+plot_evenness_TMV_resistance+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_TMV_resistance.tiff
ggsave("richness_shannon_pd_Simpson_eveness_TMV_resistance.tiff", plot = richness_shannon_pd_Simpson_eveness_TMV_resistance.tiff)
dev.off()
##############################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################
alpha_500_PD_Insect_resistance=alpha_500_PD[!is.na(alpha_500_PD$Insect_resistance),];alpha_500_PD_Insect_resistance
myco_insRes=c("yellow3","grey35","chartreuse4","#0B775E","#00AFBB", "#FB9A99");myco_insRes
##############################################################################################################################################################################################################################################
plot_richness_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Richness", fill= "Insect_resistance_2",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Insect_resistance
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Insect_resistance);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Effective_Richness", fill= "Insect_resistance_2",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Insect_resistance
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Shannon_Index", fill= "Insect_resistance_2", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Insect_resistance
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Shannon_Effective", fill= "Insect_resistance_2",  ylab = "EfF. Shannon diversity", xlab = "Insect resistance")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Insect_resistance
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "PD", fill= "Insect_resistance_2", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Insect_resistance
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Simpson_Index", fill= "Insect_resistance_2", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Insect_resistance
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Simpson_Effective", fill= "Insect_resistance_2",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Insect_resistance
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Insect_resistance=ggboxplot(alpha_500_PD_Insect_resistance, x = "Insect_resistance_2", y = "Evenness", fill= "Insect_resistance_2", ylab = "Evenness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_insRes)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Insect_resistance
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_Insect_resistance.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_Insect_resistance.tiff=(plot_richness_Insect_resistance+plot_effective_richness_Insect_resistance+plot_Shannon_Insect_resistance+plot_effective_Shannon_Insect_resistance+plot_PD_Insect_resistance+plot_Simpson_Insect_resistance+plot_effective_Simpson_Insect_resistance+plot_evenness_Insect_resistance+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Insect_resistance.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Insect_resistance.tiff", plot = richness_shannon_pd_Simpson_eveness_Insect_resistance.tiff)
dev.off()

##############################################################################################################################################################################################################################################
alpha_500_PD_seedweight_1000=alpha_500_PD[!is.na(alpha_500_PD$seedweight_1000),];alpha_500_PD_seedweight_1000
myco_weight=c("#0B775E","#FB9A99","yellow3","grey35");myco_weight

##############################################################################################################################################################################################################################################
plot_richness_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Richness", fill= "seedweight_1000",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_seedweight_1000
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~seedweight_1000);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Effective_Richness", fill= "seedweight_1000",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_seedweight_1000
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Shannon_Index", fill= "seedweight_1000", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_seedweight_1000
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Shannon_Effective", fill= "seedweight_1000",  ylab = "Eff. Shannon diversity", xlab = "Weight of 1000 seeds (grams)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_seedweight_1000
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "PD", fill= "seedweight_1000", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_seedweight_1000
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Simpson_Index", fill= "seedweight_1000", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_seedweight_1000
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Simpson_Effective", fill= "seedweight_1000",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_seedweight_1000
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_seedweight_1000=ggboxplot(alpha_500_PD_seedweight_1000, x = "seedweight_1000", y = "Evenness", fill= "seedweight_1000", ylab = "Evenness", xlab = "Weight of 1000 seeds (grams)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_weight)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_seedweight_1000
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_seedweight_1000.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_seedweight_1000.tiff=(plot_richness_seedweight_1000+plot_effective_richness_seedweight_1000+plot_Shannon_seedweight_1000+plot_effective_Shannon_seedweight_1000+plot_PD_seedweight_1000+plot_Simpson_seedweight_1000+plot_effective_Simpson_seedweight_1000+plot_evenness_seedweight_1000+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_seedweight_1000.tiff
ggsave("richness_shannon_pd_Simpson_eveness_seedweight_1000.tiff", plot = richness_shannon_pd_Simpson_eveness_seedweight_1000.tiff)
dev.off()
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
alpha_500_PD_yield=alpha_500_PD[!is.na(alpha_500_PD$yield),];alpha_500_PD_yield
myco_yield=c("chartreuse4","yellow3","#0B775E","#00AFBB","grey35");myco_yield

##############################################################################################################################################################################################################################################
plot_richness_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Richness", fill= "yield",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_yield
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~yield);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Effective_Richness", fill= "yield",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_yield
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Shannon_Index", fill= "yield", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_yield
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Shannon_Effective", fill= "yield",  ylab = "Eff. Shannon diversity", xlab = "Yield (Kg/Sq.Km)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_yield
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "PD", fill= "yield", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_yield
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Simpson_Index", fill= "yield", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_yield
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Simpson_Effective", fill= "yield",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_yield
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_yield=ggboxplot(alpha_500_PD_yield, x = "yield", y = "Evenness", fill= "yield", ylab = "Evenness", xlab = "Yield (Kg/Sq.Km)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_yield)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_yield
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_yield.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_yield.tiff=(plot_richness_yield+plot_effective_richness_yield+plot_Shannon_yield+plot_effective_Shannon_yield+plot_PD_yield+plot_Simpson_yield+plot_effective_Simpson_yield+plot_evenness_yield+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_yield.tiff
ggsave("richness_shannon_pd_Simpson_eveness_yield.tiff", plot = richness_shannon_pd_Simpson_eveness_yield.tiff)
dev.off()

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
alpha_500_PD_ovaries=alpha_500_PD[!is.na(alpha_500_PD$ovaries),];alpha_500_PD_ovaries
myco_ovaries=c("#00AFBB","#FB9A99","chartreuse4");myco_ovaries
#########################################################################################################################
plot_richness_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Richness", fill= "ovaries",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_ovaries
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~ovaries);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Effective_Richness", fill= "ovaries",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_ovaries
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Shannon_Index", fill= "ovaries", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_ovaries
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Shannon_Effective", fill= "ovaries",  ylab = "Eff. Shannon diversity", xlab = "Number of ovaries")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_ovaries
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "PD", fill= "ovaries", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_ovaries
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Simpson_Index", fill= "ovaries", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_ovaries
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Simpson_Effective", fill= "ovaries",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_ovaries
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_ovaries=ggboxplot(alpha_500_PD_ovaries, x = "ovaries", y = "Evenness", fill= "ovaries", ylab = "Evenness", xlab = "Number of ovaries")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_ovaries)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_ovaries
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_ovaries.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_ovaries.tiff=(plot_richness_ovaries+plot_effective_richness_ovaries+plot_Shannon_ovaries+plot_effective_Shannon_ovaries+plot_PD_ovaries+plot_Simpson_ovaries+plot_effective_Simpson_ovaries+plot_evenness_ovaries+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_ovaries.tiff
ggsave("richness_shannon_pd_Simpson_eveness_ovaries.tiff", plot = richness_shannon_pd_Simpson_eveness_ovaries.tiff)
dev.off()
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
alpha_500_PD_Fruit_hori_diameter=alpha_500_PD[!is.na(alpha_500_PD$Fruit_hori_diameter),];alpha_500_PD_Fruit_hori_diameter
myco_horiz=c("chartreuse4","#00AFBB","#FB9A99");myco_horiz
#########################################################################################################################
plot_richness_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Richness", fill= "Fruit_hori_diameter",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Fruit_hori_diameter
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Fruit_hori_diameter);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Effective_Richness", fill= "Fruit_hori_diameter",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Fruit_hori_diameter
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Shannon_Index", fill= "Fruit_hori_diameter", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Fruit_hori_diameter
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Shannon_Effective", fill= "Fruit_hori_diameter",  ylab = "Eff. Shannon diversity", xlab = "Berry horizontal diameter (cm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Fruit_hori_diameter
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "PD", fill= "Fruit_hori_diameter", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Fruit_hori_diameter
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Simpson_Index", fill= "Fruit_hori_diameter", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Fruit_hori_diameter
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Simpson_Effective", fill= "Fruit_hori_diameter",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Fruit_hori_diameter
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Fruit_hori_diameter=ggboxplot(alpha_500_PD_Fruit_hori_diameter, x = "Fruit_hori_diameter", y = "Evenness", fill= "Fruit_hori_diameter", ylab = "Evenness", xlab = "Berry horizontal diameter (cm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_horiz)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Fruit_hori_diameter
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_Fruit_hori_diameter.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_Fruit_hori_diameter.tiff=(plot_richness_Fruit_hori_diameter+plot_effective_richness_Fruit_hori_diameter+plot_Shannon_Fruit_hori_diameter+plot_effective_Shannon_Fruit_hori_diameter+plot_PD_Fruit_hori_diameter+plot_Simpson_Fruit_hori_diameter+plot_effective_Simpson_Fruit_hori_diameter+plot_evenness_Fruit_hori_diameter+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Fruit_hori_diameter.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Fruit_hori_diameter.tiff", plot = richness_shannon_pd_Simpson_eveness_Fruit_hori_diameter.tiff)
dev.off()
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
##########################################################################################################################
alpha_500_PD_Fruit_long_diameter=alpha_500_PD[!is.na(alpha_500_PD$Fruit_long_diameter),];alpha_500_PD_Fruit_long_diameter
myco_longdiam=c("chartreuse4","#00AFBB","#FB9A99" );myco_longdiam
#########################################################################################################################
plot_richness_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Richness", fill= "Fruit_long_diameter",  ylab = "Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") +geom_jitter(width = 0.3)+
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "A")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_richness_Fruit_long_diameter
#ggsave("plot_richness.TIFF", plot = plot_richness)
#+facet_wrap(~Fruit_long_diameter);
#Chao1_simpson.tiff=(plot6 /plot7 );Chao1_simpson.tiff
#ggsave("Chao1_simpson.tiff", plot = Chao1_simpson.tiff)

#Effective_alpha_diversity
plot_effective_richness_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Effective_Richness", fill= "Fruit_long_diameter",  ylab = "Eff. Richness", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "B")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_richness_Fruit_long_diameter
#ggsave("plot_effective_richness.TIFF", plot = plot_effective_richness)

#SHANON DIVERSITY
plot_Shannon_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Shannon_Index", fill= "Fruit_long_diameter", ylab = "Shannon Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "C")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Shannon_Fruit_long_diameter
#ggsave("plot_Shannon.TIFF", plot = plot_Shannon)


plot_effective_Shannon_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Shannon_Effective", fill= "Fruit_long_diameter",  ylab = "Eff. Shannon diversity", xlab = "Berry longitudinal diameter (cm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "D")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Shannon_Fruit_long_diameter
#ggsave("plot_effective_Shannon.TIFF", plot = plot_effective_Shannon)

#Phylogenetic diversity
plot_PD_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "PD", fill= "Fruit_long_diameter", ylab = "PD", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "E")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_PD_Fruit_long_diameter
#ggsave("plot_PD.TIFF", plot = plot_PD)


plot_Simpson_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Simpson_Index", fill= "Fruit_long_diameter", ylab = "Simpson Index", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "F")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_Simpson_Fruit_long_diameter
#ggsave("plot_Simpson.TIFF", plot = plot_Simpson)

plot_effective_Simpson_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Simpson_Effective", fill= "Fruit_long_diameter",  ylab = "Eff. Simpson diversity", xlab = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_blank(),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "G")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_effective_Simpson_Fruit_long_diameter
#ggsave("plot_effective_Simpson.TIFF", plot = plot_effective_Simpson)

#EVENNESS
#ggplots 
plot_evenness_Fruit_long_diameter=ggboxplot(alpha_500_PD_Fruit_long_diameter, x = "Fruit_long_diameter", y = "Evenness", fill= "Fruit_long_diameter", ylab = "Evenness", xlab = "Berry longitudinal diameter (cm)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ theme(legend.position="none") + 
  geom_jitter(width = 0.3)+theme(text = element_text(size = 8.5)) +scale_fill_manual(values=myco_longdiam)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "gray33", size=0.4))+
  theme(axis.title.x = element_text(face="bold", colour="gray33", size=19),
        axis.title.y = element_text(face="bold", colour="gray33", size=19),
        axis.text.x = element_text( colour="gray33", size=19),
        axis.text.y  = element_text( colour="gray33", size=19))+
  theme(text=element_text(size=19,  family="sans"))+ labs(tag = "H")+
  theme(plot.margin = margin(0, 0, 0, 0, "pt"));plot_evenness_Fruit_long_diameter
#ggsave("plot_evenness.TIFF", plot = plot_evenness)

tiff("richness_shannon_pd_Simpson_eveness_Fruit_long_diameter.tiff", units="in", width=20, height=20, res=300)
richness_shannon_pd_Simpson_eveness_Fruit_long_diameter.tiff=(plot_richness_Fruit_long_diameter+plot_effective_richness_Fruit_long_diameter+plot_Shannon_Fruit_long_diameter+plot_effective_Shannon_Fruit_long_diameter+plot_PD_Fruit_long_diameter+plot_Simpson_Fruit_long_diameter+plot_effective_Simpson_Fruit_long_diameter+plot_evenness_Fruit_long_diameter+ plot_layout(ncol = 2, byrow = FALSE));richness_shannon_pd_Simpson_eveness_Fruit_long_diameter.tiff
ggsave("richness_shannon_pd_Simpson_eveness_Fruit_long_diameter.tiff", plot = richness_shannon_pd_Simpson_eveness_Fruit_long_diameter.tiff)
dev.off()







#using the microeco package, we use the file convention there to compute for the assignment of microbiome functions of the seed bacterial community
##https://chiliubio.github.io/microeco_tutorial/explainable-class.html
#https://chiliubio.github.io/microeco_tutorial/basic-class.html#microtable-class
```{r,echo=FALSE,warning=FALSE}
setwd("/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project/Tomato-seed-microbiome/")
load("/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project/Tomato-seed-microbiome/ps2.RData")
library(microeco)
library(magrittr)
library(ggplot2)
library(aplot)
library(phyloseq)
source("phylo-meco.txt")
library(ggh4x)
library(ggradar)
library(ggnested)

mt=phyloseq2meco(ps2);mt
#phy=meco2phyloseq(mt);phy # microeco to phyloseq
class(mt$otu_table)
class(mt$tax_table)
class(mt$sample_table)
class(mt$phylo_tree)
#convert the microeco object
mt=microtable$new(sample_table = mt$sample_table, otu_table = mt$otu_table, tax_table = mt$tax_table, phylo_tree = mt$phylo_tree);mt
class(mt)
mt$sample_table
#The trans_env & trans_func classes are placed into the section ‘Explainable class’, as environmental factors & microbial functions can be generally applied to explain microbial community structure & assembly.

## Imported metadata file
#env_data <- read.table("env_data.txt", header=TRUE, sep="\t",row.names="SampleID");env_data

#mt$sample_table <- data.frame(mt$sample_table, env_data[rownames(mt$sample_table), ])
#mt$sample_table
#alternatively create another dataframe & append to the microeco object
#t1 <- trans_env$new(dataset = mt, add_data = env_data[, 1:12]);t1 
#t1$data_env


```
#trans_func class
#Ecological researchers are usually interested in the the funtional profiles of microbial communities, because functional or metabolic data is powerful to explain the structure & dynamics of microbial communities. As metagenomic sequencing is complicated & expensive, using amplicon sequencing data to predict functional profiles is an alternative choice. Several software are often used for this goal, such as PICRUSt (Langille et al. 2013), Tax4Fun (Aßhauer et al. 2015) & FAPROTAX (Stilianos Louca et al. 2016; S. Louca, Parfrey, & Doebeli 2016). These tools are great to be used for the prediction of functional profiles based on the prokaryotic communities from sequencing results. In addition, it is also important to obtain the traits or functions for each taxa, not just the whole profile of communities. FAPROTAX database is a collection of the traits & functions of prokaryotes based on the known research results published in books & literatures. We match the taxonomic information of prokaryotes against this database to predict the traits of prokaryotes on biogeochemical roles. The NJC19 database (Lim et al. 2020) is also available for animal-associated prokaryotic data, such as human gut microbiota. We also implement the FUNGuild (Nguyen et al. 2016) & FungalTraits (Põlme et al. 2020) databases to predict the fungal traits. The idea identifying prokaryotic traits & functional redundancy was initially inspired by our another study (C. Liu et al. 2022).
```{r,echo=FALSE,warning=FALSE}
t1 <- trans_func$new(mt)
t1$sample_table
t1$cal_spe_func(prok_database = "FAPROTAX")
t1$cal_spe_func_perc(abundance_weighted = TRUE)
# it is better to clone a dataset
tmp_mt <- clone(mt)
# transpose res_spe_func_perc to be a data.frame like taxonomic abundance
tmp <- as.data.frame(t(t1$res_spe_func_perc), check.names = FALSE)
# assign the table back to taxa_abund list for further analysis
tmp_mt$taxa_abund$func <- tmp
# select the "func" in taxa_abund list in trans_diff
t2 <- trans_diff$new(dataset = tmp_mt, method = "anova", group = "Insect_resistance", taxa_level = "func")
t2$plot_diff_abund(add_sig = T) + ggplot2::ylab("Relative abundance (%)") 

```

#https://tomsing1.github.io/blog/posts/tidyHeatmap/

```{r}

setwd("/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Figures/Clean_Figures/Carolina_figures/")
library(circlize)  # to create continuous color gradients
library(dplyr)  # wrangling tidy tables
library(grid)  # to define witdth / height of annotations using grid::unit()
library(ragg)  # to generate high quality graphics
library(RColorBrewer)  # to access predefined color palettes
library(tibble)  # to convert data.frames to tibbles
library(tidyHeatmap) # the main act: plotting & annotating heatmaps
library(readr)
library(edgeR)
data <- read_csv("Result_DAtest_insect.csv");data

#remove duplictated items
library(dplyr)
data <- data %>% distinct();data

### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
#The lfc.cutoff is set to 0.58; remember that we are working with log2 fold changes so this translates to an actual fold #change of 1.5 which is pretty reasonable. Let’s create vector that helps us identify the genes that meet our criteria:

threshold <- data$pval.adj < padj.cutoff & abs(data$logFC) > lfc.cutoff

# Volcano plot
ggplot(data) +
  geom_point(aes(x=logFC, y=-log10(data$pval.adj), colour=threshold)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
```


```{r}
resOE_df_ordered <- data[order(data$pval), ];resOE_df_ordered
resOE_df_ordered$genelabels <- rownames(resOE_df_ordered) %in% rownames(resOE_df_ordered[1:1000,])
#View(resOE_df_ordered)

ggplot(resOE_df_ordered) +
  geom_point(aes(x = logFC, y = -log10(pval.adj), colour = threshold)) +
  geom_text_repel(aes(x = logFC, y = -log10(pval.adj), label = ifelse(genelabels == T, rownames(resOE_df_ordered),""))) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

```



