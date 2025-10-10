setwd("/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project/Random_for") #<--- CHANGE ACCORDINGLY
#https://rpubs.com/michberr/randomforestmicrobe
#https://readingradio.github.io/J.nigra.Rmds/RF.GMW.Jnigra.html
#install.packages("devtools") 
#devtools::install_github("BioHPC/MegaR") 
#library(MegaR)
#MegaR() 
#Compute the prevalence at 0.025#
#Obtained from: https://f1000research.com/articles/5-1492/v1;#One of the reasons to filter by prevalence is to avoid spending much time analyzing taxa that were only rarely seen. This also turns out to be a useful filter of noise (taxa that are actually just artifacts of the data collection process)
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


#Training model and Random forest classifications
#x="/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project/Random_for/Table_rf.txt"
library(doParallel);packageVersion("doParallel")#‘1.0.17’
library(foreach);packageVersion("foreach")#‘1.5.2’
library(tidyverse);packageVersion("tidyverse")#‘2.0.0’
library(ggplot2);packageVersion("ggplot2")#‘3.5.0’
library(vegan);packageVersion("vegan")#‘2.6.4’
library(reshape2);packageVersion("reshape2")#‘1.4.4’
library(phyloseq);packageVersion("phyloseq")#‘1.46.0’
library(randomForest);packageVersion("randomForest")
library(doBy);packageVersion("doBy")#‘4.7.1.1’
library(plotrix);packageVersion("plotrix")#‘3.8.4’
library(ggpubr);packageVersion("ggpubr")#‘0.6.0’
#x="/Users/expeditoolimi/Documents/Southampton/Global_tomato_seed_Microbiome/Tomato100_Project/Random_for/Table_rf.txt"
Bacteria_otu <- read.table("Table.txt", header=TRUE,row.names="OTUID");Bacteria_otu
Bacteria_otu <- otu_table(Bacteria_otu, taxa_are_rows = TRUE);Bacteria_otu 

## Imported taxonomy
#Bacteria_tax <- read.csv("Taxonomy.csv",header=TRUE,row.names="OTUID");Bacteria_tax
Bacteria_tax <- read.table("Taxonomy.txt", header=TRUE, sep="\t",row.names="OTUID");Bacteria_tax
Bacteria_tax <- as.matrix(Bacteria_tax);Bacteria_tax
Bacteria_tax <- tax_table(Bacteria_tax);Bacteria_tax

## Imported metadata file
metadata <- read.table("Metadata.txt", header=TRUE, sep="\t",row.names="SampleID");metadata
#metadata <- read.csv("Metadata.csv",header=TRUE,row.names="SampleID");metadata

#rownames(metadata) <- metadata[,1]
metadata <- sample_data(metadata);metadata
#tree
#tree <- read.tree(file = "tree.nwk");tree
#tree <- phy_tree(tree)

all(rownames(metadata) %in% colnames(Bacteria_otu))
#Convert to phyloseq object
Bacteria <- merge_phyloseq(Bacteria_otu,Bacteria_tax, metadata);Bacteria#[ 16460 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Kingdom=="Bacteria");Bacteria# [ 6448 taxa and 84 samples ]
Bacteria = subset_taxa(Bacteria, Phylum!="CyanoBacteria");Bacteria#  [ 15588 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Order!="Chloroplast");Bacteria#  [ 15588 taxa and 1197 samples ]
Bacteria = subset_taxa(Bacteria, Family!="Mitochondria");Bacteria# [ 15588 taxa and 1197 samples ]
Bacteria <- prune_taxa(taxa_sums(Bacteria)>0, Bacteria);Bacteria# [ 15588 taxa and 1197 samples ]
#OTU_table_Bacteria = as(otu_table(Bacteria), "matrix")
#########################################################################################################
#Summarize 
microbiome::summarize_phyloseq(Bacteria)
#Number of singletons = 68"
#Sparsity = 0.995179763625127"
#Average number of reads = 22794.9197994987"
#Total number of reads = 27285519"
#Max. number of reads = 225892"
#Min. number of reads = 68"
taxonomy <- as.vector(phyloseq::tax_table(Bacteria))
sankey_phyloseq(Bacteria, fact = "Engkish_Name")
sums_Bacteria <- sample_sums(Bacteria);sums_Bacteria
sums_Bacteria <- data.frame(sums_Bacteria);sums_Bacteria
sums_Bacteria <- as.matrix(sums_Bacteria);sums_Bacteria
colSums(sums_Bacteria) #  27281011      
max(sums_Bacteria) #225892
min(sums_Bacteria) # 68
sums_Bacteria
#https://david-barnett.github.io/microViz/articles/web-only/tax-fixing.html #fixing the table
#https://rpubs.com/bioguo/680490 (Random forest classifications)
library(devtools);packageVersion("devtools")# ‘2.4.5’# Load the devtools package
#install_github("guokai8/microbial") # Install the package (Microbial package)
library(microbial);packageVersion("microbial")#‘0.0.22’####
#Plot relative abundances
#default normalize method is relative
#phy <- normalize(Bacteria, method = "relative");phy
#plotbar(phy,level="Phylum")
#plotbar(phy,level="Class")
#plotbar(phy,level="Order")
#plotbar(phy,level="Family")
#plotbar(phy,level="Family")
plotbeta(Bacteria, group="English_Name", color = mycola)#nice picture but limited colors
citation("microbial")
res <- biomarker(Bacteria,group="English_Name",ntree = 100);res#random forest classifition
write.csv(res, "classifier_cultivar.csv")
plotmarker(res,level="Genus")
plotmarker(res,level="Family")
plotmarker(res,level="Order")
plotmarker(res,level="Class")

#consider the factor (geographical region)
sample_data(Bacteria)
res_production <- biomarker(Bacteria,group="Production_site",ntree = 100);res_production#random forest classifition
write.csv(res_production, "classifier_productionsite.csv")
plotmarker(res_production,level="Genus")
plotmarker(res_production,level="Family")
plotmarker(res_production,level="Order")
plotmarker(res_production,level="Class")

res_origin <- biomarker(Bacteria,group="Country_of_Origin",ntree = 100);res_origin#random forest classifition
write.csv(res_origin, "classifier_origin.csv")
plotmarker(res_origin,level="Genus")
plotmarker(res_origin,level="Family")
plotmarker(res_origin,level="Order")
plotmarker(res_origin,level="Class")

#Subset the phyloseq object for the attributes: Fruit_shape, Fruit_taste, and Fruit_color
Bacteria_Fruit_taste <- subset_samples(Bacteria, Fruit_taste %in% c("Sweet", "Sour"));Bacteria_Fruit_taste#[ 15421 taxa and 1102 samples ]

#Random forest classification of Fruit_shape, Fruit_taste, and Fruit_color
res_Fruit_shape <- biomarker(Bacteria_Fruit_taste,group="Fruit_shape",ntree = 100);res_Fruit_shape#random forest classifition
write.csv(res_Fruit_shape, "classifier_Fruit_shape.csv")
plotmarker(res_Fruit_shape,level="Genus")
plotmarker(res_Fruit_shape,level="Family")
plotmarker(res_Fruit_shape,level="Order")
plotmarker(res_Fruit_shape,level="Class")

#Random forest classification of Fruit_shape, Fruit_taste, and Fruit_color
res_Fruit_taste <- biomarker(Bacteria_Fruit_taste,group="Fruit_taste",ntree = 100);res_Fruit_taste#random forest classifition
write.csv(res_Fruit_taste, "classifier_Fruit_taste.csv")
plotmarker(res_Fruit_taste,level="Genus")
plotmarker(res_Fruit_taste,level="Family")
plotmarker(res_Fruit_taste,level="Order")
plotmarker(res_Fruit_taste,level="Class")

#Random forest classification of Fruit_shape, Fruit_taste, and Fruit_color
res_Fruit_color <- biomarker(Bacteria_Fruit_taste,group="Fruit_color",ntree = 100);res_Fruit_color#random forest classifition
write.csv(res_Fruit_color, "classifier_Fruit_color.csv")
plotmarker(res_Fruit_color,level="Genus")
plotmarker(res_Fruit_color,level="Family")
plotmarker(res_Fruit_color,level="Order")
plotmarker(res_Fruit_color,level="Class")

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#Consider the seed cultivars
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
prunescale.Bacteria = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria <- taxa_sums(Bacteria)/nsamples(Bacteria) # average number of reads for each otu
sites.prune.Bacteria <- prune_taxa(tax.mean.Bacteria > prunescale.Bacteria, Bacteria)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria<- t(otu_table(sites.prune.Bacteria));predictors.Bacteria
dim(predictors.Bacteria)
# 1197 1655

# Make one column for our outcome/response variable 
response.Bacteria <- as.factor(sample_data(sites.prune.Bacteria)$English_Name_2);response.Bacteria
# Combine them into 1 data frame
rf.data.Bacteria <- data.frame(response.Bacteria, predictors.Bacteria);rf.data.Bacteria
str(rf.data.Bacteria)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria.classify <- randomForest(response.Bacteria ~., data = rf.data.Bacteria, ntree = 500, proximity=T, mtry=200)
  print(Bacteria.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria.classify <- randomForest(response.Bacteria ~., data = rf.data.Bacteria, ntree = 500, proximity=T, mtry=200)
  print(Bacteria.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria.classify$err.rate[500,1])
}

#Display and visualize results
# Look at the results
mean(boot.oob)#0.1176859
sd(boot.oob)#0.003960805
range(boot.oob)# 0.1086048 0.1286550
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 50 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria)[r.s, ])
write.csv(t.table, "RF.Bacteria.t.table.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria, data = rf.data.Bacteria[, c('response.Bacteria', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria
rf.data.Bacteria$response.Bacteria
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria)[rownames(hm),'English_Name_2']$English_Name_2)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs
plot1 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=13), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=22),
                                            axis.title.y = element_text(face="bold", colour="black", size=22),
                                            axis.text.x = element_text(face="bold", colour="black", size=13),
                                            axis.text.y  = element_text(face="bold", colour="black", size=23))+ labs(tag = "A",size=22);plot1
g1 <- ggplotGrob(plot1);g1

# Creat the bargraph
plot2 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=22),
                       axis.title.y = element_text(face="bold", colour="black", size=22),
                       axis.text.x = element_text(face="bold", colour="black", size=18),
                       axis.text.y  = element_text(face="bold", colour="black", size=20))+ labs(tag = "",size=20);plot2
g2 <- ggplotGrob(plot2);g2
library(ggpubr)
# Put them togethers
RF1=ggarrange(g1, g2, align='h', nrow=1, ncol=2, widths=c(4,2));RF1
ggsave("RF1.TIFF", plot = RF1);RF1
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#Consider the seed production site
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
prunescale.Bacteria = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria <- taxa_sums(Bacteria)/nsamples(Bacteria) # average number of reads for each otu
sites.prune.Bacteria <- prune_taxa(tax.mean.Bacteria > prunescale.Bacteria, Bacteria)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria<- t(otu_table(sites.prune.Bacteria));predictors.Bacteria
dim(predictors.Bacteria)
# 1197 1655

# Make one column for our outcome/response variable 
response.Bacteria <- as.factor(sample_data(sites.prune.Bacteria)$Production_site_2);response.Bacteria
# Combine them into 1 data frame
rf.data.Bacteria <- data.frame(response.Bacteria, predictors.Bacteria);rf.data.Bacteria
str(rf.data.Bacteria)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria.classify <- randomForest(response.Bacteria ~., data = rf.data.Bacteria, ntree = 500, proximity=T, mtry=200)
  print(Bacteria.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria.classify <- randomForest(response.Bacteria ~., data = rf.data.Bacteria, ntree = 500, proximity=T, mtry=200)
  print(Bacteria.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria.classify$err.rate[500,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1176859
sd(boot.oob)#0.003960805
range(boot.oob)# 0.1086048 0.1286550
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_pruductionsite.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria)[r.s, ])
write.csv(t.table, "RF.Bacteria.t.table_pruductionsite.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria, data = rf.data.Bacteria[, c('response.Bacteria', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria
rf.data.Bacteria$response.Bacteria
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria)[rownames(hm),'Production_site_2']$Production_site_2)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
write_csv(hm.melt,"hm.melt.csv")
hm.melt=read.csv(file="hm.melt.csv",header=T,check.names=FALSE,sep=",");Hm.melt
#class(Hm.melt)
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])
#hm.melt$Taxon =factor(hm.melt$Taxon, levels(hm.melt$Taxon)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs


# Creat the heatmap with the taxon names and OTUs
plot3 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=13), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=20),
                                            axis.title.y = element_text(face="bold", colour="black", size=20),
                                            axis.text.x = element_text(face="bold", colour="black", size=20),
                                            axis.text.y  = element_text(face="bold", colour="black", size=20))+ labs(tag = "");plot3
g3 <- ggplotGrob(plot3);g3

# Creat the bargraph
plot4 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#293352") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=20),
                       axis.title.y = element_text(face="bold", colour="black", size=20),
                       axis.text.x = element_text(face="bold", colour="black", size=20),
                       axis.text.y  = element_text(face="bold", colour="black", size=20))+ labs(tag = "");plot4
g4 <- ggplotGrob(plot4);g4
library(ggpubr)
# Put them togethers
RF2=ggarrange(g3, g4, align='h', nrow=1, ncol=2, widths=c(4,2));RF2
ggsave("RF2.TIFF", plot = RF2);RF2
##########################################################################################################################
##########################################################################################################################
#Subset the phyloseq object for the attributes: Fruit_shape, Fruit_taste, and Fruit_color
##########################################################################################################################
Bacteria_Fruit_color <- subset_samples(Bacteria, !is.na(Fruit_color));Bacteria_Fruit_color

#consider fruit color
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
prunescale.Bacteria_Fruit_color = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Fruit_color <- taxa_sums(Bacteria_Fruit_color)/nsamples(Bacteria_Fruit_color) # average number of reads for each otu
sites.prune.Bacteria_Fruit_color <- prune_taxa(tax.mean.Bacteria_Fruit_color > prunescale.Bacteria_Fruit_color, Bacteria_Fruit_color)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Fruit_color<- t(otu_table(sites.prune.Bacteria_Fruit_color));predictors.Bacteria_Fruit_color
dim(predictors.Bacteria_Fruit_color)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Fruit_color <- as.factor(sample_data(sites.prune.Bacteria_Fruit_color)$Fruit_color);response.Bacteria_Fruit_color
# Combine them into 1 data frame
rf.data.Bacteria_Fruit_color <- data.frame(response.Bacteria_Fruit_color, predictors.Bacteria_Fruit_color);rf.data.Bacteria_Fruit_color
str(rf.data.Bacteria_Fruit_color)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_color.classify <- randomForest(response.Bacteria_Fruit_color ~., data = rf.data.Bacteria_Fruit_color, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_color.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_color.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_color.classify <- randomForest(response.Bacteria_Fruit_color ~., data = rf.data.Bacteria_Fruit_color, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_color.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_color.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Fruit_color.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.2003207
sd(boot.oob)#0.004362374
range(boot.oob)#   0.1907173 0.2126582
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Fruit_color.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Fruit_color)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Fruit_color)[r.s, ])
write.csv(t.table, "RF.Bacteria_Fruit_color.t.table_Fruit_color.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Fruit_color, data = rf.data.Bacteria_Fruit_color[, c('response.Bacteria_Fruit_color', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Fruit_color
rf.data.Bacteria_Fruit_color$response.Bacteria_Fruit_color
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Fruit_color[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Fruit_color)[rownames(hm),'Fruit_color']$Fruit_color)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs
plot5 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot5
g5 <- ggplotGrob(plot5);g5

# Creat the bargraph
plot6 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot6
g6 <- ggplotGrob(plot6);g6
library(ggpubr)
# Put them togethers
RF3=ggarrange(g5, g6, align='h', nrow=1, ncol=2, widths=c(4,2));RF3
ggsave("RF3.TIFF", plot = RF3);RF3

##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_Fruit_shape <- subset_samples(Bacteria, !is.na(Fruit_shape));Bacteria_Fruit_shape

prunescale.Bacteria_Fruit_shape = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Fruit_shape <- taxa_sums(Bacteria_Fruit_shape)/nsamples(Bacteria_Fruit_shape) # average number of reads for each otu
sites.prune.Bacteria_Fruit_shape <- prune_taxa(tax.mean.Bacteria_Fruit_shape > prunescale.Bacteria_Fruit_shape, Bacteria_Fruit_shape)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Fruit_shape<- t(otu_table(sites.prune.Bacteria_Fruit_shape));predictors.Bacteria_Fruit_shape
dim(predictors.Bacteria_Fruit_shape)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Fruit_shape <- as.factor(sample_data(sites.prune.Bacteria_Fruit_shape)$Fruit_shape);response.Bacteria_Fruit_shape
# Combine them into 1 data frame
rf.data.Bacteria_Fruit_shape <- data.frame(response.Bacteria_Fruit_shape, predictors.Bacteria_Fruit_shape);rf.data.Bacteria_Fruit_shape
str(rf.data.Bacteria_Fruit_shape)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_shape.classify <- randomForest(response.Bacteria_Fruit_shape ~., data = rf.data.Bacteria_Fruit_shape, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_shape.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_shape.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_shape.classify <- randomForest(response.Bacteria_Fruit_shape ~., data = rf.data.Bacteria_Fruit_shape, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_shape.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_shape.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Fruit_shape.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Fruit_shape.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Fruit_shape)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Fruit_shape)[r.s, ])
write.csv(t.table, "RF.Bacteria_Fruit_shape.t.table_Fruit_shape.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Fruit_shape, data = rf.data.Bacteria_Fruit_shape[, c('response.Bacteria_Fruit_shape', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Fruit_shape
rf.data.Bacteria_Fruit_shape$response.Bacteria_Fruit_shape
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Fruit_shape[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Fruit_shape)[rownames(hm),'Fruit_shape']$Fruit_shape)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot7 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "C");plot7
g7 <- ggplotGrob(plot7);g7

# Creat the bargraph
plot8 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "grey35") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "D");plot8
g8 <- ggplotGrob(plot8);g8
library(ggpubr)
# Put them togethers
RF4=ggarrange(g7, g8, align='h', nrow=1, ncol=2, widths=c(4,2));RF4
ggsave("RF4.TIFF", plot = RF4);RF4

###########################################################################################################################
#TMV_resistance 

#consider fruit shape Bacteria_Fruit_taste
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_TMV_resistance <- subset_samples(Bacteria, !is.na(TMV_resistance));Bacteria_TMV_resistance

prunescale.Bacteria_TMV_resistance = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_TMV_resistance <- taxa_sums(Bacteria_TMV_resistance)/nsamples(Bacteria_TMV_resistance) # average number of reads for each otu
sites.prune.Bacteria_TMV_resistance <- prune_taxa(tax.mean.Bacteria_TMV_resistance > prunescale.Bacteria_TMV_resistance, Bacteria_TMV_resistance)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_TMV_resistance<- t(otu_table(sites.prune.Bacteria_TMV_resistance));predictors.Bacteria_TMV_resistance
dim(predictors.Bacteria_TMV_resistance)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_TMV_resistance <- as.factor(sample_data(sites.prune.Bacteria_TMV_resistance)$TMV_resistance);response.Bacteria_TMV_resistance
# Combine them into 1 data frame
rf.data.Bacteria_TMV_resistance <- data.frame(response.Bacteria_TMV_resistance, predictors.Bacteria_TMV_resistance);rf.data.Bacteria_TMV_resistance
str(rf.data.Bacteria_TMV_resistance)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_TMV_resistance.classify <- randomForest(response.Bacteria_TMV_resistance ~., data = rf.data.Bacteria_TMV_resistance, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_TMV_resistance.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_TMV_resistance.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_TMV_resistance.classify <- randomForest(response.Bacteria_TMV_resistance ~., data = rf.data.Bacteria_TMV_resistance, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_TMV_resistance.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_TMV_resistance.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_TMV_resistance.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_TMV_resistance.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_TMV_resistance)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_TMV_resistance)[r.s, ])
write.csv(t.table, "RF.Bacteria_TMV_resistance.t.table_TMV_resistance.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_TMV_resistance, data = rf.data.Bacteria_TMV_resistance[, c('response.Bacteria_TMV_resistance', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_TMV_resistance
rf.data.Bacteria_TMV_resistance$response.Bacteria_TMV_resistance
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_TMV_resistance[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_TMV_resistance)[rownames(hm),'TMV_resistance']$TMV_resistance)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot9 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "C");plot9
g9 <- ggplotGrob(plot9);g9

# Creat the bargraph
plot10 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#556670") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot10
g10 <- ggplotGrob(plot10);g10
library(ggpubr)
# Put them togethers
RF5=ggarrange(g9, g10, align='h', nrow=1, ncol=2, widths=c(4,2));RF5
ggsave("RF5.TIFF", plot = RF5);RF5

##########################################################################################################################
##########################################################################################################################
#Insect_resistance

#consider fruit shape Bacteria_Fruit_taste
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_Insect_resistance <- subset_samples(Bacteria, !is.na(Insect_resistance));Bacteria_Insect_resistance

prunescale.Bacteria_Insect_resistance = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Insect_resistance <- taxa_sums(Bacteria_Insect_resistance)/nsamples(Bacteria_Insect_resistance) # average number of reads for each otu
sites.prune.Bacteria_Insect_resistance <- prune_taxa(tax.mean.Bacteria_Insect_resistance > prunescale.Bacteria_Insect_resistance, Bacteria_Insect_resistance)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Insect_resistance<- t(otu_table(sites.prune.Bacteria_Insect_resistance));predictors.Bacteria_Insect_resistance
dim(predictors.Bacteria_Insect_resistance)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Insect_resistance <- as.factor(sample_data(sites.prune.Bacteria_Insect_resistance)$Insect_resistance);response.Bacteria_Insect_resistance
# Combine them into 1 data frame
rf.data.Bacteria_Insect_resistance <- data.frame(response.Bacteria_Insect_resistance, predictors.Bacteria_Insect_resistance);rf.data.Bacteria_Insect_resistance
str(rf.data.Bacteria_Insect_resistance)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(3)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Insect_resistance.classify <- randomForest(response.Bacteria_Insect_resistance ~., data = rf.data.Bacteria_Insect_resistance, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Insect_resistance.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Insect_resistance.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Insect_resistance.classify <- randomForest(response.Bacteria_Insect_resistance ~., data = rf.data.Bacteria_Insect_resistance, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Insect_resistance.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Insect_resistance.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Insect_resistance.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Insect_resistance.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Insect_resistance)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Insect_resistance)[r.s, ])
write.csv(t.table, "RF.Bacteria_Insect_resistance.t.table_Insect_resistance.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Insect_resistance, data = rf.data.Bacteria_Insect_resistance[, c('response.Bacteria_Insect_resistance', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Insect_resistance
rf.data.Bacteria_Insect_resistance$response.Bacteria_Insect_resistance
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Insect_resistance[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Insect_resistance)[rownames(hm),'Insect_resistance']$Insect_resistance)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot11 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot11
g11 <- ggplotGrob(plot11);g11

# Creat the bargraph
plot12 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#807182") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot12
g12 <- ggplotGrob(plot12);g12
library(ggpubr)
# Put them togethers
RF6=ggarrange(g11, g12, align='h', nrow=1, ncol=2, widths=c(4,2));RF6
ggsave("RF6.TIFF", plot = RF6);RF6
##########################################################################################################################
#classify by yeild range yield

#consider fruit shape Bacteria_Fruit_taste
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_yield <- subset_samples(Bacteria, !is.na(yield));Bacteria_yield

prunescale.Bacteria_yield = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_yield <- taxa_sums(Bacteria_yield)/nsamples(Bacteria_yield) # average number of reads for each otu
sites.prune.Bacteria_yield <- prune_taxa(tax.mean.Bacteria_yield > prunescale.Bacteria_yield, Bacteria_yield)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_yield<- t(otu_table(sites.prune.Bacteria_yield));predictors.Bacteria_yield
dim(predictors.Bacteria_yield)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_yield <- as.factor(sample_data(sites.prune.Bacteria_yield)$yield);response.Bacteria_yield
# Combine them into 1 data frame
rf.data.Bacteria_yield <- data.frame(response.Bacteria_yield, predictors.Bacteria_yield);rf.data.Bacteria_yield
str(rf.data.Bacteria_yield)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_yield.classify <- randomForest(response.Bacteria_yield ~., data = rf.data.Bacteria_yield, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_yield.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_yield.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_yield.classify <- randomForest(response.Bacteria_yield ~., data = rf.data.Bacteria_yield, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_yield.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_yield.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_yield.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_yield.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_yield)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_yield)[r.s, ])
write.csv(t.table, "RF.Bacteria_yield.t.table_yield.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_yield, data = rf.data.Bacteria_yield[, c('response.Bacteria_yield', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_yield
rf.data.Bacteria_yield$response.Bacteria_yield
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_yield[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_yield)[rownames(hm),'yield']$yield)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot13 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot13
g13 <- ggplotGrob(plot13);g13

# Creat the bargraph
plot14 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "bisque4") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot14
g14 <- ggplotGrob(plot14);g14
library(ggpubr)
# Put them togethers
RF6=ggarrange(g13, g14, align='h', nrow=1, ncol=2, widths=c(4,2));RF6
ggsave("RF6.TIFF", plot = RF6);RF6

########################################################################################################################################
#Weight of 1000seeds(g) predictors 

# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_seedweight_1000 <- subset_samples(Bacteria, !is.na(seedweight_1000));Bacteria_seedweight_1000

prunescale.Bacteria_seedweight_1000 = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_seedweight_1000 <- taxa_sums(Bacteria_seedweight_1000)/nsamples(Bacteria_seedweight_1000) # average number of reads for each otu
sites.prune.Bacteria_seedweight_1000 <- prune_taxa(tax.mean.Bacteria_seedweight_1000 > prunescale.Bacteria_seedweight_1000, Bacteria_seedweight_1000)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_seedweight_1000<- t(otu_table(sites.prune.Bacteria_seedweight_1000));predictors.Bacteria_seedweight_1000
dim(predictors.Bacteria_seedweight_1000)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_seedweight_1000 <- as.factor(sample_data(sites.prune.Bacteria_seedweight_1000)$seedweight_1000);response.Bacteria_seedweight_1000
# Combine them into 1 data frame
rf.data.Bacteria_seedweight_1000 <- data.frame(response.Bacteria_seedweight_1000, predictors.Bacteria_seedweight_1000);rf.data.Bacteria_seedweight_1000
str(rf.data.Bacteria_seedweight_1000)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_seedweight_1000.classify <- randomForest(response.Bacteria_seedweight_1000 ~., data = rf.data.Bacteria_seedweight_1000, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_seedweight_1000.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_seedweight_1000.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_seedweight_1000.classify <- randomForest(response.Bacteria_seedweight_1000 ~., data = rf.data.Bacteria_seedweight_1000, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_seedweight_1000.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_seedweight_1000.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_seedweight_1000.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_seedweight_1000.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_seedweight_1000)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_seedweight_1000)[r.s, ])
write.csv(t.table, "RF.Bacteria_seedweight_1000.t.table_seedweight_1000.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_seedweight_1000, data = rf.data.Bacteria_seedweight_1000[, c('response.Bacteria_seedweight_1000', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_seedweight_1000
rf.data.Bacteria_seedweight_1000$response.Bacteria_seedweight_1000
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_seedweight_1000[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_seedweight_1000)[rownames(hm),'seedweight_1000']$seedweight_1000)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot15 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot15
g15 <- ggplotGrob(plot15);g15

# Creat the bargraph
plot16 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#446455") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot16
g16 <- ggplotGrob(plot16);g16
library(ggpubr)
# Put them togethers
RF6=ggarrange(g15, g16, align='h', nrow=1, ncol=2, widths=c(4,2));RF6
ggsave("RF6.TIFF", plot = RF6);RF6
##########################################################################################################################
#ovaries

# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_ovaries <- subset_samples(Bacteria, !is.na(ovaries));Bacteria_ovaries

prunescale.Bacteria_ovaries = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_ovaries <- taxa_sums(Bacteria_ovaries)/nsamples(Bacteria_ovaries) # average number of reads for each otu
sites.prune.Bacteria_ovaries <- prune_taxa(tax.mean.Bacteria_ovaries > prunescale.Bacteria_ovaries, Bacteria_ovaries)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_ovaries<- t(otu_table(sites.prune.Bacteria_ovaries));predictors.Bacteria_ovaries
dim(predictors.Bacteria_ovaries)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_ovaries <- as.factor(sample_data(sites.prune.Bacteria_ovaries)$ovaries);response.Bacteria_ovaries
# Combine them into 1 data frame
rf.data.Bacteria_ovaries <- data.frame(response.Bacteria_ovaries, predictors.Bacteria_ovaries);rf.data.Bacteria_ovaries
str(rf.data.Bacteria_ovaries)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_ovaries.classify <- randomForest(response.Bacteria_ovaries ~., data = rf.data.Bacteria_ovaries, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_ovaries.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_ovaries.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_ovaries.classify <- randomForest(response.Bacteria_ovaries ~., data = rf.data.Bacteria_ovaries, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_ovaries.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_ovaries.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_ovaries.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_ovaries.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_ovaries)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_ovaries)[r.s, ])
write.csv(t.table, "RF.Bacteria_ovaries.t.table_ovaries.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_ovaries, data = rf.data.Bacteria_ovaries[, c('response.Bacteria_ovaries', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_ovaries
rf.data.Bacteria_ovaries$response.Bacteria_ovaries
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_ovaries[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_ovaries)[rownames(hm),'ovaries']$ovaries)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot17 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=29), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                                            axis.title.y = element_text(face="bold", colour="black", size=29),
                                            axis.text.x = element_text(face="bold", colour="black", size=29),
                                            axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot17
g17 <- ggplotGrob(plot17);g17

# Creat the bargraph
plot18 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#F4B5BD") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=29),
                       axis.title.y = element_text(face="bold", colour="black", size=29),
                       axis.text.x = element_text(face="bold", colour="black", size=29),
                       axis.text.y  = element_text(face="bold", colour="black", size=29))+ labs(tag = "");plot18
g18 <- ggplotGrob(plot18);g18
library(ggpubr)
# Put them togethers
RF6=ggarrange(g17, g18, align='h', nrow=1, ncol=2, widths=c(4,2));RF6
ggsave("RF6.TIFF", plot = RF6);RF6
##########################################################################################################################
#Fruit_long_diameter
#ovaries

# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_Fruit_long_diameter <- subset_samples(Bacteria, !is.na(Fruit_long_diameter));Bacteria_Fruit_long_diameter

prunescale.Bacteria_Fruit_long_diameter = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Fruit_long_diameter <- taxa_sums(Bacteria_Fruit_long_diameter)/nsamples(Bacteria_Fruit_long_diameter) # average number of reads for each otu
sites.prune.Bacteria_Fruit_long_diameter <- prune_taxa(tax.mean.Bacteria_Fruit_long_diameter > prunescale.Bacteria_Fruit_long_diameter, Bacteria_Fruit_long_diameter)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Fruit_long_diameter<- t(otu_table(sites.prune.Bacteria_Fruit_long_diameter));predictors.Bacteria_Fruit_long_diameter
dim(predictors.Bacteria_Fruit_long_diameter)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Fruit_long_diameter <- as.factor(sample_data(sites.prune.Bacteria_Fruit_long_diameter)$Fruit_long_diameter);response.Bacteria_Fruit_long_diameter
# Combine them into 1 data frame
rf.data.Bacteria_Fruit_long_diameter <- data.frame(response.Bacteria_Fruit_long_diameter, predictors.Bacteria_Fruit_long_diameter);rf.data.Bacteria_Fruit_long_diameter
str(rf.data.Bacteria_Fruit_long_diameter)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_long_diameter.classify <- randomForest(response.Bacteria_Fruit_long_diameter ~., data = rf.data.Bacteria_Fruit_long_diameter, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_long_diameter.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_long_diameter.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_long_diameter.classify <- randomForest(response.Bacteria_Fruit_long_diameter ~., data = rf.data.Bacteria_Fruit_long_diameter, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_long_diameter.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_long_diameter.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Fruit_long_diameter.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Fruit_long_diameter.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Fruit_long_diameter)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Fruit_long_diameter)[r.s, ])
write.csv(t.table, "RF.Bacteria_Fruit_long_diameter.t.table_Fruit_long_diameter.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Fruit_long_diameter, data = rf.data.Bacteria_Fruit_long_diameter[, c('response.Bacteria_Fruit_long_diameter', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Fruit_long_diameter
rf.data.Bacteria_Fruit_long_diameter$response.Bacteria_Fruit_long_diameter
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Fruit_long_diameter[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Fruit_long_diameter)[rownames(hm),'Fruit_long_diameter']$Fruit_long_diameter)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot19 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=15), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                                            axis.title.y = element_text(face="bold", colour="black", size=15),
                                            axis.text.x = element_text(face="bold", colour="black", size=15),
                                            axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot19
g19 <- ggplotGrob(plot19);g19

# Creat the bargraph
plot20 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "#B58900") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                       axis.title.y = element_text(face="bold", colour="black", size=15),
                       axis.text.x = element_text(face="bold", colour="black", size=15),
                       axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot20
g10 <- ggplotGrob(plot20);g10
library(ggpubr)
# Put them togethers
RF10=ggarrange(g19, g10, align='h', nrow=1, ncol=2, widths=c(4,2));RF10
ggsave("RF10.TIFF", plot = RF10);RF10

##########################################################################################################################
##########################################################################################################################
#Fruit_hori_diameter
# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_Fruit_hori_diameter <- subset_samples(Bacteria, !is.na(Fruit_hori_diameter));Bacteria_Fruit_hori_diameter

prunescale.Bacteria_Fruit_hori_diameter = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Fruit_hori_diameter <- taxa_sums(Bacteria_Fruit_hori_diameter)/nsamples(Bacteria_Fruit_hori_diameter) # average number of reads for each otu
sites.prune.Bacteria_Fruit_hori_diameter <- prune_taxa(tax.mean.Bacteria_Fruit_hori_diameter > prunescale.Bacteria_Fruit_hori_diameter, Bacteria_Fruit_hori_diameter)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Fruit_hori_diameter<- t(otu_table(sites.prune.Bacteria_Fruit_hori_diameter));predictors.Bacteria_Fruit_hori_diameter
dim(predictors.Bacteria_Fruit_hori_diameter)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Fruit_hori_diameter <- as.factor(sample_data(sites.prune.Bacteria_Fruit_hori_diameter)$Fruit_hori_diameter);response.Bacteria_Fruit_hori_diameter
# Combine them into 1 data frame
rf.data.Bacteria_Fruit_hori_diameter <- data.frame(response.Bacteria_Fruit_hori_diameter, predictors.Bacteria_Fruit_hori_diameter);rf.data.Bacteria_Fruit_hori_diameter
str(rf.data.Bacteria_Fruit_hori_diameter)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_hori_diameter.classify <- randomForest(response.Bacteria_Fruit_hori_diameter ~., data = rf.data.Bacteria_Fruit_hori_diameter, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_hori_diameter.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_hori_diameter.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_hori_diameter.classify <- randomForest(response.Bacteria_Fruit_hori_diameter ~., data = rf.data.Bacteria_Fruit_hori_diameter, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_hori_diameter.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_hori_diameter.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Fruit_hori_diameter.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Fruit_hori_diameter.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Fruit_hori_diameter)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Fruit_hori_diameter)[r.s, ])
write.csv(t.table, "RF.Bacteria_Fruit_hori_diameter.t.table_Fruit_hori_diameter.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Fruit_hori_diameter, data = rf.data.Bacteria_Fruit_hori_diameter[, c('response.Bacteria_Fruit_hori_diameter', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Fruit_hori_diameter
rf.data.Bacteria_Fruit_hori_diameter$response.Bacteria_Fruit_hori_diameter
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Fruit_hori_diameter[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Fruit_hori_diameter)[rownames(hm),'Fruit_hori_diameter']$Fruit_hori_diameter)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot19 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=15), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                                            axis.title.y = element_text(face="bold", colour="black", size=15),
                                            axis.text.x = element_text(face="bold", colour="black", size=15),
                                            axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot19
g19 <- ggplotGrob(plot19);g19

# Creat the bargraph
plot20 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill =  "#00A087B2") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                       axis.title.y = element_text(face="bold", colour="black", size=15),
                       axis.text.x = element_text(face="bold", colour="black", size=15),
                       axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot20
g10 <- ggplotGrob(plot20);g10
library(ggpubr)
# Put them togethers
RF10=ggarrange(g19, g10, align='h', nrow=1, ncol=2, widths=c(4,2));RF10
ggsave("RF10.TIFF", plot = RF10);RF10

##########################################################################################################################
##########################################################################################################################
#Fruittaste
Bacteria_Fruit_taste <- subset_samples(Bacteria, Fruit_taste %in% c("Sweet", "Sour"));Bacteria_Fruit_taste#[ 15588 taxa and 1102 samples ]
#ovaries

# Prunescale is the minimum average number of reads across samples that will be retained.
# All OTUs that do not average > prunescale across samples will be dropped.
Bacteria_Fruit_taste <- subset_samples(Bacteria, !is.na(Fruit_taste));Bacteria_Fruit_taste

prunescale.Bacteria_Fruit_taste = 0.05

# Prune out rare OTUs by mean relative abundance set by prunescale
tax.mean.Bacteria_Fruit_taste <- taxa_sums(Bacteria_Fruit_taste)/nsamples(Bacteria_Fruit_taste) # average number of reads for each otu
sites.prune.Bacteria_Fruit_taste <- prune_taxa(tax.mean.Bacteria_Fruit_taste > prunescale.Bacteria_Fruit_taste, Bacteria_Fruit_taste)

#Next, we are going to format a dataframe of predictors (OTUs) and responses (States)
# Make a dataframe of training data with OTUs as column and samples as rows
predictors.Bacteria_Fruit_taste<- t(otu_table(sites.prune.Bacteria_Fruit_taste));predictors.Bacteria_Fruit_taste
dim(predictors.Bacteria_Fruit_taste)
# 1102 1607

# Make one column for our outcome/response variable 
response.Bacteria_Fruit_taste <- as.factor(sample_data(sites.prune.Bacteria_Fruit_taste)$Fruit_taste);response.Bacteria_Fruit_taste
# Combine them into 1 data frame
rf.data.Bacteria_Fruit_taste <- data.frame(response.Bacteria_Fruit_taste, predictors.Bacteria_Fruit_taste);rf.data.Bacteria_Fruit_taste
str(rf.data.Bacteria_Fruit_taste)
# Set a random seed
set.seed(579383)
# My computer has 8 processor cores - adjust parameter of the following command as needed/desired
registerDoParallel(2)


# %dopar% implemented with foreach implements each loop call as an independent function call
# and then parallelizes all the function calls across the number of processers set above with
# registerDoParallel(). First, we need to set up a dataframe to store each bootstrap RF result.
boot.imp.s <- data.frame(predictors=NULL, Try=NULL, MeanDecreaseGini=NULL)
boot.imp.s <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_taste.classify <- randomForest(response.Bacteria_Fruit_taste ~., data = rf.data.Bacteria_Fruit_taste, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_taste.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_taste.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  cbind(imp.s, data.frame(Try=Try)) #)
  
}

# Same (and independent) as above but this time we are storing results from the
# out-of-bag error rates, so we car report the correct-incorrect classification
# rate for each group, and what types of misclassifications were made.
boot.oob <- NULL
boot.oob <- foreach (Try=1:100, .combine=rbind) %dopar% {
  
  Bacteria_Fruit_taste.classify <- randomForest(response.Bacteria_Fruit_taste ~., data = rf.data.Bacteria_Fruit_taste, ntree = 100, proximity=T, mtry=200)
  print(Bacteria_Fruit_taste.classify)
  
  # Make a data frame with predictor names and their importance
  imp.s <- importance(Bacteria_Fruit_taste.classify)
  imp.s <- data.frame(predictors = rownames(imp.s), imp.s)
  
  c(boot.oob, Bacteria_Fruit_taste.classify$err.rate[100,1])
}

#Display and visualize results

# Look at the results
mean(boot.oob)#0.1998276
sd(boot.oob)#0.003660103
range(boot.oob)#  0.1905626 0.2087114
# Create a dataframe with summary statistics from our bootstrap runs
summary.boot.s <- summaryBy(MeanDecreaseGini ~ predictors, data= boot.imp.s, FUN = c(mean, sd, std.error));summary.boot.s
write.csv(summary.boot.s,"summary.boot.s_Fruit_taste.csv")
# Look at it
head(summary.boot.s)
str(summary.boot.s)
# Order the predictor levels by importance
imp.s.sort <- arrange(summary.boot.s, desc(MeanDecreaseGini.mean));imp.s.sort
imp.s.sort$predictors <- factor(imp.s.sort$predictors, levels = imp.s.sort$predictors)
# Select the top 10 predictors
imp.s.10 <- imp.s.sort[1:5, ];imp.s.10

# What are those OTUs?
# otunames.s = otunames of predictors
# r.s = A logical containing which rows from phyloseq object were recovered as top 40 predictors
otunames.s <- imp.s.10$predictors;otunames.s
r.s <- rownames(tax_table(Bacteria_Fruit_taste)) %in% otunames.s;r.s

# Use r.s to pull out the taxonomic information for those OTUS in otunames.s
t.table <- as.data.frame(tax_table(Bacteria_Fruit_taste)[r.s, ])
write.csv(t.table, "RF.Bacteria_Fruit_taste.t.table_Fruit_taste.csv", quote=F,row.names=T)


#VISUALIZATION HERE
# Reformat the taxonomic names and paste them to OTU names for axis labels
# (highest determined taxonomic rate, get rid of special characters, etc.)
axislabels <- with(t.table, paste(rownames(t.table), gsub("[sgfocp]__","",Genus,perl=T) %>% gsub("_"," ",.,perl=T) %>% gsub("[\\s]?(unclassified|ge)[\\s]?","",.,perl=T)))
names(axislabels) <- rownames(t.table)
i <- grep("(uncultured)|([01234567890-]$)+", axislabels, perl=T)
axislabels[i]<-paste(names(axislabels)[i], t.table$Family[i])

# Now we sum the abundance of each OTU by state so we can order the
# RF predictor OTUs by which state they were most abundant in, followed
# by their importance value. This makes the heatmap easier to interpret.
sums<-summaryBy(. ~ response.Bacteria_Fruit_taste, data = rf.data.Bacteria_Fruit_taste[, c('response.Bacteria_Fruit_taste', as.character(otunames.s))], FUN = sum, keep.names = T);sums
state.organize <- sapply(as.character(otunames.s), FUN = function(x) which.max(sums[,x]));state.organize
sums$response.Bacteria_Fruit_taste
rf.data.Bacteria_Fruit_taste$response.Bacteria_Fruit_taste
# hm is a dataframe of the abundance of the predictors, standardized and log transformed
hm <- rf.data.Bacteria_Fruit_taste[,as.character(otunames.s)] %>% decostand(method='max', MARGIN=2) %>% decostand(method='log');hm
hm$sample <- gsub('_',' ',sample_data(Bacteria_Fruit_taste)[rownames(hm),'Fruit_taste']$Fruit_taste)

hm.melt <- melt(hm[c(order(state.organize), 11)]);hm.melt
names(hm.melt) <- c("Tree","Taxon","StandardizedAbundance")
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
#https://www.royfrancis.com/a-guide-to-elegant-tiled-heatmaps-in-r-2019/ (heatmaps)
# Turn predictors into a factor (re-ordered the same as the hm.melt object so they match up)
imp.s.10$predictors <- factor(imp.s.10$predictors, levels(imp.s.10$predictors)[order(state.organize)])

# Creat the heatmap with the taxon names and OTUs

plot19 <- ggplot(hm.melt, aes(x=Tree, y=Taxon, fill=StandardizedAbundance)) + 
  geom_tile() + scale_fill_gradient(low="white", high="salmon") + 
  theme(axis.text.x = element_text(angle=90, size=15), legend.position = "none", axis.text.y = element_text(hjust=0))  + 
  scale_y_discrete(labels=axislabels)+theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                                            axis.title.y = element_text(face="bold", colour="black", size=15),
                                            axis.text.x = element_text(face="bold", colour="black", size=15),
                                            axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot19
g19 <- ggplotGrob(plot19);g19

# Creat the bargraph
plot20 <- ggplot(imp.s.10, aes(x = predictors, y = MeanDecreaseGini.mean)) +
  geom_bar(stat = "identity", fill = "coral3") +
  ylab("Bootstrap Mean \u00b1 SD Gini Index") +
  theme_bw() +
  geom_errorbar( aes(x=predictors, ymin=MeanDecreaseGini.mean-MeanDecreaseGini.sd, ymax=MeanDecreaseGini.mean+MeanDecreaseGini.sd), width=0.75, alpha=0.9, size=0.5) +
  coord_flip() + theme(axis.title.x = element_text(face="bold", colour="black", size=15),
                       axis.title.y = element_text(face="bold", colour="black", size=15),
                       axis.text.x = element_text(face="bold", colour="black", size=15),
                       axis.text.y  = element_text(face="bold", colour="black", size=15))+ labs(tag = "");plot20
g10 <- ggplotGrob(plot20);g10
library(ggpubr)
# Put them togethers
RF10=ggarrange(g19, g10, align='h', nrow=1, ncol=2, widths=c(4,2));RF10
ggsave("RF10.TIFF", plot = RF10);RF10
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
