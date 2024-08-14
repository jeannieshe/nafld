# import libraries ####
library(fgsea)
library(gage)
library(EGSEAdata)
library(topGO)
library(org.Hs.eg.db)
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(limma)
library(edgeR)
library(pheatmap)
library(ggplot2)

# load pathway data and desired dataset ####
egsea.data(species = "human",returnInfo = TRUE)
kegg_list <-  kegg.pathways$human$kg.sets
gse <- getGEO("GSE168285", GSEMatrix = TRUE)
show(gse)

samp_df <- pData(gse[[1]])
rownames(samp_df) <- samp_df$description
title_df <- samp_df[,c('title')]
samp_df <- samp_df[,c('title','treatment:ch1')] #just keep the fields we want
colnames(samp_df) <- c('title','treatment')
samp_df <- samp_df %>% column_to_rownames('title')
head(samp_df)
write.csv(samp_df, file = "C:...GSE168285samp_df.csv", row.names = TRUE)

# load in expression data (raw counts) ####
getGEOSuppFiles("GSE168285")
system("gunzip C:...GSE168285/GSE168285_raw_counts.txt.gz")
count_file <- "C:...GSE168285/GSE168285_raw_counts.txt"
count_df <- read.delim2(count_file)
count_df <- count_df %>% column_to_rownames('rowname')
head(count_df)
count_df_numeric <- apply(count_df,c(1,2),as.numeric)
write.csv(count_df_numeric, file = "C:/../GSE168285count_df.csv", row.names = TRUE)

gene_names <- read.delim2("C:/..GSE168285/GSE168285_gene_annotation.txt")
count_df_annotated <- count_df # has the gene names on the last column
count_df_annotated$gene_symbols <- gene_names$GeneSymbol

# differential gene expression ####
#instead of running DESeq, since the values are decimals, and we aren't sure if they are averaged technical replicates, 
#we will run limma and prepare the dataset with voom

# FatLean_samples_boolean <- (samp_df$treatment == "Fat") | (samp_df$treatment == "Lean")
# count_df_FatLean <- count_df_numeric[,FatLean_samples_boolean]
# samp_df_FatLean <- head(samp_df, 24)

d0 <- DGEList(count_df_numeric)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
counts_voomed <- voom(d)
counts_voomed_with_design <- voom(d, design)

# create an ExpressionSet data object to prepare for limma
experimentData <- new("MIAME",name="NAFL MPS Media Conditions")
metadata <- data.frame(labelDescription=
                         c("treatment"),
                       row.names=c("treatment"))
pData <- samp_df
phenoData <- new("AnnotatedDataFrame",
                 data=pData, varMetadata=metadata)

exprs <- data.matrix(counts_voomed_with_design)
exampleSet <- ExpressionSet(assayData=exprs,
                            phenoData=phenoData,
                            experimentData=experimentData)

# creating the design matrix with only Lean as the comparison
group <- as.factor(samp_df$treatment)
group <- relevel(group, "Lean") #this is the phenotype we want to use as the reference phenotype
design <- model.matrix(~group)

# run the limma DESeq fit
fit <- lmFit(counts_voomed_with_design, design)
fit <- eBayes(fit, robust=TRUE, trend=TRUE)

# filter out for just the Fat media condition
coefficient_of_interest <- "groupFat"
top.table <- topTable(fit, coef="groupFat", n = Inf)
saveRDS(top.table, "C:/..GSE168285LeanFat_limmaFinal.rds")
# now we can use the t column in place of the stat column for running GSEA

# change ENSEMBL IDs to gene names for legibility
heatmap_dat <- subset_counts
rownames(heatmap_dat) <- mapIds(org.Hs.eg.db, keys = rownames(heatmap_dat), keytype = "ENSEMBL", column="SYMBOL")

pheatmap(heatmap_dat, annotation_col=samp_df[colnames(subset_counts),c("cell_type","treatment")], scale="row", cluster_cols=T, cluster_rows=T,
         show_rownames=T, show_colnames=T, angle_col=45)

# gene set expression analysis ####
top.table_just_t <- top.table %>% dplyr::select(t)
measurements <- top.table_just_t

# prepare genesets for GO terms
genes <- factor(x = rep(1,nrow(top.table)),levels = c(0,1))
names(genes) <- rownames(top.table)

GOobject <- new("topGOdata",ontology = "BP", allGenes = genes, annot=annFUN.org, mapping="org.Hs.eg.db", 
                ID = "ENSEMBL", nodeSize = 10)
term.genes <- genesInTerm(GOobject, GOobject@graph@nodes)
prior_genesets <- term.genes

all_groups <- names(top.table)
all_genesets <- NULL

for (i in 1:length(all_groups)) {
  coefficient_of_interest <- all_groups[i]
  top.table <- topTable(fit, coeff=coefficient_of_interest, n = Inf)
  top.table_just_t <- top.table %>% dplyr::select(t)
  measurements <- top.table_just_t
  
  print(coefficient_of_interest)
  print("running fgsea for enrichment space")
  n_permutations <- 5000
  genesets_list <-apply(measurements,MARGIN = 2,fgsea,pathways = prior_genesets,
                        minSize=10,
                        maxSize=500,
                        nperm = n_permutations)
  current_geneset <- genesets_list[[1]]
  all_genesets[[i]] <- current_geneset
  print("fgsea finished, preparing outputs")
  
}

# prepare output of gsea geneset
NES <- genesets_list[[1]]$NES
padj <- genesets_list[[1]]$padj
pval <- genesets_list[[1]]$pval

for (i in 2:length(genesets_list)) {
  
  NES <- cbind(NES,genesets_list[[i]]$NES)
  padj <- cbind(padj,genesets_list[[i]]$padj)
  pval <- cbind(pval,genesets_list[[i]]$pval)
}
colnames(NES) <- names(genesets_list)
rownames(NES) <- genesets_list[[1]]$pathway
colnames(pval) <- names(genesets_list)
rownames(pval) <- genesets_list[[1]]$pathway

fat_geneset <- genesets_list[[1]]
saveRDS(GOterm_gsea_fat_geneset, "C:/../GSE168285FatLean_GOtermGSEA.rds")


# graphing data ####
# define the colors for adjusted p-values
color_scale <- scale_color_gradient(low = "blue", high = "red")

# map GO IDs to pathway names
top_pathways$PathwayName <- mapIds(
  org.Hs.eg.db,
  keys = top_pathways$pathway,
  column = "PATH",
  keytype = "GO",
  multiVals = "first"
)

for (i in 1:length(top_pathways$pathway)) {
  top_pathways$PathwayName[i] <- Term(top_pathways$pathway[i])
}

ggplot(top_pathways, aes(x = NES, y = PathwayName, size = size, color = GeneRatio)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  color_scale +
  labs(
    x = "NES",
    y = "Pathway",
    size = "Number of Genes",
    color = "GeneRatio",
    title = "Dot Plot of Top 35 Differentially Expressed Pathways"
  ) +
  theme_minimal()

saveRDS(dotplot_FatLean, "C:/./Downloads/GSE168285DotPlot1.0.rds")

# sort the pathway data by EnrichmentScore in descending order
fat_geneset$absNES <- abs(fat_geneset$NES)
sorted_pathways <- fat_geneset[order(-fat_geneset$absNES), ]

# take top 25 pathways
top_pathways <- head(sorted_pathways, 35)
total_genes <- 58333
top_pathways$GeneRatio <- top_pathways$size / total_genes

# creating bar plot instead and ranking pathways by NES score
# reorder the levels of PathwayName based on NES scores
top_pathways$PathwayName <- factor(top_pathways$PathwayName, levels = top_pathways$PathwayName[order(top_pathways$NES)])

# plotting the bar plot with the ordered PathwayName
bar_plot <- ggplot(top_pathways, aes(x = NES, y = PathwayName, fill = GeneRatio)) +
  geom_bar(stat = "identity") +
  scale_fill_continuous() + # Assuming GeneRatio is continuous. Change scale_fill_continuous() as needed.
  labs(
    x = "NES",
    y = "Pathway",
    fill = "GeneRatio",
    title = "Top 35 Differentially Expressed Pathways"
  ) +
  theme_minimal()
ggsave("C:/./Downloads/pathwayenrichmentplot.png", plot = bar_plot, width = 7, height = 8, dpi = 300)

# preparing data for clustering in python ####
cpmed <- cpm(count_df_numeric,log=TRUE,prior.count = 2)
scaled <- scale(t(cpmed),scale=F) # this only subtracts means, does not fully zscore
# if you set scale=T then the data is zscored and fully normalized 
hist(scaled)
rownames_to_column()

write.csv(scaled, file = "C:/..GSE168285NormalizedCounts.csv", row.names = TRUE)
