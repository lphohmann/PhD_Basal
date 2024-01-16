# Script: Visualizing SCAN-B differential gene expression results
# Author: Lennart Hohmann
# Date: 16.01.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "SCANB"
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2,
               ggrepel,
               tidyverse,
               DESeq2,
               edgeR,
               RColorBrewer,
               pheatmap,
               limma,
               clusterProfiler,
               org.Hs.eg.db,
               enrichplot)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/raw/gene_count_matrix-4.3.csv"
infile.4 <- "./data/SCANB/2_transcriptomic/processed/DE_result.RData"
# output paths
#outfile.1 <- #paste0(data.path,"....RData") # excel file with enrichment results
plot.file <- paste0(output.path,cohort,"_DE.pdf")
#txt.file <- paste0(output.path,cohort,"_DE.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# sample ids
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# annotation
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% unname(unlist(sampleIDs)),]

# count data
count.dat <- read.csv(infile.3)
count.dat <- count.dat[,c("gene_id",intersect(anno$GEX.assay, colnames(count.dat)))]

# correct gene and sampleIDs
count.dat$gene_id <- gsub(".*\\|", "", count.dat$gene_id)
count.dat <- count.dat[!duplicated(count.dat[c("gene_id")]), ]
rownames(count.dat) <- count.dat$gene_id
count.dat$gene_id <- NULL
names(count.dat) <- gsub("\\..*", "", colnames(count.dat))

# remove samples for which no count data is available
anno <- anno[anno$Sample %in% names(count.dat),]

# DE res
res <- loadRData(infile.4)

#######################################################################
# Visualize results
#######################################################################

# Vulcano plots
# Basal_vs_LumA
volcanoPlot.Basal_vs_LumA <- ggplot(res$table,
                                    aes(x = logFC.Basal_vs_LumA, y = -log10(FDR),
                                        color = ifelse(FDR < 0.05 & abs(logFC.Basal_vs_LumA) > 1,
                                                       "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, Log"[10]*"")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue")) +
  geom_text_repel(aes(x = logFC.Basal_vs_LumA, y = -log10(FDR), 
                      label = rownames(res$table[order(
                        -abs(res$table$logFC.Basal_vs_LumA)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res$table[order(-abs(
                    res$table$logFC.Basal_vs_LumA)), ][1:10,])
print(volcanoPlot.Basal_vs_LumA)

# Basal_vs_LumB
volcanoPlot.Basal_vs_LumB <- ggplot(res$table,
                                    aes(x = logFC.Basal_vs_LumB, y = -log10(FDR),
                                        color = ifelse(FDR < 0.05 & abs(logFC.Basal_vs_LumB) > 1,
                                                       "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, Log"[10]*"")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue")) +
  geom_text_repel(aes(x = logFC.Basal_vs_LumB, y = -log10(FDR), 
                      label = rownames(res$table[order(
                        -abs(res$table$logFC.Basal_vs_LumB)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res$table[order(-abs(
                    res$table$logFC.Basal_vs_LumB)), ][1:10,])
print(volcanoPlot.Basal_vs_LumB)

# Heatmaps
# Basal_vs_LumA
# Basal_vs_LumB

#######################################################################
# Enrichment analyses
#######################################################################

genes <- rownames(sigRes.Basal_vs_LumA)
#GO SEA
goSEA <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)
#Convert symbols to Entrez IDs
Entrez.ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db,
  drop = FALSE)$ENTREZID

#KEGG SEA
keggSEA <- enrichKEGG(
  gene = Entrez.ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

#GO GSEA
allFC <- sort(res$table$logFC.Basal_vs_LumA, decreasing = TRUE)
names(allFC) <- rownames(res$table)
goGSEA <- gseGO(
  gene = allFC,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)

#Visulization
cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 0.5)
dotplot(goSEA)
goSEA <- pairwise_termsim(goSEA)
treeplot(goSEA)
