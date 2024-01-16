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
infile.1 <- "./data/Parameters/color_palette.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/gene_count_matrix-4.3_processed.RData"
infile.4 <- "./data/SCANB/2_transcriptomic/processed/DE_result.RData"
# output paths
#outfile.1 <- #paste0(data.path,"....RData") # excel file with enrichment results
plot.file <- paste0(output.path,cohort,"_DE_visualization.pdf")
#txt.file <- paste0(output.path,cohort,"_DE_visualization.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load palette
color.palette <- loadRData(infile.1)[c("LumA","LumB","Basal")]

# processed count data
normCounts <- loadRData(infile.3)

# annotation # ggf. add prolif, SR, IR metagene scores for annatation tracks
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% colnames(assay(normCounts)), c("Sample","NCN.PAM50")]

# DE res
res <- loadRData(infile.4)

#######################################################################
# Check results
#######################################################################

# check results
sigRes.Basal_vs_LumA <- subset(res$table, FDR < 0.05 & abs(logFC.Basal_vs_LumA) > 1)
sigRes.Basal_vs_LumB <- subset(res$table, FDR < 0.05 & abs(logFC.Basal_vs_LumB) > 1)
nrow(sigRes.Basal_vs_LumA)
nrow(sigRes.Basal_vs_LumB)
# or use more stringed logFC cutoff e.g., log2(3) = 1.584963?

#######################################################################
# Visualize results: Vulcano plots
#######################################################################

# Basal_vs_LumA
volcanoPlot.Basal_vs_LumA <- ggplot(res$table,
                                    aes(x = logFC.Basal_vs_LumA, y = -log10(FDR),
                                        color = ifelse(FDR < 0.05 & 
                                                         abs(logFC.Basal_vs_LumA) > 1,
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
                                        color = ifelse(FDR < 0.05 & 
                                                         abs(logFC.Basal_vs_LumB) > 1,
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

#######################################################################
# Visualize results: Heatmaps
#######################################################################

# Basal_vs_LumA
anno.Basal_vs_LumA <- anno[anno$NCN.PAM50 %in% c("Basal","LumA"),]
normCounts.Basal_vs_LumA <- assay(normCounts)[
  rownames(sigRes.Basal_vs_LumA), anno.Basal_vs_LumA$Sample]
sampleDist <- cor(normCounts.Basal_vs_LumA, method = "spearman")
pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(PAM50 = as.factor(anno.Basal_vs_LumA$NCN.PAM50),
                                     row.names = anno.Basal_vs_LumA$Sample),
         annotation_colors = list(PAM50 = color.palette), 
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0)

# Basal_vs_LumB LumA
anno.Basal_vs_LumB <- anno[anno$NCN.PAM50 %in% c("Basal","LumB"),]
normCounts.Basal_vs_LumB <- assay(normCounts)[
  rownames(sigRes.Basal_vs_LumB), anno.Basal_vs_LumB$Sample]
sampleDist <- cor(normCounts.Basal_vs_LumB, method = "spearman")
pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(PAM50 = as.factor(anno.Basal_vs_LumB$NCN.PAM50),
                                     row.names = anno.Basal_vs_LumB$Sample),
         annotation_colors = list(PAM50 = color.palette), 
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0)

# Basal_vs_All: include genes that are distinct for both comparisons
normCounts.Basal_vs_All <- assay(normCounts)[
  intersect(rownames(sigRes.Basal_vs_LumA), rownames(sigRes.Basal_vs_LumB)),]
sampleDist <- cor(normCounts.Basal_vs_All, method = "spearman")
pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(PAM50 = as.factor(anno$NCN.PAM50),
                                     row.names = anno$Sample),
         annotation_colors = list(PAM50 = color.palette), # change order?
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0)


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
