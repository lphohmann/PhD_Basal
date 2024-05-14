# Script: Visualizing SCAN-B differential gene expression results from FPKM data
# Author: Lennart Hohmann
# Date: 16.01.2024
#TODO:  add metagene annotation tracks to heatmap
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
               VennDiagram,
               ggrepel,
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
infile.0 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.1 <- "./data/Parameters/color_palette.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.4 <- "./data/SCANB/2_transcriptomic/processed/DE_result_counts.RData"
infile.5 <- "./data/SCANB/2_transcriptomic/processed/DE_result_FPKM.RData"
infile.6 <- "./data/SCANB/2_transcriptomic/processed/Metagene_scores.RData"
# output paths
#outfile.1 <- #paste0(data.path,"....RData") # excel file with enrichment results
plot.file <- paste0(output.path,cohort,"_DE_visualization_FPKM.pdf")
#txt.file <- paste0(output.path,cohort,"_DE_visualization_FPKM.txt")
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
sample.ids <- loadRData(infile.0)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# processed count data
fpkm.dat <- loadRData(infile.3)
fpkm.dat <- fpkm.dat[colnames(fpkm.dat) %in% unname(unlist(sample.ids))]

# annotation # ggf. add prolif, SR, IR metagene scores for annatation tracks
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% unname(unlist(sample.ids)), c("Sample","NCN.PAM50")]

# add metagene scores to the annotation
mg.scores <- as.data.frame(t(loadRData(infile.6)))
mg.scores$Sample <- rownames(mg.scores)
rownames(mg.scores) <- NULL
anno <- merge(anno, mg.scores, by = "Sample", all.x = TRUE)

# apply the function to each numeric column using apply
anno[, sapply(anno, is.numeric)] <- apply(anno[, sapply(anno, is.numeric)], 2, add_labels)

# DE res
res.fpkm <- loadRData(infile.5) 
res.counts <- loadRData(infile.4) 

#######################################################################
# Check results
#######################################################################
#colnames(res.fpkm)
#colnames(res.counts)

# check results
# try different cutoffs: select genes where the differential expression 
# is statistically significant and large enough to be considered biologically meaningful
FC_cutoff <- log2(2) # log2(3)

# ordered vectors of DEGs (ordered based on logFC for later GSEA)
fpkm.degs.luma.df <- res.fpkm[res.fpkm$Basal.LumA.padj <= 0.05 & abs(res.fpkm$Basal.LumA.logFC) > FC_cutoff,]
fpkm.degs.luma <- rownames(fpkm.degs.luma.df)[order(fpkm.degs.luma.df$Basal.LumA.logFC,decreasing=TRUE)]
fpkm.degs.lumb.df <- res.fpkm[res.fpkm$Basal.LumB.padj <= 0.05 & abs(res.fpkm$Basal.LumB.logFC) > FC_cutoff,]
fpkm.degs.lumb <- rownames(fpkm.degs.lumb.df)[order(fpkm.degs.lumb.df$Basal.LumB.logFC,decreasing=TRUE)]

# fpkm.degs.luma <- rownames(
#   res.fpkm[res.fpkm$Basal.LumA.padj <= 0.05 & abs(res.fpkm$Basal.LumA.logFC) > FC_cutoff,])
# fpkm.degs.lumb <- rownames(
#   res.fpkm[res.fpkm$Basal.LumB.padj <= 0.05 & abs(res.fpkm$Basal.LumB.logFC) > FC_cutoff,])
# counts.degs.luma <- rownames(
#   res.fpkm[res.counts$FDR.LumA <= 0.05 & abs(res.counts$logFC.LumA) > FC_cutoff,])
# counts.degs.lumb <- rownames(
#   res.fpkm[res.counts$FDR.LumB <= 0.05 & abs(res.counts$logFC.LumB) > FC_cutoff,])

#length(intersect(fpkm.degs.luma,fpkm.degs.lumb)) 
#length(intersect(counts.degs.luma,counts.degs.lumb)) 

#######################################################################
# Venn diagram
#######################################################################

venn.plot <- venn.diagram(
  x = list(LumA = fpkm.degs.luma, LumB = fpkm.degs.lumb),
  category.names = c("LumA DEGs", "LumB DEGs"),
  filename = NULL
)

plot.list <- append(plot.list, list(venn.plot))

#######################################################################
# Visualize results: Vulcano plots
#######################################################################

# Basal_vs_LumA
volcanoPlot.Basal_vs_LumA <- ggplot(res.fpkm, aes(x = Basal.LumA.logFC, 
                                             y = -log10(Basal.LumA.padj),
                                             color = ifelse(
                                               Basal.LumA.padj < 0.05 & 
                                                 abs(Basal.LumA.logFC) > FC_cutoff, 
                                               "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, -Log"[10]*"")) +
  #ylim(c(0,240)) +
  geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue")) + 
  geom_text_repel(aes(x = Basal.LumA.logFC, y = -log10(Basal.LumA.padj), 
                      label = rownames(res.fpkm[order(
                        -abs(res.fpkm$Basal.LumA.logFC)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res.fpkm[order(-abs(
                    res.fpkm$Basal.LumA.logFC)), ][1:10,])

#print(volcanoPlot.Basal_vs_LumA)
plot.list <- append(plot.list, list(volcanoPlot.Basal_vs_LumA))

# Basal_vs_LumB
volcanoPlot.Basal_vs_LumB <- ggplot(res.fpkm, aes(x = Basal.LumB.logFC, 
                                                  y = -log10(Basal.LumB.padj),
                                                  color = ifelse(
                                                    Basal.LumB.padj < 0.05 & 
                                                      abs(Basal.LumB.logFC) > FC_cutoff, 
                                                    "darkred", "grey"))) +
  geom_point() +
  xlab(expression("Fold Change, Log"[2]*"")) +
  ylab(expression("Adjusted P value, -Log"[10]*"")) +
  #ylim(c(0,240)) +
  geom_vline(xintercept = c(-FC_cutoff, FC_cutoff), linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("darkred", "grey", "steelblue")) + 
  geom_text_repel(aes(x = Basal.LumB.logFC, y = -log10(Basal.LumB.padj), 
                      label = rownames(res.fpkm[order(
                        -abs(res.fpkm$Basal.LumB.logFC)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res.fpkm[order(-abs(
                    res.fpkm$Basal.LumB.logFC)), ][1:10,])

#print(volcanoPlot.Basal_vs_LumB)
plot.list <- append(plot.list, list(volcanoPlot.Basal_vs_LumB))

#######################################################################
# Visualize results: Heatmaps
#######################################################################
mg.colors <- c("<= -2"="#2e4053","-1 to -2"="#5d6d7e",
               "-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1",
               "0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00")


# Basal_vs_LumA
anno.Basal_vs_LumA <- anno[anno$NCN.PAM50 %in% c("Basal","LumA"),]
dat.Basal_vs_LumA <- fpkm.dat[
  fpkm.degs.luma, anno.Basal_vs_LumA$Sample]
sampleDist <- cor(dat.Basal_vs_LumA, method = "spearman")
plot <- pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(Proliferation = anno.Basal_vs_LumA$Mitotic_progression,
                                     SteroidResponse = anno.Basal_vs_LumA$SR,
                                     ImmuneResponse = anno.Basal_vs_LumA$IR,
                                     PAM50 = as.factor(anno.Basal_vs_LumA$NCN.PAM50),
                                     row.names = anno.Basal_vs_LumA$Sample), 
         annotation_colors = list(Proliferation = mg.colors,
                                  SteroidResponse = mg.colors,
                                  ImmuneResponse = mg.colors,
                                  PAM50 = color.palette),  
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0,
         treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

# Basal_vs_LumB LumA
anno.Basal_vs_LumB <- anno[anno$NCN.PAM50 %in% c("Basal","LumB"),]
dat.Basal_vs_LumB <- fpkm.dat[
  fpkm.degs.lumb, anno.Basal_vs_LumB$Sample]
sampleDist <- cor(dat.Basal_vs_LumB, method = "spearman")
plot <- pheatmap(sampleDist,
                 clustering_distance_rows = as.dist(1 - sampleDist),
                 clustering_distance_cols = as.dist(1 - sampleDist),
                 annotation_col = data.frame(Proliferation = anno.Basal_vs_LumB$Mitotic_progression,
                                             SteroidResponse = anno.Basal_vs_LumB$SR,
                                             ImmuneResponse = anno.Basal_vs_LumB$IR,
                                             PAM50 = as.factor(anno.Basal_vs_LumB$NCN.PAM50),
                                             row.names = anno.Basal_vs_LumB$Sample), 
                 annotation_colors = list(Proliferation = mg.colors,
                                          SteroidResponse = mg.colors,
                                          ImmuneResponse = mg.colors,
                                          PAM50 = color.palette),  
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 treeheight_row = 0,
                 treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

# Basal_vs_All: include genes that are distinct for both comparisons
dat.Basal_vs_Both <- fpkm.dat[
  intersect(fpkm.degs.luma,fpkm.degs.lumb),]
sampleDist <- cor(dat.Basal_vs_Both, method = "spearman")
plot <- pheatmap(sampleDist,
                 clustering_distance_rows = as.dist(1 - sampleDist),
                 clustering_distance_cols = as.dist(1 - sampleDist),
                 annotation_col = data.frame(Proliferation = anno$Mitotic_progression,
                                             SteroidResponse = anno$SR,
                                             ImmuneResponse = anno$IR,
                                             PAM50 = as.factor(anno$NCN.PAM50),
                                             row.names = anno$Sample), 
                 annotation_colors = list(Proliferation = mg.colors,
                                          SteroidResponse = mg.colors,
                                          ImmuneResponse = mg.colors,
                                          PAM50 = color.palette),  
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 treeheight_row = 0,
                 treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

#######################################################################
# Singluar enrichment analyses (SEA) for Basal-specific DEGs
#######################################################################

# basal specific genes and their entrez ids
genes <- intersect(fpkm.degs.luma, fpkm.degs.lumb)
#Convert symbols to Entrez IDs
Entrez.ids <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = c("ENTREZID"),
  OrgDb = org.Hs.eg.db,
  drop = FALSE)$ENTREZID

#---------------------------------
# SEA based on different databases

#GO SEA
goSEA <- enrichGO(
  gene = genes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  #subontologies: BP for Biological Process, MF for Molecular Function, and CC for Cellular Component
  ont = "BP", #BP, MF or CC
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

#Visulization
plot <- cnetplot(goSEA, colorEdge = TRUE, cex_label_gene = 0.5)
plot.list <- append(plot.list, list(plot))

plot <- dotplot(goSEA)
plot.list <- append(plot.list, list(plot))

goSEA <- pairwise_termsim(goSEA)
plot <- treeplot(goSEA)
plot.list <- append(plot.list, list(plot))


# #---------------------------------
# 
# #KEGG SEA
# keggSEA <- enrichKEGG(
#   gene = Entrez.ids,
#   organism = "hsa",
#   keyType = "ncbi-geneid",
#   pvalueCutoff = 0.05,
#   qvalueCutoff = 0.05)
# 
# #Visulization
# cnetplot(keggSEA, colorEdge = TRUE, cex_label_gene = 0.5)
# dotplot(keggSEA)
# keggSEA <- pairwise_termsim(keggSEA)
# treeplot(keggSEA)


#######################################################################
# Gene set enrichment analyses (GSEA), for LumA and LumB comparisons
#######################################################################

#GO GSEA
fpkm.degs.luma.sorted <- sort(fpkm.degs.luma.df$Basal.LumA.logFC, decreasing = TRUE)
names(fpkm.degs.luma.sorted) <- fpkm.degs.luma

fpkm.degs.lumb.sorted <- sort(fpkm.degs.lumb.df$Basal.LumB.logFC, decreasing = TRUE)
names(fpkm.degs.lumb.sorted) <- fpkm.degs.lumb

# comp 1
goGSEA <- gseGO(
  gene = fpkm.degs.luma.sorted,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  nPermSimple = 100000,
  eps=0
)

#Visulization
plot <- cnetplot(goGSEA, colorEdge = TRUE, cex_label_gene = 0.5)
plot.list <- append(plot.list, list(plot))
plot <- dotplot(goGSEA,x = "Count",showCategory = 5)
plot.list <- append(plot.list, list(plot))

goGSEA <- pairwise_termsim(goGSEA)
plot <- treeplot(goGSEA)
plot.list <- append(plot.list, list(plot))

# comp 2
goGSEA <- gseGO(
  gene = fpkm.degs.lumb.sorted,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  nPermSimple = 100000,
  eps=0
)

#Visulization
plot <- cnetplot(goGSEA, colorEdge = TRUE, cex_label_gene = 0.5)
plot.list <- append(plot.list, list(plot))
plot <- dotplot(goGSEA,x = "Count",showCategory = 5)
plot.list <- append(plot.list, list(plot))

goGSEA <- pairwise_termsim(goGSEA)
plot <- treeplot(goGSEA)
plot.list <- append(plot.list, list(plot))

#######################################################################
# Until here for now, continue tomorrow.
#######################################################################

pdf(file = plot.file, onefile = TRUE)
for(i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
}
dev.off()
