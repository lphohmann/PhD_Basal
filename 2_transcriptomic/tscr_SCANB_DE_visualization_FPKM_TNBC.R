# Script: Visualizing SCAN-B differential gene expression results from FPKM data TNBC
# Author: Lennart Hohmann
# Date: 06.05.2024
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
infile.0 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.1 <- "./data/Parameters/TNBC_color_palette.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/All_LogScaled_gex.RData"
#infile.4 <- "./data/SCANB/2_transcriptomic/processed/DE_result_counts.RData"
infile.5 <- "./data/SCANB/2_transcriptomic/processed/DE_result_FPKM_TNBC.RData"
infile.6 <- "./data/SCANB/2_transcriptomic/processed/Metagene_scores_TNBC.RData"
# output paths
#outfile.1 <- #paste0(data.path,"....RData") # excel file with enrichment results
plot.file <- paste0(output.path,cohort,"_DE_visualization_FPKM_TNBC.pdf")
#txt.file <- paste0(output.path,cohort,"_DE_visualization_FPKM_TNBC.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load palette
color.palette <- loadRData(infile.1)[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")]
sample.ids <- loadRData(infile.0)[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")]

# processed count data
fpkm.dat <- loadRData(infile.3)
fpkm.dat <- fpkm.dat[colnames(fpkm.dat) %in% unname(unlist(sample.ids))]

# annotation # ggf. add prolif, SR, IR metagene scores for annatation tracks
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% unname(unlist(sample.ids)), c("Sample","NCN.PAM50")]
anno$Subtype <- ifelse(anno$Sample %in% unlist(sample.ids["TNBC_NonBasal"]),
                       "TNBC_NonBasal", 
                       ifelse(anno$Sample %in% unlist(sample.ids["TNBC_Basal"]),
                              "TNBC_Basal", "ERpHER2n_Basal"))

# add metagene scores to the annotation
mg.scores <- as.data.frame(t(loadRData(infile.6)))
mg.scores$Sample <- rownames(mg.scores)
rownames(mg.scores) <- NULL
anno <- merge(anno, mg.scores, by = "Sample", all.x = TRUE)

# apply the function to each numeric column using apply
anno[, sapply(anno, is.numeric)] <- apply(anno[, sapply(anno, is.numeric)], 2, add_labels)

# DE res
res.fpkm <- loadRData(infile.5) 
#res.counts <- loadRData(infile.4) 

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
fpkm.degs.bas.df <- res.fpkm[res.fpkm$Basal.tnbc.basal.padj <= 0.05 & abs(res.fpkm$Basal.tnbc.basal.logFC) > FC_cutoff,]
fpkm.degs.bas <- rownames(fpkm.degs.bas.df)[order(fpkm.degs.bas.df$Basal.tnbc.basal.logFC,decreasing=TRUE)]

fpkm.degs.nonbas.df <- res.fpkm[res.fpkm$Basal.tnbc.nonbasal.padj <= 0.05 & abs(res.fpkm$Basal.tnbc.nonbasal.logFC) > FC_cutoff,]
fpkm.degs.nonbas <- rownames(fpkm.degs.nonbas.df)[order(fpkm.degs.nonbas.df$Basal.tnbc.nonbasal.logFC,decreasing=TRUE)]

length(intersect(fpkm.degs.bas,fpkm.degs.nonbas)) 

#######################################################################
# Venn diagram
#######################################################################

plot <- venn.diagram(
  x = list(TNBC_Basal = fpkm.degs.bas, TNBC_NonBasal = fpkm.degs.nonbas),
  category.names = c("TNBC_Basal DEGs", "TNBC_NonBasal DEGs"),
  filename = NULL
)

plot.list <- append(plot.list, list(plot))

#######################################################################
# Visualize results: Vulcano plots
#######################################################################
# Basal.tnbc.basal.logFC
# Basal.tnbc.nonbasal.logFC
# Basal.tnbc.basal.padj
# Basal.tnbc.nonbasal.padj

# Basal_vs_tnbc basal
volcanoPlot.Basal_vs_TNBCbasal <- ggplot(res.fpkm, aes(x = Basal.tnbc.basal.logFC, 
                                             y = -log10(Basal.tnbc.basal.padj),
                                             color = ifelse(
                                               Basal.tnbc.basal.padj < 0.05 & 
                                                 abs(Basal.tnbc.basal.logFC) > FC_cutoff, 
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
  geom_text_repel(aes(x = Basal.tnbc.basal.logFC, y = -log10(Basal.tnbc.basal.padj), 
                      label = rownames(res.fpkm[order(
                        -abs(res.fpkm$Basal.tnbc.basal.logFC)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res.fpkm[order(-abs(
                    res.fpkm$Basal.tnbc.basal.logFC)), ][1:10,])

#print(volcanoPlot.Basal_vs_TNBCbasal)
plot.list <- append(plot.list, list(volcanoPlot.Basal_vs_TNBCbasal))

# Basal_vs_ tnbc non basal
volcanoPlot.Basal_vs_TNBCnonbasal <- ggplot(res.fpkm, aes(x = Basal.tnbc.nonbasal.logFC, 
                                                  y = -log10(Basal.tnbc.nonbasal.padj),
                                                  color = ifelse(
                                                    Basal.tnbc.nonbasal.padj < 0.05 & 
                                                      abs(Basal.tnbc.nonbasal.logFC) > FC_cutoff, 
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
  geom_text_repel(aes(x = Basal.tnbc.nonbasal.logFC, y = -log10(Basal.tnbc.nonbasal.padj), 
                      label = rownames(res.fpkm[order(
                        -abs(res.fpkm$Basal.tnbc.nonbasal.logFC)), ][1:10,]),
                      size = 2, color = "steelblue"),
                  data = res.fpkm[order(-abs(
                    res.fpkm$Basal.tnbc.nonbasal.logFC)), ][1:10,])

#print(volcanoPlot.Basal_vs_TNBCnonbasal)
plot.list <- append(plot.list, list(volcanoPlot.Basal_vs_TNBCnonbasal))

#######################################################################
# Visualize results: Heatmaps
#######################################################################
mg.colors <- c("<= -2"="#2e4053","-1 to -2"="#5d6d7e",
               "-0.5 to -1"="#aeb6bf","-0.5 to 0.5"="#ecf0f1",
               "0.5 to 1"="#edbb99","1 to 2"="#dc7633",">= 2"="#ba4a00")

# Basal_vs_TNBCbasal
anno.Basal_vs_TNBCbasal <- anno[anno$Subtype %in% c("ERpHER2n_Basal","TNBC_Basal"),]
dat.Basal_vs_TNBCbasal <- fpkm.dat[
  fpkm.degs.bas, anno.Basal_vs_TNBCbasal$Sample]
sampleDist <- cor(dat.Basal_vs_TNBCbasal, method = "spearman")
plot <- pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(Proliferation = anno.Basal_vs_TNBCbasal$Mitotic_progression,
                                     SteroidResponse = anno.Basal_vs_TNBCbasal$SR,
                                     ImmuneResponse = anno.Basal_vs_TNBCbasal$IR,
                                     Subtype = as.factor(anno.Basal_vs_TNBCbasal$Subtype),
                                     row.names = anno.Basal_vs_TNBCbasal$Sample), 
         annotation_colors = list(Proliferation = mg.colors,
                                  SteroidResponse = mg.colors,
                                  ImmuneResponse = mg.colors,
                                  Subtype = color.palette),  
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0,
         treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

# Basal_vs_TNBCnonbasal
anno.Basal_vs_TNBCnonbasal <- anno[anno$Subtype %in% c("ERpHER2n_Basal","TNBC_NonBasal"),]
dat.Basal_vs_TNBCnonbasal <- fpkm.dat[
  fpkm.degs.nonbas, anno.Basal_vs_TNBCnonbasal$Sample]
sampleDist <- cor(dat.Basal_vs_TNBCnonbasal, method = "spearman")
plot <- pheatmap(sampleDist,
                 clustering_distance_rows = as.dist(1 - sampleDist),
                 clustering_distance_cols = as.dist(1 - sampleDist),
                 annotation_col = data.frame(Proliferation = anno.Basal_vs_TNBCnonbasal$Mitotic_progression,
                                             SteroidResponse = anno.Basal_vs_TNBCnonbasal$SR,
                                             ImmuneResponse = anno.Basal_vs_TNBCnonbasal$IR,
                                             Subtype = as.factor(anno.Basal_vs_TNBCnonbasal$Subtype),
                                             row.names = anno.Basal_vs_TNBCnonbasal$Sample), 
                 annotation_colors = list(Proliferation = mg.colors,
                                          SteroidResponse = mg.colors,
                                          ImmuneResponse = mg.colors,
                                          Subtype = color.palette),  
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 treeheight_row = 0,
                 treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

# Basal_vs_All: include genes that are distinct for both comparisons
dat.Basal_vs_Both <- fpkm.dat[
  intersect(fpkm.degs.bas,fpkm.degs.nonbas),]
sampleDist <- cor(dat.Basal_vs_Both, method = "spearman")
plot <- pheatmap(sampleDist,
                 clustering_distance_rows = as.dist(1 - sampleDist),
                 clustering_distance_cols = as.dist(1 - sampleDist),
                 annotation_col = data.frame(Proliferation = anno$Mitotic_progression,
                                             SteroidResponse = anno$SR,
                                             ImmuneResponse = anno$IR,
                                             Subtype = as.factor(anno$Subtype),
                                             row.names = anno$Sample), 
                 annotation_colors = list(Proliferation = mg.colors,
                                          SteroidResponse = mg.colors,
                                          ImmuneResponse = mg.colors,
                                          Subtype = color.palette),  
                 show_rownames = FALSE,
                 show_colnames = FALSE,
                 treeheight_row = 0,
                 treeheight_col = 0)
plot.list <- append(plot.list, list(plot))

#######################################################################
# Singluar enrichment analyses (SEA) for Basal-specific DEGs
#######################################################################

# basal specific genes and their entrez ids
genes <- intersect(fpkm.degs.bas, fpkm.degs.nonbas)
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

plot <- dotplot(goSEA,x = "Count",showCategory = 5)
plot.list <- append(plot.list, list(plot))

plot <- dotplot(goSEA,x = "GeneRatio",showCategory = 5)
plot.list <- append(plot.list, list(plot))
#print(plot)
goSEA <- pairwise_termsim(goSEA)
plot <- treeplot(goSEA)
plot.list <- append(plot.list, list(plot))

#######################################################################
# Gene set enrichment analyses (GSEA)
#######################################################################

#GO GSEA

fpkm.degs.bas.sorted <- sort(fpkm.degs.bas.df$Basal.tnbc.basal.logFC, decreasing = TRUE)
names(fpkm.degs.bas.sorted) <- fpkm.degs.bas

fpkm.degs.nonbas.sorted <- sort(fpkm.degs.nonbas.df$Basal.tnbc.nonbasal.logFC, decreasing = TRUE)
names(fpkm.degs.nonbas.sorted) <- fpkm.degs.nonbas

# comp 1
goGSEA <- gseGO(
 gene = fpkm.degs.bas.sorted,
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


# goGSEA <- pairwise_termsim(goGSEA)
# plot <- treeplot(goGSEA)
# plot.list <- append(plot.list, list(plot))

# comp 2
goGSEA <- gseGO(
  gene = fpkm.degs.nonbas.sorted,
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


# goGSEA <- pairwise_termsim(goGSEA)
# plot <- treeplot(goGSEA)
# plot.list <- append(plot.list, list(plot))


#######################################################################
# Until here for now, continue tomorrow.
#######################################################################

pdf(file = plot.file, onefile = TRUE)
for(i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
}
dev.off()
