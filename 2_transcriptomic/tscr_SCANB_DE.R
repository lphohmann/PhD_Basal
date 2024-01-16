# Script: Differential gene expression analysis in SCAN-B
# Author: Lennart Hohmann
# Date: 15.01.2024
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
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
# output paths
outfile.1 <- paste0(data.path,"DE_result.RData")
plot.file <- paste0(output.path,cohort,"_DE.pdf")
txt.file <- paste0(output.path,cohort,"_DE.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# sample ids
sampleIDs <- loadRData("./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData")[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# annotation
anno <- loadRData("./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% unname(unlist(sampleIDs)),]

# count data
count.dat <- read.csv("./data/SCANB/2_transcriptomic/raw/gene_count_matrix-4.3.csv")
count.dat <- count.dat[,c("gene_id",intersect(anno$GEX.assay, colnames(count.dat)))]

# correct gene and sampleIDs
count.dat$gene_id <- gsub(".*\\|", "", count.dat$gene_id)
count.dat <- count.dat[!duplicated(count.dat[c("gene_id")]), ]
rownames(count.dat) <- count.dat$gene_id
count.dat$gene_id <- NULL
names(count.dat) <- gsub("\\..*", "", colnames(count.dat))

# remove samples for which no count data is available
anno <- anno[anno$Sample %in% names(count.dat),]

#######################################################################
# Define count and sample tables
#######################################################################

# countTable: samples column, genes rows
countTable <- count.dat
# sampleTable: sampleID PAM50 columns
sampleTable <- anno[c("Sample","NCN.PAM50")]
sampleTable$NCN.PAM50 <- as.factor(sampleTable$NCN.PAM50)

#######################################################################
# Filter low counts
#######################################################################

# calc. mean log2 counts per million (CPM) for gene 
# cpm() normalizes for different sequencing depths to make samples comparable
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1)) 
hist(meanLog2CPM)
sum(meanLog2CPM <= 1) #This line calculates the number of rows in the count table where the mean log2 CPM is less than or equal to 1. It counts how many rows have a mean log2 CPM value indicating low expression.
countTable <- countTable[meanLog2CPM > 1, ] #This line subsets the original countTable by keeping only those rows where the mean log2 CPM is greater than 1. It filters out rows with low expression based on the previously calculated mean log2 CPM values.
dim(countTable)

#######################################################################
# QC and normalizing
#######################################################################

#Prepare data for QC
dds <- DESeqDataSetFromMatrix(as.matrix(countTable),
                              design = ~ 0 + NCN.PAM50,
                              colData = sampleTable)

print(dds)

#Normalize
normCounts <- vst(dds, blind = TRUE) # make sure to handle the heteroscedastic data by normalising away the mean/variance relationship (higher count values associated with higher variance)
assay(normCounts)[1:5, 1:5]

#Distribution
hist(assay(normCounts))

#######################################################################
# Visualization of overall gene expression patterns
#######################################################################

#Sample heatmap
sampleDist <- cor(assay(normCounts), method = "spearman")
sampleColor <- brewer.pal(4, "Accent")[1:3]
names(sampleColor) <- levels(sampleTable$NCN.PAM50)
pheatmap(sampleDist,
         clustering_distance_rows = as.dist(1 - sampleDist),
         clustering_distance_cols = as.dist(1 - sampleDist),
         annotation_col = data.frame(PAM50 = sampleTable$NCN.PAM50,
                                     row.names = sampleTable$Sample),
         annotation_colors = list(PAM50 = sampleColor),
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0,
         filename = "./output/2_transcriptomic/heatmap.png")

#Sample PCA
# pcaRes <- prcomp(t(assay(normCounts)))
# varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
# pcaDF <- data.frame(
#   PC1 = pcaRes$x[, 1],
#   PC2 = pcaRes$x[, 2],
#   Disease = sampleTable$NCN.PAM50,
#   Sample = sampleTable$Sample)
# pcaPlot <- ggplot(
#   data = pcaDF,
#   mapping = aes(x = PC1, y = PC2, color = NCN.PAM50, label = Sample)) +
#   geom_point(size = 3) +
#   geom_text_repel(size = 4) +
#   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
#   scale_color_manual(values = brewer.pal(3, "Accent"))
# print(pcaPlot)
# ggsave("./output/2_transcriptomic/pca_all.png")

#Remove outliers?
# countTable <- subset(countTable, select = -C6_T_FF)
# sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# sampleTable <- droplevels(sampleTable)
# colnames(countTable)


#######################################################################
# Differential gene expression analysis
#######################################################################

#Step 1: Define design matrix
designMatrix <- model.matrix(~ 0 + NCN.PAM50, data = sampleTable)
designMatrix[1:3,1:3]

#Step 2: Define contrast matrix (Basal vs. LumA and Basal vs. LumB)
contrastMatrix <- makeContrasts(Basal_vs_LumA = NCN.PAM50Basal - NCN.PAM50LumA,
                                Basal_vs_LumB = NCN.PAM50Basal - NCN.PAM50LumB,
                                levels = designMatrix)
head(contrastMatrix)

#Step 3: Fit model
dge <- DGEList(countTable)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, designMatrix, robust = TRUE)
fit <- glmQLFit(dge, designMatrix, robust = TRUE)

#Step 4: Perform hypothesis testing
res <- glmQLFTest(fit, contrast = contrastMatrix)
res <- topTags(res, n = nrow(countTable))
sigRes.Basal_vs_LumA <- subset(res$table, FDR < 0.05 & abs(logFC.Basal_vs_LumA) > 1)
sigRes.Basal_vs_LumB <- subset(res$table, FDR < 0.05 & abs(logFC.Basal_vs_LumB) > 1)
nrow(sigRes.Basal_vs_LumA)
nrow(sigRes.Basal_vs_LumB)
save(res,file=outfile.1)

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
