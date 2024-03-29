# Script: Differential gene expression analysis in SCAN-B based on count data
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
               clusterProfiler)
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
# output paths
outfile.1 <- paste0(data.path,"DE_result_counts.RData")
outfile.2 <- paste0(data.path,"gene_count_matrix-4.3_processed.RData")
#plot.file <- paste0(output.path,cohort,"_DE.pdf")
#txt.file <- paste0(output.path,cohort,"_DE.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
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

#######################################################################
# Define count and sample tables
#######################################################################

# countTable: samples column, genes rows
countTable <- count.dat
# sampleTable: sampleID PAM50 columns
sampleTable <- anno[c("Sample","NCN.PAM50","LibraryProtocol")]
sampleTable$NCN.PAM50 <- as.factor(sampleTable$NCN.PAM50)
sampleTable$LibraryProtocol <- as.factor(sampleTable$LibraryProtocol)
# critical that count matrix sample data are in the same order and format
countTable <- countTable[, sampleTable$Sample]
identical(colnames(countTable),sampleTable$Sample) # TRUE
#sampleTable[1:5,]
#countTable[1:5,1:5]

#######################################################################
# Filter low counts
#######################################################################
# remove non- or very low-expressed genes, which are not interesting to 
# include in the analysis. Genes at the lower end of the count
# distribution are susceptible to noise since a single read can make a 
# large difference in terms of mean gene expression. 
# Inspect the CPM histogram and choose a value that gets rid of the left 
# peak in the distribution.This peak represents non- and/or very 
# low-expressed genes. Check that a reasonable number of genes remain 
# for the downstream analysis.

# calc. mean log2 counts per million (CPM) for gene 
# cpm() normalizes for different sequencing depths to make samples comparable
meanLog2CPM <- rowMeans(log2(cpm(countTable) + 1)) 
# distribution of the mean log2 CPM values across genes
hist(meanLog2CPM) 
# how many genes have a mean log2 CPM value <= 1 indicating low expression
sum(meanLog2CPM <= 1) #5451
# filters out genes with low expression
countTable <- countTable[meanLog2CPM > 1, ]
dim(countTable)

#######################################################################
# QC and normalizing: Note that the normalized counts should only be 
# used for QC, not the statistical analysis.
#######################################################################
# correct the counts for sequencing depth (number of counts in each sample) 
# and to remove the effect of heteroscedasticity (that variance increases 
# with increasing mean expression).

# prepare data for QC
dds <- DESeqDataSetFromMatrix(as.matrix(countTable),
                              # 
                              design = ~ LibraryProtocol + NCN.PAM50, 
                              colData = sampleTable)

print(dds)

# normalize using variance stabilizing transformation 
normCounts <- vst(dds, 
                  blind = TRUE) # performed without considering sample groups
# extract the normalized counts
assay(normCounts)[1:5, 1:5]

# distribution of variance-stabilized counts
hist(assay(normCounts)) #high requency peaks of 7000000, due to higher sample size?

#######################################################################
# Visualization of overall gene expression patterns
#######################################################################

# Sample heatmap
# sampleDist <- cor(assay(normCounts), method = "spearman")
# sampleColor <- brewer.pal(4, "Accent")[1:3]
# names(sampleColor) <- levels(sampleTable$NCN.PAM50)
# pheatmap(sampleDist,
#          clustering_distance_rows = as.dist(1 - sampleDist),
#          clustering_distance_cols = as.dist(1 - sampleDist),
#          annotation_col = data.frame(PAM50 = sampleTable$NCN.PAM50,
#                                      row.names = sampleTable$Sample),
#          annotation_colors = list(PAM50 = sampleColor),
#          show_rownames = FALSE,
#          show_colnames = FALSE,
#          treeheight_row = 0,
#          filename = "./output/2_transcriptomic/heatmap.png")

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
#   mapping = aes(x = PC1, y = PC2, color = Disease, label = Sample)) +
#   geom_point(size = 3) +
#   geom_text_repel(size = 4) +
#   labs(x = paste0("PC1 (", varExp[1], " %)"), y = paste0("PC2 (", varExp[2], " %)")) +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 12), legend.text = element_text(size = 10)) +
#   scale_color_manual(values = brewer.pal(3, "Accent"))
# print(pcaPlot)
#ggsave("./output/2_transcriptomic/pca_all.png")

#Remove outliers?
# countTable <- subset(countTable, select = -C6_T_FF)
# sampleTable <- subset(sampleTable, subset = sample != "C6_T_FF")
# sampleTable <- droplevels(sampleTable)
# colnames(countTable)

#######################################################################
# Differential gene expression analysis 
#######################################################################

# Step 1: Define design matrix
designMatrix <- model.matrix(~ LibraryProtocol + NCN.PAM50, data = sampleTable)

rownames(designMatrix) <- colnames(countTable)
designMatrix[1:3,]
# Step 2: Define contrast matrix (Basal vs. LumA and Basal vs. LumB)
# contrastMatrix <- makeContrasts(Basal_vs_LumA = Intercept - NCN.PAM50LumA,
#                                 Basal_vs_LumB = Intercept - NCN.PAM50LumB,
#                                 levels = designMatrix)
# head(contrastMatrix)

#------
# Step 2.5: Correct for different LibraryProtocols (only for plotting)
assay(normCounts) <- removeBatchEffect(assay(normCounts), 
                                       normCounts$LibraryProtocol, 
                                       design=designMatrix)

# save processed count data
save(normCounts, file=outfile.2)
#------

# Step 3: Fit model
dge <- DGEList(countTable, group=sampleTable$NCN.PAM50) # creates a DGEList object
dge <- calcNormFactors(dge) # calculate normalization factors for the raw count data
dge <- estimateDisp(dge, designMatrix, robust = TRUE) # estimate dispersion (variability)
# dge object contains normalized counts and dispersion estimates
fit <- glmQLFit(dge, designMatrix, robust = TRUE) # fit generalized linear model to the data.

# Step 4: Perform hypothesis testing for both groups
qlf.Basal_vs_LumA <- glmQLFTest(fit, coef=4)
res.Basal_vs_LumA <- topTags(qlf.Basal_vs_LumA, n = nrow(countTable))
summary(decideTests(qlf.Basal_vs_LumA))

qlf.Basal_vs_LumB <- glmQLFTest(fit, coef=5)
res.Basal_vs_LumB <- topTags(qlf.Basal_vs_LumB, n = nrow(countTable))
summary(decideTests(qlf.Basal_vs_LumB))

# Step 5: Compile both comps into one df
df.a <- res.Basal_vs_LumA$table[c("logFC","PValue","FDR")]
names(df.a) <- c("logFC.LumA","PValue.LumA","FDR.LumA")
df.a <- df.a[order(rownames(df.a)), ]

df.b <- res.Basal_vs_LumB$table[c("logFC","PValue","FDR")]
names(df.b) <- c("logFC.LumB","PValue.LumB","FDR.LumB")
df.b <- df.b[order(rownames(df.b)), ]
#identical(rownames(df.a),rownames(df.b)) # TRUE
res <- cbind(df.a, df.b)

# Step 6: Save results
save(res, file=outfile.1)

# check results
sigRes.Basal_vs_LumA <- subset(res, FDR.LumA < 0.05 & abs(logFC.LumA) > 1)
sigRes.Basal_vs_LumB <- subset(res, FDR.LumB < 0.05 & abs(logFC.LumB) > 1)
nrow(sigRes.Basal_vs_LumA)
nrow(sigRes.Basal_vs_LumB)
