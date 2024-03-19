# Script: UMAP based on FPKM and Count data in SCAN-B; and PCA
# Author: Lennart Hohmann
# Date: 25.01.2024
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
#source("./scripts/3_WGS/src/wgs_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(umap,
               DESeq2,
               edgeR,
               limma,
               ggplot2)
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
infile.4 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.5 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_UMAP.pdf")
#txt.file <- paste0(output.path,cohort,"_UMAP.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# sample IDs
sampleIDs <- unname(unlist(loadRData(infile.5)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]))

# load palette
color.palette <- loadRData(infile.1)[c("LumA","LumB","Basal")]

# processed count data
normCounts <- loadRData(infile.3)

# processed count data
FPKM.dat <- loadRData(infile.4)
FPKM.dat <- FPKM.dat[sampleIDs]

# annotation # ggf. add prolif, SR, IR metagene scores for annatation tracks
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE, c("Sample","NCN.PAM50")]
anno <- anno[anno$Sample %in% sampleIDs, ]
#anno <- anno[anno$Sample %in% colnames(assay(normCounts)), c("Sample","NCN.PAM50")] # here
#all.equal(anno$Sample,colnames(assay(normCounts))) # same order

#######################################################################
# UMAP based on FPKM data; all genes
#######################################################################

# checks same order
all.equal(anno$Sample,colnames(FPKM.dat)) # no
FPKM.dat <- FPKM.dat[anno$Sample]
all.equal(anno$Sample,colnames(FPKM.dat)) # yes

# each row is a sample
umap.dat <- umap(t(FPKM.dat))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$NCN.PAM50[match(rownames(umap.dat.plot),anno$Sample)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("FPKM UMAP all filtered genes n=",nrow(FPKM.dat)), 
     xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = names(color.palette), col = color.palette, 
       pch = 16)

plot <- recordPlot()
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP based on count data; all filtered genes
#######################################################################

anno <- anno[anno$Sample %in% colnames(assay(normCounts)), c("Sample","NCN.PAM50")] # here
all.equal(anno$Sample,colnames(assay(normCounts))) # same order

# transpose data so each row is a sample
umap.dat <- umap(t(assay(normCounts)))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$NCN.PAM50[match(rownames(umap.dat.plot),anno$Sample)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("Count UMAP all filtered genes n=",nrow(assay(normCounts))), 
     xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = names(color.palette), col = color.palette, 
       pch = 16)

plot <- recordPlot()
plot.list <- append(plot.list, list(plot))

#######################################################################
# PCA based on count data
#######################################################################

#Sample PCA
pcaRes <- prcomp(t(assay(normCounts)))
varExp <- round(pcaRes$sdev^2 / sum(pcaRes$sdev^2) * 100)
pcaDF <- data.frame(
  PC1 = pcaRes$x[, 1],
  PC2 = pcaRes$x[, 2],
  Group = anno$NCN.PAM50,
  Sample = anno$Sample)

# Scatter plot of PC1 and PC2
plot(pcaDF$PC1, pcaDF$PC2, 
     col = color.palette[factor(pcaDF$Group, levels = names(color.palette))],
     pch = 16, cex = 2, 
     xlab = paste0("PC1 (", varExp[1], " %)"), 
     ylab = paste0("PC2 (", varExp[2], " %)"), 
     main = paste0("Count PCA all filtered genes n=",nrow(assay(normCounts))), )

# Add legend
legend("topright", legend = names(color.palette), col = color.palette, 
       pch = 16)

plot <- recordPlot()
plot.list <- append(plot.list, list(plot))

#######################################################################
#######################################################################

pdf(file = plot.file, onefile = TRUE)
for(i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()
