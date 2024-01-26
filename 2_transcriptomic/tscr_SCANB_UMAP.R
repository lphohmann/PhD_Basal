# Script: UMAP based on ... in SCAN-B
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
               limma)
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
plot.file <- paste0(output.path,cohort,"_UMAP.pdf")
#txt.file <- paste0(output.path,cohort,"_UMAP.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()
#-------------------
start.time <- Sys.time()

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
#all.equal(anno$Sample,colnames(assay(normCounts))) # same order
# DE res
res <- loadRData(infile.4)

#######################################################################
# DE results
#######################################################################

# add bonferroni correction as well
res$Bonf.LumA <- p.adjust(res$PValue.LumA, method = "bonferroni")
res$Bonf.LumB <- p.adjust(res$PValue.LumB, method = "bonferroni")

# check results
FC_cutoff <- log2(2) # log2(3)
padj_method <- "FDR" #"Bonf"
DEGs.Basal_vs_LumA <- rownames(subset(
  res, get(paste0(padj_method,".LumA")) < 0.05 & abs(logFC.LumA) > FC_cutoff))
DEGs.Basal_vs_LumB <- rownames(subset(
  res, get(paste0(padj_method,".LumB")) < 0.05 & abs(logFC.LumB) > FC_cutoff))

#######################################################################
# UMAP based on genes that are distinct for both comparisons
#######################################################################

# Basal_vs_All: include genes that are distinct for both comparisons
normCounts.Basal_vs_All <- assay(normCounts)[
  intersect(DEGs.Basal_vs_LumA, DEGs.Basal_vs_LumB),]

# transpose data so each row is a sample
umap.dat <- umap(t(normCounts.Basal_vs_All))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$NCN.PAM50[match(rownames(umap.dat.plot),anno$Sample)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("UMAP Basal-like core DEGs n=", nrow(normCounts.Basal_vs_All)), 
     xlab = "UMAP1", ylab = "UMAP2")
legend("bottomright", legend = names(color.palette), col = color.palette, 
       pch = 16)

plot <- recordPlot()
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP based on all filtered genes
#######################################################################

# transpose data so each row is a sample
umap.dat <- umap(t(assay(normCounts)))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$NCN.PAM50[match(rownames(umap.dat.plot),anno$Sample)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("UMAP all filtered genes n=",nrow(assay(normCounts))), 
     xlab = "UMAP1", ylab = "UMAP2")
legend("topright", legend = names(color.palette), col = color.palette, 
       pch = 16)

plot <- recordPlot()
plot.list <- append(plot.list, list(plot))

#######################################################################
# UMAP based on top 5000 varying genes
#######################################################################

# transpose data so each row is a sample
umap.dat <- umap(t(assay(normCounts)))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$NCN.PAM50[match(rownames(umap.dat.plot),anno$Sample)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("UMAP all filtered genes n=",nrow(assay(normCounts))), 
     xlab = "UMAP1", ylab = "UMAP2")
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
