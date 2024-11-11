# Script: UMAP based on FPKM and Count data in METABRIC
# Author: Lennart Hohmann
# Date: 08.11.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "METABRIC"
#-------------------
# packages
source("./scripts/src/general_functions.R")
#source("./scripts/3_WGS/src/wgs_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(umap)#,DESeq2,edgeR,limma)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("DESeq2")
#-------------------
# set/create output directories
# for plots
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/Parameters/color_palette.RData"
infile.2 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
#infile.3 <- "./data/METABRIC/2_transcriptomic/processed/gene_count_matrix-4.3_processed.RData"
infile.4 <- "./data/METABRIC/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.5 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
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

FPKM.dat <- loadRData(infile.4)
FPKM.dat <- FPKM.dat[,colnames(FPKM.dat) %in% sampleIDs]
#na_data_rows <- FPKM.dat[rowSums(is.na(FPKM.dat)) > 0, ] # 6 na values
FPKM.dat <- na.omit(FPKM.dat)

# annotation # ggf. add prolif, SR, IR metagene scores for annatation tracks
anno <- loadRData(infile.2)
anno <- anno[, c("METABRIC_ID","PAM50")]
anno <- anno[anno$METABRIC_ID %in% colnames(FPKM.dat), ]
#anno <- anno[anno$Sample %in% colnames(assay(normCounts)), c("Sample","NCN.PAM50")] # here
#all.equal(anno$Sample,colnames(assay(normCounts))) # same order

#######################################################################
# UMAP based on FPKM data; all genes
#######################################################################

# checks same order
all.equal(anno$METABRIC_ID,colnames(FPKM.dat)) # no
FPKM.dat <- FPKM.dat[anno$METABRIC_ID]
all.equal(anno$METABRIC_ID,colnames(FPKM.dat)) # yes

# each row is a sample
umap.dat <- umap(t(FPKM.dat)) # cont here
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$PAM50[match(row.names(umap.dat.plot),anno$METABRIC_ID)]

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
# UMAP based on FPKM data; 5000 most varying
#######################################################################

# checks same order
all.equal(anno$METABRIC_ID,colnames(FPKM.dat)) # yes
gene.vars <- apply(FPKM.dat,1,var,na.rm = TRUE)
top.genes <- names(gene.vars[order(-gene.vars)][1:5000])
FPKM.dat.filt <- FPKM.dat[top.genes,] # filt here

# each row is a sample
umap.dat <- umap(t(FPKM.dat.filt))
umap.dat.plot <- as.data.frame(umap.dat$layout)
# add metadata
umap.dat.plot$PAM50 <- anno$PAM50[match(rownames(umap.dat.plot),anno$METABRIC_ID)]

# plot
plot(umap.dat.plot$V1, umap.dat.plot$V2, 
     col = color.palette[factor(umap.dat.plot$PAM50, levels = names(color.palette))],
     pch = 16,
     main = paste0("FPKM UMAP top most varying genes n=",nrow(FPKM.dat.filt)), 
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
