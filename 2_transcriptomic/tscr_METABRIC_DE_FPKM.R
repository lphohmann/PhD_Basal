# Script: Differential gene expression analysis in METABRIC based on FPKM data
# Author: Lennart Hohmann
# Date: 11.11.2024
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
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
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
infile.1 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/METABRIC/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
# output paths
outfile.1 <- paste0(data.path,"DE_result_FPKM.RData")
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

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]
#lapply(sampleIDs, length)

# load gex
gex.df <- loadRData(infile.2)
gex.df <- gex.df[, colnames(gex.df) %in% unname(unlist(sampleIDs))]

# Filter sampleIDs based on IDs in gex data
sampleIDs <- lapply(sampleIDs, function(ids) ids[ids %in% colnames(gex.df)])
#lapply(sampleIDs, length)

#######################################################################
# 4. Test each gene
#######################################################################

# with apply
de.res <- apply(gex.df,1,function(x) {
  basal.dat <- as.numeric(x[unname(unlist(sampleIDs["ERpHER2n_Basal"]))])
  luma.dat <- as.numeric(x[unname(unlist(sampleIDs["ERpHER2n_LumA"]))])
  lumb.dat <- as.numeric(x[unname(unlist(sampleIDs["ERpHER2n_LumB"]))])
  Basal.LumA.pval <- wilcox.test(basal.dat, luma.dat)$p.value
  Basal.LumA.logFC <- mean(basal.dat)-mean(luma.dat)
  Basal.LumB.pval <- wilcox.test(basal.dat, lumb.dat)$p.value
  Basal.LumB.logFC <- mean(basal.dat)-mean(lumb.dat)
  return(c(Basal.LumA.pval,Basal.LumA.logFC,Basal.LumB.pval,Basal.LumB.logFC))
})

row.names(de.res) <- c("Basal.LumA.pval","Basal.LumA.logFC","Basal.LumB.pval","Basal.LumB.logFC")
de.res <- as.data.frame(t(de.res))

# adjust p value 
de.res$Basal.LumA.padj <- p.adjust(de.res$Basal.LumA.pval, "fdr")
de.res$Basal.LumB.padj <- p.adjust(de.res$Basal.LumB.pval, "fdr")

#View(de.res)
#a.sig <- de.res[de.res$Basal.LumA.padj<=0.05,]
#a.sig <- a.sig[abs(a.sig$Basal.LumA.logFC)>=1,]

#sum(de.res$Basal.LumA.padj<=0.05)

#######################################################################
# 4. save
#######################################################################

# save the file
save(de.res,file = outfile.1)
# 
# 
# 
# 
# ## check
#x <- loadRData("./data/METABRIC/2_transcriptomic/processed/DE_result_FPKM.RData")
#View(x[c("ESR1","AR","PGR","STAT3","STAT6","ABAT","AGR2","CA12","DNALI1","FOXA1","GATA3","SLC44A4","TBC1D9","XBP1"),])
# "ESR1","AR","PGR","STAT3","STAT6","ABAT","AGR2","CA12","DNALI1","FOXA1","GATA3","SLC44A4","TBC1D9","XBP1"

