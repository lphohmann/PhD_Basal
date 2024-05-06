# Script: Differential gene expression analysis in SCAN-B based on FPKM data in TNBC
# Author: Lennart Hohmann
# Date: 06.05.2024
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
pacman::p_load(biomaRt)
#BiocManager::install("biomaRt")
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.2 <- "./data/SCANB/2_transcriptomic/processed/All_LogScaled_gex.RData"
# output paths
outfile.1 <- paste0(data.path,"DE_result_FPKM_TNBC.RData")
#plot.file <- paste0(output.path,cohort,"_DE_TNBC.pdf")
#txt.file <- paste0(output.path,cohort,"_DE_TNBC.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")]

# load gex
gex.df <- loadRData(infile.2)
gex.df <- gex.df[, unname(unlist(sampleIDs))]

#######################################################################
# 4. Test each gene
#######################################################################

# with apply
de.res <- apply(gex.df,1,function(x) {
  basal.dat <- as.numeric(x[unname(unlist(sampleIDs["ERpHER2n_Basal"]))])
  tnbc.basal.dat <- as.numeric(x[unname(unlist(sampleIDs["TNBC_Basal"]))])
  tnbc.nonbasal.dat <- as.numeric(x[unname(unlist(sampleIDs["TNBC_NonBasal"]))])
  Basal.tnbc.basal.pval <- wilcox.test(basal.dat, tnbc.basal.dat)$p.value
  Basal.tnbc.basal.logFC <- mean(basal.dat)-mean(tnbc.basal.dat)
  Basal.tnbc.nonbasal.pval <- wilcox.test(basal.dat, tnbc.nonbasal.dat)$p.value
  Basal.tnbc.nonbasal.logFC <- mean(basal.dat)-mean(tnbc.nonbasal.dat)
  return(c(Basal.tnbc.basal.pval,Basal.tnbc.basal.logFC,Basal.tnbc.nonbasal.pval,Basal.tnbc.nonbasal.logFC))
})

row.names(de.res) <- c("Basal.tnbc.basal.pval","Basal.tnbc.basal.logFC","Basal.tnbc.nonbasal.pval","Basal.tnbc.nonbasal.logFC")
de.res <- as.data.frame(t(de.res))

# adjust p value 
de.res$Basal.tnbc.basal.padj <- p.adjust(de.res$Basal.tnbc.basal.pval, "fdr")
de.res$Basal.tnbc.nonbasal.padj <- p.adjust(de.res$Basal.tnbc.nonbasal.pval, "fdr")

#View(de.res)
a.sig <- de.res[de.res$Basal.tnbc.basal.padj<=0.05,]
a.sig <- a.sig[abs(a.sig$Basal.tnbc.basal.logFC)>=1,]

b.sig <- de.res[de.res$Basal.tnbc.nonbasal.padj<=0.05,]
b.sig <- a.sig[abs(a.sig$Basal.tnbc.nonbasal.logFC)>=1,]

#length(de.res$Basal.LumA.padj<=0.05)

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
# x <- loadRData("./data/SCANB/2_transcriptomic/processed/DE_result_FPKM.RData")
# View(x[c("ESR1","AR","PGR","STAT3","STAT6","ABAT","AGR2","CA12","DNALI1","FOXA1","GATA3","SLC44A4","TBC1D9","XBP1"),])
# "ESR1","AR","PGR","STAT3","STAT6","ABAT","AGR2","CA12","DNALI1","FOXA1","GATA3","SLC44A4","TBC1D9","XBP1"

