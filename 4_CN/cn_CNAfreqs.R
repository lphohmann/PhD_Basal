# Script: Create summarized CNA files for SCANB/BASIS and Calc. CNA Frequencies in SCAN-B
# Author: Lennart Hohmann
# Date: 11.02.2024
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
#pacman::p_load()
#-------------------
# set/create output directories
# for plots
output.path <- "./output/4_CN/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"
infile.3 <- "./data/BASIS/4_CN/raw/CNA_genelevel_all.RData"
infile.4 <- "./data/BASIS/4_CN/raw/LumA_CollectedFrequencyData.RData"
infile.5 <- "./data/BASIS/4_CN/raw/LumB_CollectedFrequencyData.RData"
infile.6 <- "./data/BASIS/4_CN/processed/ProbeGeneMap.RData.RData"

#infile.3 <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv" # move ot basal project
# output paths
outfile.1 <- paste0(data.path,"CNA_genelevel_all.RData")
outfile.2 <- paste0(data.path,"CNA_GLFreqs_all.RData")
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load SCANB data
#######################################################################

# load Basal ids
basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load ASCAT gene data
ascat.list <- loadRData(infile.2)
names(ascat.list) <- gsub("\\..*", "", names(ascat.list))
ascat.list <- ascat.list[names(ascat.list) %in% basal.ids]

# get data format: gene sample1 sample2 ... 
ascat.df.scanb <- do.call(rbind, lapply(ascat.list, function(x) x$CNA))
ascat.df.scanb <- t(ascat.df.scanb)
ascat.df.scanb <- cbind(ascat.list[[1]][c("gene","chr","start","end")],ascat.df.scanb)
#View(ascat.df.scanb)

# get freqs
gene.CNA.freqs <- ascat.df.scanb
# calc loss/gain freqs per group
gene.CNA.freqs$freqloss.Basal <- apply(
  gene.CNA.freqs[,5:ncol(gene.CNA.freqs)], 1, function(x) (
    length(which(x==-1))/ncol(gene.CNA.freqs[,5:ncol(gene.CNA.freqs)]))*-100) # i add a minus to make it easier for plotting

gene.CNA.freqs$freqgain.Basal <- apply(
  gene.CNA.freqs[,5:ncol(gene.CNA.freqs)], 1, function(x) (
    length(which(x==1))/ncol(gene.CNA.freqs[,5:ncol(gene.CNA.freqs)]))*100)

gene.CNA.freqs <- gene.CNA.freqs[c("gene","chr","start","end","freqloss.Basal","freqgain.Basal")]

#######################################################################
# load BASIS data and save files
#######################################################################

luma.ids
lumb.ids

# load ASCAT gene data
#ascat.list <- loadRData(infile.2)
#names(ascat.list) <- gsub("\\..*", "", names(ascat.list))
#ascat.list <- ascat.list[names(ascat.list) %in% basal.ids]

# get data format: gene sample1 sample2 ... 


#######################################################################
# load BASIS data and save files
#######################################################################

ascat.df.basis <- loadRData(infile.3)
common.genes <- intersect(ascat.df.scanb$gene,ascat.df.basis$gene)
ascat.df.scanb <- ascat.df.scanb[ascat.df.scanb$gene %in% common.genes,]
ascat.df.basis <- ascat.df.basis[ascat.df.basis$gene %in% common.genes,]
ascat.df.basis <- ascat.df.basis[, !names(ascat.df.basis) %in% c("chr", "centerPos", "Genome_pos")]
ascat.df.all <- merge(ascat.df.scanb,ascat.df.basis,by="gene")
save(ascat.df.all,file=outfile.1)

# freq data
cn.luma <- loadRData(infile.4)
cn.luma <- do.call("cbind", list(cn.luma$fData,"LumA_Gain"=cn.luma$CN_Gain,"LumA_Loss"=cn.luma$CN_Loss))
cn.lumb <- loadRData(infile.5)
cn.lumb <- do.call("cbind", list(cn.lumb$fData,"LumB_Gain"=cn.lumb$CN_Gain,"LumB_Loss"=cn.lumb$CN_Loss))
cn.basis <- merge(cn.luma,cn.lumb,by=c("reporterId","chromosome","centerPosition"))
basis.map <- loadRData(infile.6)
cn.basis.mapped <- merge(cn.basis,basis.map,by="reporterId") 

# filter NA rows
cn.basis.mapped <- cn.basis.mapped[which(!is.na(cn.basis.mapped$Gene_symbol)),]

# 
common.genes <- intersect(gene.CNA.freqs$gene,cn.basis.mapped$Gene_symbol)
gene.CNA.freqs <- gene.CNA.freqs[gene.CNA.freqs$gene %in% common.genes,]
cn.basis.mapped <- cn.basis.mapped[cn.basis.mapped$Gene_symbol %in% common.genes,]
cn.basis.mapped <- cn.basis.mapped[!duplicated(cn.basis.mapped$Gene_symbol), ]
names(cn.basis.mapped)[names(cn.basis.mapped) == "Gene_symbol"] <- "gene"

gene.CNA.freqs.all <- merge(gene.CNA.freqs,cn.basis.mapped[c("gene","LumA_Gain",
                                                             "LumA_Loss","LumB_Gain","LumB_Loss")],
                            by="gene")
names(gene.CNA.freqs.all)[names(gene.CNA.freqs.all) == "freqloss.Basal"] <- "Basal_Loss"
names(gene.CNA.freqs.all)[names(gene.CNA.freqs.all) == "freqgain.Basal"] <- "Basal_Gain"


gene.CNA.freqs.all$LumA_Loss <- gene.CNA.freqs.all$LumA_Loss*(-1)
gene.CNA.freqs.all$LumB_Loss <- gene.CNA.freqs.all$LumB_Loss*(-1)

save(gene.CNA.freqs.all,file=outfile.2)

