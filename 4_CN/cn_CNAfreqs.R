# Script: Create summarized CNA files for SCANB and Calc. CNA Frequencies in SCAN-B
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

#infile.3 <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv" # move ot basal project
# output paths
outfile.1 <- paste0(data.path,"CNA_Basal.RData")
outfile.2 <- paste0(data.path,"CNA_Basal_GLFreqs.RData")
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load Basal ids
basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load ASCAT gene data
ascat.list <- loadRData(infile.2)
names(ascat.list) <- gsub("\\..*", "", names(ascat.list))
ascat.list <- ascat.list[names(ascat.list) %in% basal.ids]

#######################################################################
# output: GainFreq and LossFreq for each gene
#######################################################################

# get data format: gene sample1 sample2 ... 
ascat.df <- do.call(rbind, lapply(ascat.list, function(x) x$CNA))
ascat.df <- t(ascat.df)
ascat.df <- cbind(ascat.list[[1]][c("gene","chr","start","end")],ascat.df)
#View(ascat.df)
save(ascat.df,file=outfile.1)

# get data format: Gene # GainFreq # LossFreq
# calc loss/gain freqs per group
ascat.df$freqloss.Basal <- apply(
  ascat.df[,5:ncol(ascat.df)], 1, function(x) (
    length(which(x==-1))/ncol(ascat.df[,5:ncol(ascat.df)]))*-100) # i add a minus to make it easier for plotting

ascat.df$freqgain.Basal <- apply(
  ascat.df[,5:ncol(ascat.df)], 1, function(x) (
    length(which(x==1))/ncol(ascat.df[,5:ncol(ascat.df)]))*100)

gene.CNA.freqs <- ascat.df[c("gene","chr","start","end","freqloss.Basal","freqgain.Basal")]

save(gene.CNA.freqs,file=outfile.2)

# add the basis data to this for later plotting
# only keep common genes

x <- loadRData("./data/BASIS/4_CN/raw/LumA_CollectedFrequencyData.RData")
View(x)
