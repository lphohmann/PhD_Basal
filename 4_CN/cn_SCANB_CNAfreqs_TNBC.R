# Script: IN TNBC: Create summarized CNA files for SCANB and Calc. CNA Frequencies 
# Author: Lennart Hohmann
# Date: 22.03.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
#cohort <- "SCANB"
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/processed/ASCAT_genelevel_TNBC.RData"
infile.9 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.10 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
# output paths
outfile.1 <- paste0(data.path,"CNA_genelevel_TNBC.RData")
outfile.2 <- paste0(data.path,"CNA_GLFreqs_TNBC.RData")
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
sample.ids <- loadRData(infile.1)
names(sample.ids)
#"ERpHER2n_Basal" "TNBC_Basal" "TNBC_NonBasal" 

# load ASCAT gene data
ascat.list.scanb <- loadRData(infile.2)
# map PD ids to S ids
id.key <- loadRData(infile.10)[c("PD_ID","External_ID_sample")]
names(ascat.list.scanb) <- id.key$External_ID_sample[match(names(ascat.list.scanb),id.key$PD_ID)]
ascat.list.tnbc <- ascat.list.scanb[names(ascat.list.scanb) %in% unname(unlist(sample.ids[c("TNBC_Basal","TNBC_NonBasal")]))]

# get data format: gene sample1 sample2 ... 
ascat.df.tnbc <- do.call(rbind, lapply(ascat.list.tnbc, function(x) x$CNA))
ascat.df.tnbc <- t(ascat.df.tnbc)
ascat.df.tnbc <- cbind(ascat.list.tnbc[[1]][c("gene","chr","start","end")],ascat.df.tnbc)
#View(ascat.df.tnbc)

# get freqs
CNA.freqs.tnbc.basal <- ascat.df.tnbc[colnames(ascat.df.tnbc) %in% 
                                        c(c("gene","chr","start","end"),unname(unlist(sample.ids[c("TNBC_Basal")])))]
CNA.freqs.tnbc.nonbasal <- ascat.df.tnbc[colnames(ascat.df.tnbc) %in% 
                                           c(c("gene","chr","start","end"),unname(unlist(sample.ids[c("TNBC_NonBasal")])))]


# calc loss/gain freqs per group
CNA.freqs.tnbc.basal$freqloss.tnbc.Basal <- apply(
  CNA.freqs.tnbc.basal[,5:ncol(CNA.freqs.tnbc.basal)], 1, function(x) (
    length(which(x==-1))/ncol(CNA.freqs.tnbc.basal[,5:ncol(CNA.freqs.tnbc.basal)]))*-100)
CNA.freqs.tnbc.basal$freqgain.tnbc.Basal <- apply(
  CNA.freqs.tnbc.basal[,5:ncol(CNA.freqs.tnbc.basal)], 1, function(x) (
    length(which(x==1))/ncol(CNA.freqs.tnbc.basal[,5:ncol(CNA.freqs.tnbc.basal)]))*100)
CNA.freqs.tnbc.nonbasal$freqloss.tnbc.NonBasal <- apply(
  CNA.freqs.tnbc.nonbasal[,5:ncol(CNA.freqs.tnbc.nonbasal)], 1, function(x) (
    length(which(x==-1))/ncol(CNA.freqs.tnbc.nonbasal[,5:ncol(CNA.freqs.tnbc.nonbasal)]))*-100)
CNA.freqs.tnbc.nonbasal$freqgain.tnbc.NonBasal <- apply(
  CNA.freqs.tnbc.nonbasal[,5:ncol(CNA.freqs.tnbc.nonbasal)], 1, function(x) (
    length(which(x==1))/ncol(CNA.freqs.tnbc.nonbasal[,5:ncol(CNA.freqs.tnbc.nonbasal)]))*100)

CNA.freqs.tnbc <- cbind(CNA.freqs.tnbc.basal,
                        CNA.freqs.tnbc.nonbasal[,5:ncol(CNA.freqs.tnbc.nonbasal)])

CNA.freqs.tnbc <- CNA.freqs.tnbc[c("gene","chr","start","end",
                                   "freqloss.tnbc.Basal","freqgain.tnbc.Basal",
                                   "freqloss.tnbc.NonBasal","freqgain.tnbc.NonBasal")]

#######################################################################
# prepare and save files 
#######################################################################

# exclude rows with NA
ascat.df.tnbc <- ascat.df.tnbc[complete.cases(ascat.df.tnbc), ]

# only include common genes in frequencies
CNA.freqs.tnbc <- CNA.freqs.tnbc[CNA.freqs.tnbc$gene %in% ascat.df.tnbc$gene, ]

# save
subtype.samples <- list("TNBC.Basal"=unname(unlist(sample.ids["TNBC_Basal"])),
                        "TNBC.NonBasal"=unname(unlist(sample.ids["TNBC_NonBasal"])))
# make into one df
ascat.save.file <- list("ascat.df.tnbc" = ascat.df.tnbc,
                        "subtype.samples" = subtype.samples)
save(ascat.save.file,file=outfile.1)
save(CNA.freqs.tnbc,file=outfile.2)
