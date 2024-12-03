# Script: TNBC part - Plotting and testing expression of metagenes in SCAN-B
# Author: Lennart Hohmann
# Date: 08.01.2024
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
pacman::p_load(readxl)
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
#infile.1 <- "./data/SCANB/5_TNBC_NatMed/ASCAT_InSilico_SCAN_B_TNBC_EasySegments.RData"
#infile.2 <- "./data/SCANB/5_TNBC_NatMed/driver.df.scanb.complete.csv"
#infile.3 <- "./data/SCANB/5_TNBC_NatMed/Final_fData_Object.RData"
#infile.4 <- "./data/SCANB/5_TNBC_NatMed/SCANB_TNBC_collected_kataegis.RData"
infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.7 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.8 <- "./data/Parameters/TNBC_color_palette.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_TILcounts_TNBC.pdf")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.7)[c("TNBC_Basal","TNBC_NonBasal" )]

# load palette
color.palette <- loadRData(infile.8)

# load gex
til.dat <- loadRData(infile.5)
til.dat <- til.dat[which(til.dat$External_ID_sample %in% unname(unlist(sampleIDs))), 
                   c("External_ID_sample","sTIL","TILs")]


til.dat$Subtype <- sapply(til.dat$External_ID_sample, function(sampleID) {
  for (sublist_name in names(sampleIDs)) {
    if (sampleID %in% sampleIDs[[sublist_name]]) {
      return(sublist_name)
    }
  }
  return(NA)
})

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))

# subtype data
tnbc.basal.dat <- til.dat[which(til.dat$Subtype=="TNBC_Basal"),]
tnbc.nonbasal.dat <- til.dat[which(til.dat$Subtype=="TNBC_NonBasal"),]

# plot
bp <- boxplot(list(TNBC_Basal=tnbc.basal.dat$TILs,
              TNBC_NonBasal=tnbc.nonbasal.dat$TILs), 
  col = color.palette[c("TNBC_Basal","TNBC_NonBasal")], 
  names = c("TNBC_Basal","TNBC_NonBasal"),
  ylab = "TIL counts",
  main = "TIL counts")
axis(3,at=1:length(bp$n),labels=bp$n)

# plot
bp <- boxplot(list(TNBC_Basal=tnbc.basal.dat$sTIL,
                   TNBC_NonBasal=tnbc.nonbasal.dat$sTIL), 
              col = color.palette[c("TNBC_Basal","TNBC_NonBasal")], 
              names = c("TNBC_Basal","TNBC_NonBasal"),
              ylab = "sTIL counts",
              main = "sTIL counts")
axis(3,at=1:length(bp$n),labels=bp$n)
#######################################################################

par(mfrow = c(1, 1))
dev.off()

