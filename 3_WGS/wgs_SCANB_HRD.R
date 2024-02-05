# Script: HRD in SCAN-B
# Author: Lennart Hohmann
# Date: 31.01.2024
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
pacman::p_load(readxl)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/3_WGS/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.4 <- "./data/Parameters/color_palette.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_waterfall_ERpHER2nBasal.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# check which samples passed QC
wgs.sum <- read_excel(infile.2, sheet = "Summary")
#View(wgs.sum)
sign.rearr <- read_excel(infile.2, sheet = "RearrangmentSigs")
sign.mut <- read_excel(infile.2, sheet = "SBSsigs")

# convert SCANB to proportion per sample
mut.scanb <- as.data.frame(t(apply(t(mut.scanb),2,function(x) (x/sum(x, na.rm=TRUE))))) 

basis.anno <- loadRData()
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & basis.anno$PAM50_AIMS %in% c("LumA","LumB"),]
mut.basis$MutSigProportions <- as.data.frame(mut.basis$MutSigProportions)

# check which signatures are present in both cohorts
common.mut.sigs <- intersect(colnames(mut.basis$MutSigProportions),colnames(mut.scanb))

# exclude non-shared ones
mut.scanb <- mut.scanb[,which(names(mut.scanb) %in% common.mut.sigs)]
mut.basis$MutSigProportions <- mut.basis$MutSigProportions[,which(
  names(mut.basis$MutSigProportions) %in% common.mut.sigs)]
# also check mapping file from johan
mut.basis <- as.data.frame(do.call(cbind, mut.basis))
names(mut.basis) <- gsub("MutSigProportions.","", names(mut.basis))

# put together in an object ready for plotting 
mut.data <- rbind(mut.basis,mut.scanb) # add subtype


#######################################################################
# plot
#######################################################################


#######################################################################
# stats
#######################################################################



#######################################################################
#######################################################################
