# Script: Mutational & Rearrangement signatures in SCAN-B
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

wgs.sum <- read_excel(infile.2, sheet = "Summary")
sign.rearr <- read_excel(infile.2, sheet = "RearrangmentSigs")
sign.mut <- read_excel(infile.2, sheet = "SBSsigs")
wgs.hrd <- read_excel(infile.2, sheet = "HRDetect")

basis.anno <- loadRData(infile.3)
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),]

