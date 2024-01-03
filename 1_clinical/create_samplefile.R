# Script: sampleID included in analyses for SCANB
# Author: Lennart Hohmann
# Date: 01.01.2024
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
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, openxlsx)
#-------------------
# input paths
infile.1 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
outfile.ERpHER2n_sampleIDs <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData" 

################################################################################
# ERpHER2n
################################################################################
all.samples <- loadRData(infile.1) %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>% 
  filter(NCN.PAM50 %in% c("LumA","LumB","Her2","Basal"))

ERpHER2n_sampleIDs <- list(
  "ERpHER2n_Basal"= all.samples[which(all.samples$NCN.PAM50=="Basal"),]$Sample,
  "ERpHER2n_HER2E"= all.samples[which(all.samples$NCN.PAM50=="Her2"),]$Sample,
  "ERpHER2n_LumA"= all.samples[which(all.samples$NCN.PAM50=="LumA"),]$Sample,
  "ERpHER2n_LumB"= all.samples[which(all.samples$NCN.PAM50=="LumB"),]$Sample)

save(ERpHER2n_sampleIDs, file = outfile.ERpHER2n_sampleIDs)

# define color palette for plotting also
color.palette <- c("Her2"="#d334eb","LumA"="#2176d5","LumB"="#34c6eb","Basal"="#ba0606")
save(color.palette,file="./data/Parameters/color_palette.RData")
################################################################################