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
infile.2 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
# output paths
outfile.ERpHER2n_sampleIDs <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData" 
outfile.color.palette <- "./data/Parameters/color_palette.RData"
outfile.TNBC_sampleIDs <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData" 
outfile.color.palette_TNBC <- "./data/Parameters/TNBC_color_palette.RData"


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
save(color.palette, file = outfile.color.palette)

################################################################################
################################################################################
# TNBC data

tnbc.ids <- loadRData(infile.2)[c("PD_ID","External_ID_sample")]
#length(unique(tnbc.ids$External_ID_sample)) #235

anno <- loadRData(infile.1)
anno <- anno[anno$Follow.up.cohort == TRUE,]
tnbc.anno <- anno[anno$Sample %in% tnbc.ids$External_ID_sample,]
tnbc.anno <- tnbc.anno[!duplicated(tnbc.anno$Sample),]
tbl <- table(tnbc.anno$NCN.PAM50)

tnbc.anno$TNBC_Group <- ifelse(tnbc.anno$NCN.PAM50 == "Basal", "TNBC_Basal", "TNBC_NonBasal")

TNBC_sampleIDs <- list(
  "ERpHER2n_Basal"= ERpHER2n_sampleIDs[["ERpHER2n_Basal"]],
  "TNBC_Basal"= tnbc.anno[which(tnbc.anno$TNBC_Group=="TNBC_Basal"),]$Sample,
  "TNBC_NonBasal"= tnbc.anno[which(tnbc.anno$TNBC_Group=="TNBC_NonBasal"),]$Sample)

save(TNBC_sampleIDs, file = outfile.TNBC_sampleIDs)


# define color palette for plotting also
color.palette_TNBC <- c("TNBC_NonBasal"="#fec44f","TNBC_Basal"="#FF2400","ERpHER2n_Basal"="#ba0606")
save(color.palette_TNBC, file = outfile.color.palette_TNBC)
