# Script: other things 
# Author: Lennart Hohmann
# Date: 09.01.2024
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

################################################################################
# Exporting basal ids; Ids with reference year
################################################################################
# input paths
infile.1 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
outfile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_sampleIDs.RData" 
outfile.2 <- "./data/SCANB/0_GroupSamples/ERpHER2n_referenceYear.RData" 

all.samples <- loadRData(infile.1) %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>%
  filter(NCN.PAM50 %in% c("LumA","LumB","Her2","Basal"))

all.samples <- all.samples[c("Sample","NCN.PAM50","ReferenceYear")]
all.samples$ReferenceYear_2010_to_2012 <- ifelse(
  all.samples$ReferenceYear %in% c(2010,2011,2012), 1, 0)

#table(all.samples$NCN.PAM50,all.samples$ReferenceYear_2010_to_2012)

ERpHER2nBasal_sampleIDs <- all.samples[all.samples$NCN.PAM50=="Basal",]$Sample
save(ERpHER2nBasal_sampleIDs, file = outfile.1)

ERpHER2n_referenceYear <- all.samples
#View(ERpHER2n_referenceYear)
save(ERpHER2n_referenceYear, file = outfile.2)

################################################################################


kat <- loadRData("../Project_HER2E/data/SCANB/4_CN/processed/CN_kat_all_samples.RData")
df <- as.data.frame(names(kat))
df$PAM50 <- all.samples[match(df$`names(kat)`,all.samples$Sample),]$NCN.PAM50
View(df)
table(is.na(df$PAM50))
