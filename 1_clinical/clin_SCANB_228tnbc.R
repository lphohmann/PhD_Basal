# Script: check sample # methylation cohort tnbc scanb 
# Author: Lennart Hohmann
# Date: 19.02.2025
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

fu.ids <- loadRData("./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData")
fu.ids <- fu.ids[fu.ids$Follow.up.cohort==TRUE,]
tnbc.anno <- loadRData("./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData")
fu.ids$PDid <- tnbc.anno$PD_ID[match(fu.ids$Sample,tnbc.anno$External_ID_sample)]
fu.ids <- fu.ids[!is.na(fu.ids$PDid),]

epi.dat <- read.table("./data/SCANB/5_TNBC_NatMed/nmfClusters_groupVariablesWithAnno_3groupBasal_2groupNonBasal_byAtac_bySd.txt",sep="\t")
epi.dat <- epi.dat[epi.dat$Sample %in% fu.ids$PDid,]
table(epi.dat$NMF_atacDistal)
