# Script: sampleID included in analyses for SCANB and METABRIC
# Author: Lennart Hohmann
# Date: 01.01.2024
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
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, openxlsx)
#-------------------
# input paths
infile.1 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.2 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.3 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"

# for johan
infile.4 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGS_sampleIDs.RData"

# output paths
outfile.color.palette <- "./data/Parameters/color_palette.RData"
outfile.color.palette_TNBC <- "./data/Parameters/TNBC_color_palette.RData"
outfile.scanb.ERpHER2n_sampleIDs <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData" 
outfile.scanb.TNBC_sampleIDs <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData" 

outfile.metabric.ERpHER2n_sampleIDs <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
outfile.metabric.TNBC_sampleIDs <- "./data/METABRIC/0_GroupSamples/TNBC_sampleIDs.RData" 

################################################################################
# color palettes
################################################################################

# define color palette for plotting also
color.palette <- c("LumA"="#2176d5","LumB"="#34c6eb","Basal"="#ba0606") #,"Her2"="#d334eb"
save(color.palette, file = outfile.color.palette)

# define color palette for plotting also
color.palette_TNBC <- c("TNBC_NonBasal"="#fec44f","TNBC_Basal"="#FF2400","ERpHER2n_Basal"="#ba0606")
save(color.palette_TNBC, file = outfile.color.palette_TNBC)

################################################################################
# SCAN-B
################################################################################

################################################################################
# ERpHER2n

all.samples <- loadRData(infile.1) %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(ER=="Positive") %>% 
  filter(HER2=="Negative") %>% 
  filter(NCN.PAM50 %in% c("LumA","LumB","Basal"))

# check treatments
#View(all.samples[which(all.samples$NCN.PAM50=="Her2"),])

ERpHER2n_sampleIDs <- list(
  "ERpHER2n_Basal"= all.samples[which(all.samples$NCN.PAM50=="Basal"),]$Sample,
  #"ERpHER2n_HER2E"= all.samples[which(all.samples$NCN.PAM50=="Her2"),]$Sample,
  "ERpHER2n_LumA"= all.samples[which(all.samples$NCN.PAM50=="LumA"),]$Sample,
  "ERpHER2n_LumB"= all.samples[which(all.samples$NCN.PAM50=="LumB"),]$Sample)

save(ERpHER2n_sampleIDs, file = outfile.scanb.ERpHER2n_sampleIDs)

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

save(TNBC_sampleIDs, file = outfile.scanb.TNBC_sampleIDs)

################################################################################
# for johan
infile.4 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGS_sampleIDs.RData"

wgs.basal.scanb <- loadRData(infile.4)
all.basal.scanb <- ERpHER2n_sampleIDs[["ERpHER2n_Basal"]]

df.basal <- data.frame(
  SampleID = all.basal.scanb,
  PAM50 = "Basal",
  WGS = ifelse(all.basal.scanb %in% wgs.basal.scanb$Lund.tumour.id, 1, 0),
  TumourID = ifelse(all.basal.scanb %in% wgs.basal.scanb$Lund.tumour.id, 
                    wgs.basal.scanb$Tumour[match(all.basal.scanb, 
                                                 wgs.basal.scanb$Lund.tumour.id)], NA)
)

save(df.basal,file="./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGSanno.RData")

################################################################################
# METABRIC
################################################################################

################################################################################
# ERpHER2n
metabric.all <- loadRData(infile.3)
metabric.all <- metabric.all[!is.na(metabric.all$METABRIC_ID),]
#metabric.all$Chemotherapy[is.na(metabric.all$Chemotherapy)] <- 0
#metabric.all$Endocrine[is.na(metabric.all$Endocrine)] <- 0
erp.dat <- metabric.all[which(metabric.all$HER2_amp == "no" & 
                                metabric.all$ER_IHC_status =="pos"),]
#table(is.na(erp.dat$ClinGroup))
table(erp.dat$PAM50)
ERpHER2n_sampleIDs.mb <- list(
  "ERpHER2n_Basal"= erp.dat[which(erp.dat$PAM50=="Basal"),]$METABRIC_ID,
  #"ERpHER2n_HER2E"= erp.dat[which(erp.dat$PAM50=="Her2"),]$METABRIC_ID,
  "ERpHER2n_LumA"= erp.dat[which(erp.dat$PAM50=="LumA"),]$METABRIC_ID,
  "ERpHER2n_LumB"= erp.dat[which(erp.dat$PAM50=="LumB"),]$METABRIC_ID)

save(ERpHER2n_sampleIDs.mb, file = outfile.metabric.ERpHER2n_sampleIDs)

################################################################################
# TNBC data

tnbc.dat <- metabric.all[which(metabric.all$ClinGroup=="TNBC"),] # but NA for some samples

# define groups
table(tnbc.dat$PAM50)
tnbc.dat$TNBC_Group <- ifelse(tnbc.dat$PAM50 == "Basal", "TNBC_Basal", "TNBC_NonBasal")
table(tnbc.dat$TNBC_Group)

TNBC_sampleIDs.metabric <- list(
  "ERpHER2n_Basal"= ERpHER2n_sampleIDs.mb[["ERpHER2n_Basal"]],
  "TNBC_Basal"= tnbc.dat[which(tnbc.dat$TNBC_Group=="TNBC_Basal"),]$METABRIC_ID,
  "TNBC_NonBasal"= tnbc.dat[which(tnbc.dat$TNBC_Group=="TNBC_NonBasal"),]$METABRIC_ID)

save(TNBC_sampleIDs.metabric, file = outfile.metabric.TNBC_sampleIDs)


