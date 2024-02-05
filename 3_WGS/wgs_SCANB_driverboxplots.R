# Script: Driver mutation boxplots in SCAN-B/BASIS
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
pacman::p_load(readxl,
               reshape2)
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
infile.3 <- "./data/SCANB/3_WGS/raw/Project2_Basal_like_Drivers_22Jan24_ForJohan.xlsx"
infile.4 <- "./data/BASIS/3_WGS/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx"
infile.5 <- "data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"
infile.6 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_boxplot_ERpHER2nBasal.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- unname(unlist(loadRData(infile.1)[c("ERpHER2n_Basal")]))

# driver data
driv.indel <- read_excel(infile.3, sheet = "PindelDrivers")
driv.point <- read_excel(infile.3, sheet = "CavemanDrivers")
driv.rearr <- read_excel(infile.3, sheet = "BRASS_drivers")
#driv.amp <- NULL # need ID key, so ask Johan

wgs.sum <- read_excel(infile.2, sheet = "Summary")
sign.rearr <- read_excel(infile.2, sheet = "RearrangmentSigs")
sign.mut <- read_excel(infile.2, sheet = "SBSsigs")
wgs.hrd <- read_excel(infile.2, sheet = "HRDetect")

# load IDkey and correct sampleIDs -> ask Johan for key 




# driver genes from BASIS LumA/LumB
driv.basis <- as.data.frame(read_excel(infile.4, 
                                       sheet = "COMBINED_EVENTS"))
View(driv.basis)
basis.anno <- loadRData(infile.6)
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & basis.anno$PAM50_AIMS %in% c("LumA","LumB"),]
table(basis.anno$PAM50_AIMS)


