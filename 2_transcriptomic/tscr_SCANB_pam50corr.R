# Script: PAM50 correlations in SCAN-B
# Author: Lennart Hohmann
# Date: 26.01.2024
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
#pacman::p_load()
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
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/pam50ann_correlations_summarizedByMean_REL4.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_pam50corr.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- unlist(loadRData(infile.1)[c("ERpHER2n_Basal")])

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% sampleIDs, ]

# load data
corr.data <- loadRData(infile.3)
corr.data <- corr.data[corr.data$Assay %in% anno$GEX.assay,]
corr.data <- corr.data[!duplicated(corr.data$Assay),]

#######################################################################
# Plot
#######################################################################
# what to plot now

# LumA vs LumB
#plot(x=corr.data$meanLumA,y=corr.data$meanLumB, pch = 16)
#boxplot(corr.data[c("meanBasal","meanHer2","meanLumA","meanLumB","meanNormal")])

#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
bp <- boxplot(corr.data[c("meanBasal","meanHer2","meanLumA","meanLumB","meanNormal")],
        ylab="Centroid correlation values",
        main="PAM50 centroid correlations in Basal-like cases")
axis(3,at=1:length(bp$n),labels=bp$n)

dev.off()