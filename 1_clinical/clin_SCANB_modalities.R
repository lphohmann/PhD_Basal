# Script: overview data modalities SCAN-B
# Author: Lennart Hohmann
# Date: 02.03.2025
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_PredictRecurrence/")
#-------------------
# packages
#source("./scripts/src/")
#if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
#-------------------
# set/create output directories
#output.path <- "./output/"
#dir.create(output.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
outfile.1 <- paste0(output.path,"_modalitiesVenn.pdf") #.txt
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()
#######################################################################

# function: loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path){
  load(file.path)
  get(ls()[ls() != "file.path"])
} 

#######################################################################
# clinical data
#######################################################################

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]

#View(anno)
anno <- anno[c("Sample","ER","PR","HER2","LN.spec",
     "NHG","Size.mm","TreatGroup","DRFi_days",
     "Age",  "OSbin","OS","RFIbin","RFI",
     "DRFIbin","DRFI","NCN.PAM50")]               

#######################################################################
# RNA-seq data
#######################################################################

#######################################################################
# DNA methylation data (n=499)
#######################################################################

#######################################################################
# WGS variant calling data (n=?) 499 -> is published?
#######################################################################
