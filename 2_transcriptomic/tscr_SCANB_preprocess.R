# Script: Preprocessing data for gene expression analyses in SCAN-B
# includes: mean centering across all ER+ cases, z-transforming, convert to hugo gene names
# Author: Lennart Hohmann
# Date: 04.01.2024
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
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               ggplot2,
               grid,
               gridExtra,
               ggplotify,
               janitor)
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
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData" 
infile.3 <- "./data/SCANB/2_transcriptomic/raw/genematrix_noNeg.Rdata" 

# output paths
outfile.1 <- paste0(data.path,"ERp_LogScaled_gex.RData") # scaled gex
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load annotation data and select 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE & anno$ER == "Positive", ]

# load and prep gene expr. data
gex <- as.data.frame(loadRData(infile.3))
# correct colnames
gex <- gex[anno$GEX.assay]
names(gex) <- anno$Sample[match(colnames(gex),anno$GEX.assay)]
# correct gene names
gex$Ensemble_ID <- gsub("\\..*$", "", rownames(gex))
rownames(gex) <- NULL
hgnc.table <- ensemble_to_hgnc(gex$Ensemble_ID)
gex$Hgnc_ID <- hgnc.table$hgnc_symbol[match(gex$Ensemble_ID,hgnc.table$ensembl_gene_id)]
gex <- gex[!duplicated(gex$Hgnc_ID), ] # keep 1st row of each symbol
gex <- gex[!is.na(gex$Hgnc_ID),]
rownames(gex) <- gex$Hgnc_ID
gex$Hgnc_ID <- NULL
gex$Ensemble_ID <- NULL

#######################################################################
# log-transform and scale
#######################################################################

# log transform FPKM data
gex <- as.data.frame(log2(gex + 1))
# z-transform
gex <- as.data.frame(t(apply(
  gex, 1, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))) # for some rows there may be 0 variance so i have to handle these cases

#######################################################################
# save
#######################################################################

save(gex,file=outfile.1)
