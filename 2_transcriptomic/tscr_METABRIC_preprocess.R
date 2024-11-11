# Script: Preprocessing data for gene expression analyses in METABRIC
# the data is already log-transformed and z-transformed so it just brings it into the same format as SCAN-B
# Author: Lennart Hohmann
# Date: 08.11.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "METABRIC"
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
data.path <- "./data/METABRIC/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.2 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData" 
infile.3 <- "./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt" 
# output paths
outfile.1 <- paste0(data.path,"ERp_LogScaled_gex.RData") # scaled gex
outfile.2 <- paste0(data.path,"All_LogScaled_gex.RData")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# all samples (follow.up.cohort)
#######################################################################

# load annotation data and select 
anno <- loadRData(infile.2)
anno <- anno[!is.na(anno$METABRIC_ID),]

# load and prep gene expr. data
gex <- as.data.frame(read.table(infile.3, sep="\t"))
names(gex) <- gex[1, ]
gex <- gex[-1, ]

gex$Entrez_Gene_Id <- NULL
#gex[gex == ""] <- NA
gex <- gex[!is.na(gex$Hugo_Symbol), ]
gex <- gex[!duplicated(gex$Hugo_Symbol), ]
row.names(gex) <- gex$Hugo_Symbol
gex$Hugo_Symbol <- NULL
gex <- gex[, colnames(gex) %in% anno$METABRIC_ID]
gex[] <- lapply(gex, as.numeric)
#df_with_na <- gex[apply(gex, 1, function(row) any(is.na(row))), ]
#View(df_with_na)
gex <- na.omit(gex)

#######################################################################
# save
#######################################################################

save(gex,file=outfile.2)

#######################################################################
# save erp samples
#######################################################################

# load annotation data and select 
anno <- loadRData(infile.2)
anno <- anno[!is.na(anno$METABRIC_ID) & anno$ER_IHC_status =="pos",]

gex <- gex[,colnames(gex) %in% anno$METABRIC_ID]
save(gex,file=outfile.1)
