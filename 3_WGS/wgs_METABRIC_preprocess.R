# Script: Preprocess METABRIC mut data into correct format
# Author: Lennart Hohmann
# Date: 11.11.2024
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
#pacman::p_load()
#-------------------
# set/create output directories
# for data
data.path <- "./data/METABRIC/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.0 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.1 <- "./data/METABRIC/3_WGS/raw/data_mutations_extended.txt"
infile.2 <- "./data/METABRIC/4_CN/raw/data_CNA.txt"
infile.3 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
# output paths
#outfile.1 <- #paste0(data.path,"....RData")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load annotation data and select 
anno <- loadRData(infile.3)
anno <- anno[!is.na(anno$METABRIC_ID),]

# load and prep gene expr. data
cn.data <- read.delim(infile.2,header = TRUE, sep = "\t", dec = ".")
colnames(cn.data)[2:ncol(cn.data)] <- gsub("\\.", "-", colnames(cn.data)[2:ncol(cn.data)])

cn.data$Entrez_Gene_Id <- NULL
cn.data <- cn.data[!is.na(cn.data$Hugo_Symbol), ]
cn.data <- cn.data[!duplicated(cn.data$Hugo_Symbol), ]
row.names(cn.data) <- cn.data$Hugo_Symbol
cn.data$Hugo_Symbol <- NULL
cn.data <- cn.data[, colnames(cn.data) %in% anno$METABRIC_ID]
cn.data[] <- lapply(cn.data, as.numeric)
cn.data <- na.omit(cn.data)
#View(cn.data)

#sampleIDs.erp <- loadRData(infile.0)
#cn.data.erp <- cn.data[,names(cn.data) %in% c(unname(unlist(sampleIDs.erp))]

x <- loadRData("./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nBasal.RData")
View(x)

x <- loadRData("./data/SCANB/3_WGS/processed/")
View(x)

mut.data <- as.data.frame(read.delim(infile.1, header = FALSE, sep = "\t", dec = "."))
colnames(mut.data) <- mut.data[2,]
mut.data <- mut.data[-c(1, 2), ]
mut.data <- mut.data[c("Hugo_Symbol", 
                       "Variant_Classification",
                       "Tumor_Sample_Barcode")]
mut.data <- mut.data[mut.data$Tumor_Sample_Barcode %in% anno$METABRIC_ID,]

View(cn.data)

