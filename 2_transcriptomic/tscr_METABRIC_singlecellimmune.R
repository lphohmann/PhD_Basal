# Script: Danenberg immune part in metabric
# Author: Lennart Hohmann
# Date: 30.05.2024
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
# for plots
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
infile.2 <- "./data/METABRIC/Danenberg/MBTMEIMCPublic/SingleCells.csv"

# output paths
outfile.1 <- paste0(data.path,"SingleCells_PTcounts.RData")# total counts of cell phenotypes per sample
outfile.2 <- paste0(data.path,"SingleCells_PTprop.RData")# ratio of cell phenotypes per sample (based on total cells)
outfile.3 <- paste0(data.path,"SingleCells_PTmedians.RData")# nested list with one df for each cell phenotype, having one row per sample with the mean and meadian of each of the markers (calculated over all the cells within the sample)
#plot.file <- paste0(output.path,cohort,"_singlecellimmune.pdf")
#txt.file <- paste0(output.path,cohort,"_singlecellimmune.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

#anno <- loadRData(infile.1)
#anno$Chemotherapy[is.na(anno$Chemotherapy)] <- 0
#anno$Endocrine[is.na(anno$Endocrine)] <- 0
#erp.dat <- anno[anno$HER2_amp == "no" & anno$ER_IHC_status =="pos",]
#tnbc.dat <- anno[anno$HER2_amp == "no" & anno$ER_IHC_status =="neg",]

# danenberg 
cell.dat <- read.table(infile.2,sep=",",header=TRUE)
#table(cell.dat$cellPhenotype)

# # check how many case numbers
# her2 <- cell.dat[cell.dat$metabric_id %in% 
#                    erp.dat$METABRIC_ID[erp.dat$PAM50 =="Her2"],]
# basal <- cell.dat[cell.dat$metabric_id %in% 
#                     erp.dat$METABRIC_ID[erp.dat$PAM50 =="Basal"],]

#######################################################################
# only include 1 image per sample
#######################################################################

# only include 1 image per sample (select first one, or one with most cells?)
#x <- cell.dat[c("ImageNumber","metabric_id")]
#x <- x[!duplicated(x), ]
#table(table(x$metabric_id))
# approach 1: keep first image for each 
# find first occurence of each id and the corresponding image number
first_occurrences <- !duplicated(cell.dat$metabric_id)
# keep only the first occurrences
first_occurrences <- cell.dat[first_occurrences, c("metabric_id","ImageNumber")]
cell.dat.filt <- cell.dat[cell.dat$ImageNumber %in% first_occurrences$ImageNumber,]
#View(cell.dat.filt)
#######################################################################
# cell phenotype proportions
#######################################################################
#ratio of cell phenotypes per sample (based on total cells)
count.df <- data.frame(Phenotype = unique(cell.dat.filt$cellPhenotype))

# df with one row per sample with total counts of cell phenotypes

for (sample.id in unique(cell.dat.filt$metabric_id)) {
  print(sample.id)
  #sample.id = unique(cell.dat.filt$metabric_id)[1]
  s.dat <- cell.dat.filt[cell.dat.filt$metabric_id == sample.id,]
  # cell counts
  tbl <- as.data.frame(table(s.dat$cellPhenotype))
  count.df[[sample.id]] <- tbl$Freq[match(count.df$Phenotype,tbl$Var1)]
}

save(count.df, file=outfile.1)

#######################################################################
# df with one row per sample with ratio of cell phenotypes
#######################################################################

proportion.df <- count.df[1]
# Divide each value by the sum of its respective column
proportion.df[2:ncol(count.df)] <- apply(count.df[2:ncol(count.df)], 2, 
                       function(x) x / sum(x,na.rm = TRUE))
names(proportion.df) <- names(count.df)

save(proportion.df, file=outfile.2)

#######################################################################
# nested list with one df for each cell phenotype, having one row per sample with the mean and median of each of the markers (calculated over all the cells within the sample)
#######################################################################
res.list <- list()
for (pt in unique(cell.dat.filt$cellPhenotype)) {
  #pt = unique(cell.dat.filt$cellPhenotype)[1]
  print(pt)
  pt.dat <- cell.dat.filt[cell.dat.filt$cellPhenotype==pt,]
  pt.list <- split(pt.dat,pt.dat$metabric_id)
  
  medians <- lapply(pt.list,function(x) {
    sample.id <- x$metabric_id[1]
    print(sample.id)
    sample.medians <- apply(x[5:ncol(x)], 2, median)
    label <- c("metabric_id","cellPhenotype",
               paste0(names(sample.medians),".median"))
    y <- c(sample.id,pt,as.vector(sample.medians))
    names(y) <- label
    return(y)
  })
  
  medians.df <- do.call(rbind, medians)
  res.list <- append(res.list, list(medians.df))
}
names(res.list) <- unique(cell.dat.filt$cellPhenotype)

save(res.list, file=outfile.3)

