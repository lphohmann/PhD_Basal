# Script: Driver mutation waterfallplot in SCAN-B TNBC groups
# Author: Lennart Hohmann
# Date: 10.05.2024
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
               reshape2,
               GenVisR,
               RColorBrewer,
               data.table)
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.4 <- "./data/SCANB/5_TNBC_NatMed/driver.df.scanb.complete.csv"
infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_waterfall_TNBC.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data for TNBC Basal samples
#######################################################################

# driver genes from tnbc groups
tnbc.samples <- loadRData(infile.1)[c("TNBC_Basal","TNBC_NonBasal")]
driv.tnbc <- read.table(infile.4,sep=",",header=TRUE)
tnbc.anno <- loadRData(infile.5)
#View(tnbc.anno)
tnbc.anno$Subtype <- sapply(tnbc.anno$External_ID_sample, function(sampleID) {
  for (sublist_name in names(tnbc.samples)) {
    if (sampleID %in% tnbc.samples[[sublist_name]]) {
      return(sublist_name)
    }
  }
  return(NA)
})

# convert sample name
driv.tnbc$SampleID <- tnbc.anno$External_ID_sample[match(driv.tnbc$Sample,tnbc.anno$PD_ID)]
driv.tnbc$Subtype <- tnbc.anno$Subtype[match(driv.tnbc$SampleID,tnbc.anno$External_ID_sample)]

tnbc.anno <- tnbc.anno[tnbc.anno$External_ID_sample %in% unlist(tnbc.samples),]
driv.tnbc <- driv.tnbc[driv.tnbc$SampleID %in% unlist(tnbc.samples),]
#View(driv.tnbc)

# get into right format to match scanb basal

driv.tnbc$variant_class <- paste0(driv.tnbc$Type,"_",driv.tnbc$Effect)
# not included in the basal drivers
driv.tnbc <- driv.tnbc[driv.tnbc$variant_class != "loss_NA", ]
# indel framesift, indel inframe
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class %in% c("Del_frameshift","Ins_frameshift"),
                                  "indel_frameshift",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class %in% c("Del_inframe","Ins_inframe"),
                                  "indel_inframe",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "amp_NA",
                                  "CN_amplification",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "Sub_missense",
                                  "point_missense",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "Sub_nonsense",
                                  "point_nonsense",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.deletion_NA",
                                  "rearr_deletion",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.inversion_NA",
                                  "rearr_inversion",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.tandem-duplication_NA",
                                  "rearr_tandem-duplication",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.translocation_NA",
                                  "rearr_translocation",driv.tnbc$variant_class)

driv.tnbc <- driv.tnbc[c("SampleID","Gene","variant_class","Subtype")]
names(driv.tnbc) <- c("sample","gene","mutation","subtype")

#######################################################################
# plot waterfall
#######################################################################

#myHierarchy<- data.table("mutation"=unique(driv.tnbc$mutation), 
#                         color=brewer.pal(length(unique(driv.tnbc$mutation)), "Set1"))

myHierarchy <- data.table(
  mutation = c("indel_frameshift", "point_nonsense", "point_missense", "indel_ess_splice",
               "point_ess_splice", 
               "indel_inframe", "CN_amplification", "rearr_deletion", "rearr_translocation", 
               "rearr_inversion", "rearr_tandem-duplication", "point_silent"),
  color = c("#E41A1C", "#FF7F00",  "#33A02C",  "#40E0D0", "#1F78B4",  "#ADD8E6",  "#6A3D9A",  
            "#F781BF",  "#89a832",  "#A65628",  "#FFFF00",  "#A6A6A6"))

plot.dat <- driv.tnbc[driv.tnbc$subtype=="TNBC_Basal",c("sample","gene","mutation")]

# add samples to the list that have no dirvers alterations, but are still part of the group
missing_samples <- setdiff(tnbc.samples[["TNBC_Basal"]], plot.dat$sample)
plot.dat <- rbind(plot.dat, data.frame(sample = missing_samples, gene = NA, mutation = NA))

plot <- GenVisR::Waterfall(plot.dat,
                           mutationHierarchy = myHierarchy,
                           plotA = NULL,
                           recurrence = 0,
                           geneMax = 10,
                           gridOverlay = FALSE,
                           sampleNames = TRUE,
                           drop = TRUE) 
#drawPlot(plot)
# append to list
plot.list <- append(plot.list,list(plot))

# 2
plot.dat <- driv.tnbc[driv.tnbc$subtype=="TNBC_NonBasal",c("sample","gene","mutation")]
missing_samples <- setdiff(tnbc.samples[["TNBC_NonBasal"]], plot.dat$sample)
plot.dat <- rbind(plot.dat, data.frame(sample = missing_samples, gene = NA, mutation = NA))

plot <- GenVisR::Waterfall(plot.dat,
                           mutationHierarchy = myHierarchy,
                           plotA = NULL,
                           recurrence = 0,
                           geneMax = 10,
                           gridOverlay = FALSE,
                           sampleNames = FALSE,
                           drop = TRUE) 

# append to list
plot.list <- append(plot.list,list(plot))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE, height = 10/2, width = 15/2)

for (i in 1:length(plot.list)) {
  drawPlot(plot.list[[i]])
  #grid::grid.newpage()
  #grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()
