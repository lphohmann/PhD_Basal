# Script: Driver mutation waterfallplot in METABRIC
# Author: Lennart Hohmann
# Date: 30.01.2024
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
#source("./scripts/3_WGS/src/wgs_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,
               reshape2,
               GenVisR,
               data.table)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/3_WGS/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
infile.2 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_waterfall_ERpHER2nBasal.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load METABRIC data for ERpHER2nBasal samples
#######################################################################

driv.df <- loadRData(infile.1)
driv.df$PAM50 <- NULL
sample.ids <- loadRData(infile.2)
# make mut classsificaton unifrom

driv.df$variant_class[driv.df$variant_class == "CN_amp"] <- "CN_amplification"
driv.df$variant_class[driv.df$variant_class == "Frame_Shift_Del"] <- "indel_frameshift"
driv.df$variant_class[driv.df$variant_class == "Frame_Shift_Ins"] <- "indel_frameshift"
driv.df$variant_class[driv.df$variant_class == "Nonsense_Mutation"] <- "point_nonsense"
driv.df$variant_class[driv.df$variant_class == "Missense_Mutation"] <- "point_missense"
driv.df$variant_class[driv.df$variant_class == "In_Frame_Del"] <- "indel_inframe"
driv.df$variant_class[driv.df$variant_class == "In_Frame_Ins"] <- "indel_inframe"
driv.df$variant_class[driv.df$variant_class == "Silent"] <- "point_silent"
driv.df$variant_class[driv.df$variant_class == "Splice_Site"] <- "point_ess_splice"
driv.df$variant_class[driv.df$variant_class == "Splice_Region"] <- "point_ess_splice"

#######################################################################
# plot waterfall
#######################################################################
names(driv.df) <- c("sample", "gene", "mutation")

#myHierarchy <- data.table("mutation"=unique(driv.df$mutation), 
#                         color=brewer.pal(12, "Paired")[1:length(unique(driv.df$mutation))])
                         #color=brewer.pal(length(unique(driv.df$mutation)), "Set1"))

myHierarchy <- data.table(
  mutation = c("indel_frameshift", "point_nonsense", "point_missense", "indel_ess_splice",
               "point_ess_splice", 
               "indel_inframe", "CN_amplification", "rearr_deletion", "rearr_translocation", 
               "rearr_inversion", "rearr_tandem-duplication", "point_silent"),
  color = c("#E41A1C", "#FF7F00",  "#33A02C",  "#40E0D0", "#1F78B4",  "#ADD8E6",  "#6A3D9A",  
            "#F781BF",  "#89a832",  "#A65628",  "#FFFF00",  "#A6A6A6"))

# add samples to the list that have no dirvers alterations, but are still part of the group
missing_samples <- setdiff(sample.ids[["ERpHER2n_Basal"]], driv.df$sample)
driv.df <- rbind(driv.df, data.frame(sample = missing_samples, gene = NA, mutation = NA))

plot <- GenVisR::Waterfall(driv.df,
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
