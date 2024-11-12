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

#######################################################################
# plot waterfall
#######################################################################
names(driv.df) <- c("sample", "gene", "mutation")

myHierarchy <- data.table("mutation"=unique(driv.df$mutation), 
                         color=brewer.pal(12, "Paired")[1:length(unique(driv.df$mutation))])
                         #color=brewer.pal(length(unique(driv.df$mutation)), "Set1"))

plot <- GenVisR::Waterfall(driv.df,
                           mutationHierarchy = myHierarchy,
                           plotA = NULL,
                           recurrence = 0,
                           geneMax = 10,
                           gridOverlay = FALSE,
                           sampleNames = FALSE,
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
