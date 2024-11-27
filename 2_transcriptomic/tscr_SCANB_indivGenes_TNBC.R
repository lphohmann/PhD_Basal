# Script: TNBC - Plotting and testing expression of indiv. selected genes in SCAN-B
# Author: Lennart Hohmann
# Date: 08.01.2024
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
pacman::p_load(ggplot2,
               tidyverse)
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.2 <- "./data/Parameters/TNBC_color_palette.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/All_LogScaled_gex.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_indivGenes_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_indivGenes_TNBC.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "TNBC_NonBasal", "TNBC_Basal")]

# load palette
color.palette <- loadRData(infile.2)[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")]

# load gex
gex.df <- loadRData(infile.3)
gex.df <- gex.df[, unname(unlist(sampleIDs))]

#######################################################################
# analyses
#######################################################################

gene.vec <- c("ESR1","AR","PGR","STAT3","STAT6","ABAT","AGR2","CA12","DNALI1","FOXA1","GATA3","SLC44A4","TBC1D9","XBP1","FOXC1")
#setdiff(gene.vec, rownames(gex.df))
#gene <- "ESR1"
for (gene in gene.vec) {
  print(gene)
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  gene.gex <- gex.df[gene,] 
  
  # subtype data
  basal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  tnbc.basal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["TNBC_Basal"]))]))
  tnbc.nonbasal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["TNBC_NonBasal"]))]))
  
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  tnbc.basal.stats <- get_stats(tnbc.basal.dat)
  tnbc.nonbasal.stats <- get_stats(tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c("ERpHER2n_Basal\n",capture.output(basal.stats), "\n",
                               "TNBC_Basal\n",capture.output(tnbc.basal.stats), "\n",
                               "TNBC_NonBasal\n",capture.output(tnbc.nonbasal.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  tnbc.basal.res <- wilcox.test(basal.dat, tnbc.basal.dat)
  tnbc.nonbasal.res <- wilcox.test(basal.dat, tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c(capture.output(tnbc.basal.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.nonbasal.res), "\n###########################################\n"))
  
  # plot c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")
  plot.par <- list(
                   data = list(TNBC_NonBasal=tnbc.nonbasal.dat,
                               TNBC_Basal=tnbc.basal.dat,
                               ERpHER2n_Basal=basal.dat), 
                   col = color.palette, 
                   names = names(color.palette),
                   ylab = "mRNA expression (log2)",
                   main = gene)
  # boxplot(plot_parameters$data, 
  #         col = plot_parameters$col,
  #         names = plot_parameters$names,
  #         ylab = plot_parameters$ylab,
  #         main = plot_parameters$main)
  # plot <- recordPlot()
  # plot.new()
  plot.parameters <- append(plot.parameters, list(plot.par))
  #plot.list <- append(plot.list, list(plot))
}
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
          col = plot.parameters[[i]]$col,
          names = plot.parameters[[i]]$names,
          ylab = plot.parameters[[i]]$ylab,
          main = plot.parameters[[i]]$main)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
