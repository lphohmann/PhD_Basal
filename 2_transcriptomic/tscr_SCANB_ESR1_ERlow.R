# Script: boxplot esr1 expression in tnbc and er low
# Author: Lennart Hohmann
# Date: 27.03.2025
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
# pacman::p_load(ggplot2,
#                tidyverse)
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
infile.4 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.5 <- "./data/Parameters/color_palette.RData"
infile.6 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.7 <- "./data/SCANB/2_transcriptomic/processed/All_LogUnscaled_gex.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_ESR1_ERlow.pdf")
txt.file <- paste0(output.path,cohort,"_ESR1_ERlow.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs.1 <- loadRData(infile.1)[c("ERpHER2n_Basal", "TNBC_NonBasal", "TNBC_Basal")]
sampleIDs.2 <- loadRData(infile.4)[c("ERpHER2n_LumB", "ERpHER2n_LumA")]
sampleIDs <- c(sampleIDs.1, sampleIDs.2)

# load palette
color.palette <- c(
  TNBC_1t10 = "#9c27b0",   # Vibrant Purple  
  TNBC_u1 = "#6a0dad",      # Deep Purple  
  ERpHER2n_Basal = "#ba0606",  # Red (unchanged)  
  ERpHER2n_LumB = "#34c6eb",   # Light Blue (unchanged)  
  ERpHER2n_LumA = "#2176d5"    # Dark Blue (unchanged)  
)

# load gex
gex.df.scaled <- loadRData(infile.3)
gex.df.scaled <- gex.df.scaled[, unname(unlist(sampleIDs))]

gex.df.unscaled <- loadRData(infile.7)
gex.df.unscaled <- gex.df.unscaled[, unname(unlist(sampleIDs))]

# load the tnbc cohort anno file
tnbc.anno <- loadRData(infile.6)
tnbc.anno <- tnbc.anno[tnbc.anno$External_ID_sample %in% unlist(sampleIDs[c("TNBC_NonBasal", "TNBC_Basal")]),]
tnbc.anno <- tnbc.anno[c("External_ID_sample","ERperc_group")]
tnbc.anno <- tnbc.anno[!is.na(tnbc.anno$ERperc_group),] # 1 sample NA

tnbc_u1.ids <- tnbc.anno$External_ID_sample[tnbc.anno$ERperc_group == "<1%"]
tnbc_1t10.ids <- tnbc.anno$External_ID_sample[tnbc.anno$ERperc_group == "1-10%"]

#######################################################################
# analyses
#######################################################################
gex.list <- list("scaled"=gex.df.scaled,"unscaled"=gex.df.unscaled)
for (i in 1:length(gex.list)) {
  
  gex.df <- gex.list[[i]]
  gene <- "ESR1"
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  gene.gex <- gex.df[gene,] 
  
  # subtype data
  basal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  luma.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumA"]))]))
  lumb.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumB"]))]))
  tnbc.u1.dat <- as.numeric(as.vector(gene.gex[, tnbc_u1.ids]))
  tnbc.1t10.dat <- as.numeric(as.vector(gene.gex[, tnbc_1t10.ids]))
  
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  luma.stats <- get_stats(luma.dat)
  lumb.stats <- get_stats(lumb.dat)
  tnbc.u1.stats <- get_stats(tnbc.u1.dat)
  tnbc.1t10.stats <- get_stats(tnbc.1t10.dat)
  
  txt.out <- append(txt.out, c("ERpHER2n_Basal\n",capture.output(basal.stats), "\n",
                               "ERpHER2n_LumA\n",capture.output(luma.stats), "\n",
                               "ERpHER2n_LumB\n",capture.output(lumb.stats), "\n",
                               "TNBC_Basal\n",capture.output(tnbc.u1.stats), "\n",
                               "TNBC_NonBasal\n",capture.output(tnbc.1t10.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- wilcox.test(basal.dat, luma.dat)
  lumb.res <- wilcox.test(basal.dat, lumb.dat)
  tnbc.u1.res <- wilcox.test(basal.dat, tnbc.u1.dat)
  tnbc.1t10.res <- wilcox.test(basal.dat, tnbc.1t10.dat)
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.u1.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.1t10.res), "\n###########################################\n"))
  
  # plot c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")
  plot.par <- list(
    data = list(ERpHER2n_LumA=luma.dat,
                ERpHER2n_LumB=lumb.dat,
                ERpHER2n_Basal=basal.dat,
                TNBC_1t10=tnbc.1t10.dat,
                TNBC_u1=tnbc.u1.dat
                ), 
    col = color.palette[c("ERpHER2n_LumA","ERpHER2n_LumB",
                          "ERpHER2n_Basal","TNBC_1t10","TNBC_u1"
                          )], 
    names = c("ERpHER2n_LumA","ERpHER2n_LumB",
              "ERpHER2n_Basal","TNBC_1t10","TNBC_u1"),
    ylab = paste0(names(gex.list)[i]," mRNA expression (log2)"),
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
                main = plot.parameters[[i]]$main,
                las = 2)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
