# Script: Plotting and testing expression of metagenes in SCAN-B
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
               tidyverse,
               readxl)
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
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.4 <- "./data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX"
# output paths
plot.file <- paste0(output.path,cohort,"_metagenes.pdf")
txt.file <- paste0(output.path,cohort,"_metagenes.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# metagene definitions
metagene.def <- as.data.frame(read_excel(infile.4)) 
colnames(metagene.def) <- c("Module","Gene","Entrez_ID")

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# load palette
color.palette <- loadRData(infile.2)[c("LumA","LumB","Basal")]

# load gex
gex.df <- loadRData(infile.3)
gex.df <- gex.df[, unname(unlist(sampleIDs))]

#######################################################################
# calc. sample metagene scores
#######################################################################

metagenes <- unique(metagene.def$Module)
# for each sample calc. a score (mean) for each metagene
mg.scores <- apply(gex.df, 2, function(x) {
  res <- c()
  for (i in seq_along(metagenes)) {
    mg.genes <- metagene.def[metagene.def$Module==metagenes[i],"Gene"]
    mg.gex <- x[mg.genes]
    res[i] <- mean(mg.gex,na.rm=TRUE)
  }
  return(res)
})

rownames(mg.scores) <- metagenes

#######################################################################
# test and plot
#######################################################################





#######################################################################
# analyses
#######################################################################

gene.vec <- c("ERBB2","ESR1","FGFR4","TSPAN6","TNMD","DPM1","SCYL3","C1orf112","FGR")

for (gene in gene.vec) {
  
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  gene.gex <- gex.df[gene,] # misses pam 50 anno col
  
  # subtype data
  basal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  luma.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumA"]))]))
  lumb.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumB"]))]))
  
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  luma.stats <- get_stats(luma.dat)
  lumb.stats <- get_stats(lumb.dat)
  
  txt.out <- append(txt.out, c("Basal\n",capture.output(basal.stats), "\n",
                               "LumA\n",capture.output(luma.stats), "\n",
                               "LumB\n",capture.output(lumb.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- wilcox.test(basal.dat, luma.dat)
  lumb.res <- wilcox.test(basal.dat, lumb.dat)
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
                   data = list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
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
  boxplot(plot.parameters[[i]]$data,
          col = plot.parameters[[i]]$col,
          names = plot.parameters[[i]]$names,
          ylab = plot.parameters[[i]]$ylab,
          main = plot.parameters[[i]]$main)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
