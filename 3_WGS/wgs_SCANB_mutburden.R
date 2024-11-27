# Script: plot TMB SCAN-B vs BASIS
# Author: Lennart Hohmann
# Date: 25.11.2024
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
#source("./scripts/3_WGS/src/")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl)
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
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGSanno.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/3_WGS/raw/HRDproject_table_S1.xlsx" 
# output paths
plot.file <- paste0(output.path,cohort,"_TMB.pdf")
txt.file <- paste0(output.path,cohort,"_TMB.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################
color.palette <- loadRData(infile.2)
sample.IDs <- loadRData(infile.1)
scanb.dat <- read_excel(infile.3,sheet="Table_S1")
scanb.dat <- scanb.dat[which(scanb.dat$Specimen %in% sample.IDs$SampleID),c("Specimen","PAM50_NCN","Caveman.counts_CLPM.0_ASMD140_final","Pindel.counts_QUAL250_Repeats_10_final")]

# basis mut data
basis.dat <- loadRData("data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData")
basis.dat <- basis.dat[which(basis.dat$ClinicalGroup == "ERposHER2neg" & 
                               basis.dat$PAM50_AIMS %in% c("LumA","LumB")),c("sample_name","PAM50_AIMS","nbrSubs","nbrIndels")]

names(basis.dat) <- names(scanb.dat)
# combined
mut.counts <- rbind(scanb.dat,basis.dat)
#View(mut.counts)
mut.counts$N_mut <- mut.counts$Caveman.counts_CLPM.0_ASMD140_final + mut.counts$Pindel.counts_QUAL250_Repeats_10_final

# Mutations Per Megabase (Mut/Mb)
mut.counts$Mut_per_Mb <- mut.counts$N_mut / genomic_size_mb


# plot

txt.out <- append(txt.out, c("\n","Mutational burden","\n",
                             "\n###########################################\n"))

# subtype data
basal.dat <- mut.counts[mut.counts$PAM50_NCN == "Basal",]$N_mut
lumb.dat <- mut.counts[mut.counts$PAM50_NCN == "LumB",]$N_mut
luma.dat <- mut.counts[mut.counts$PAM50_NCN == "LumA",]$N_mut

# here

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
  ylab = "N_mut",
  main = "TMB")
# boxplot(plot_parameters$data, 
#         col = plot_parameters$col,
#         names = plot_parameters$names,
#         ylab = plot_parameters$ylab,
#         main = plot_parameters$main)
# plot <- recordPlot()
# plot.new()
plot.parameters <- append(plot.parameters, list(plot.par))
#plot.list <- append(plot.list, list(plot))

#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                ylim = c(0, 30000),
                main = plot.parameters[[i]]$main)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)