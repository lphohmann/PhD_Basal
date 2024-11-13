# Script: Driver mutation frequecies in METABRIC
# Author: Lennart Hohmann
# Date: 12.11.2024
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
               reshape2)
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
infile.1 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nAll.RData"
infile.7 <- "./data/Parameters/color_palette.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_mutfreqs.pdf")
txt.file <- paste0(output.path,cohort,"_mutfreqs.txt")
#-------------------
# storing objects
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load palette
color.palette <- loadRData(infile.7)[c("LumA","LumB","Basal")]

driv.dat <- loadRData(infile.1)
#View(driv.dat)

# total sample counts to calc. mut freqs
pam50.counts <- table(
  driv.dat[!duplicated(driv.dat[,c("sample")]),]$PAM50)

# need to calc based on all samples that have mut info, not only ones with driver alterations
#Basal 45  LumA 614  LumB 397
pam50.counts[] <- c(45,614,397)

#######################################################################
# plot and stats
#######################################################################

# genes to plot/test
genes <- c("TP53","MYC","PIK3CA","BRCA2","ESR1")
for (gene in genes) {
  #print(gene)
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  # gene data
  gene.dat <- as.data.frame(table(driv.dat[driv.dat$gene==gene,]$PAM50))
  row.names(gene.dat) <- gene.dat$Var1
  gene.dat$Var1 <- NULL
  gene.dat$Not_mutated <- pam50.counts[row.names(gene.dat)] - gene.dat$Freq
  names(gene.dat) <- c("Mutated","Not_Mutated")
  gene.mut.freqs <- (gene.dat$Mutated / pam50.counts[row.names(gene.dat)])*100
  
  # statistics
  txt.out <- append(txt.out, c("\n",capture.output(gene.mut.freqs), "\n",
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- if(sum(is.na(gene.dat[c("Basal","LumA"),]))==0) {
    fisher.test(gene.dat[c("Basal","LumA"),]) } else {"NA"}
  lumb.res <- if(sum(is.na(gene.dat[c("Basal","LumB"),]))==0) {
    fisher.test(gene.dat[c("Basal","LumB"),]) } else {"NA"}
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    height = gene.mut.freqs[c("LumA","LumB","Basal")], 
    names=names(color.palette), 
    col = color.palette,
    ylim=c(0,100),
    main=gene,
    ylab="Mutation frequency (%)")
  plot.parameters <- append(plot.parameters, list(plot.par))
}

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- barplot(height=plot.parameters[[i]]$height, 
                names=plot.parameters[[i]]$names, 
                col =plot.parameters[[i]]$col,
                ylim=plot.parameters[[i]]$ylim,
                main=plot.parameters[[i]]$main,
                ylab=plot.parameters[[i]]$ylab)
  axis(3,at=bp,labels=pam50.counts[names(color.palette)])
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
