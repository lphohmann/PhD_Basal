# Script: TNBC part - Plotting and testing expression of metagenes in SCAN-B
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
pacman::p_load(readxl)
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
#infile.1 <- "./data/SCANB/5_TNBC_NatMed/ASCAT_InSilico_SCAN_B_TNBC_EasySegments.RData"
#infile.2 <- "./data/SCANB/5_TNBC_NatMed/driver.df.scanb.complete.csv"
#infile.3 <- "./data/SCANB/5_TNBC_NatMed/Final_fData_Object.RData"
#infile.4 <- "./data/SCANB/5_TNBC_NatMed/SCANB_TNBC_collected_kataegis.RData"
#infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.6 <- "./data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX"
infile.7 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.8 <- "./data/Parameters/TNBC_color_palette.RData"
infile.9 <- "./data/SCANB/2_transcriptomic/processed/All_LogScaled_gex.RData" # Scaled mRNA expression
# output paths
outfile.1 <- "./data/SCANB/2_transcriptomic/processed/Metagene_scores_TNBC.RData"
plot.file <- paste0(output.path,cohort,"_metagenes_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_metagenes_TNBC.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# metagene definitions
metagene.def <- as.data.frame(read_excel(infile.6)) 
colnames(metagene.def) <- c("Module","Gene","Entrez_ID")

# load sampleIDs
sampleIDs <- loadRData(infile.7)

# load palette
color.palette <- loadRData(infile.8)

# load gex
gex.df <- loadRData(infile.9)
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
mg.scores <- as.data.frame(mg.scores)
save(mg.scores,file=outfile.1)

#View(mg.scores)

#######################################################################
# test and plot
#######################################################################

for (metagene in metagenes) {
  
  txt.out <- append(txt.out, c("\n",metagene,"\n",
                               "\n###########################################\n"))
  
  metagene.gex <- mg.scores[metagene,] 
  
  # subtype data
  erp.basal.dat <- as.numeric(as.vector(
    metagene.gex[, unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  tnbc.basal.dat <- as.numeric(as.vector(
    metagene.gex[, unname(unlist(sampleIDs["TNBC_Basal"]))]))
  tnbc.nonbasal.dat <- as.numeric(as.vector(
    metagene.gex[, unname(unlist(sampleIDs["TNBC_NonBasal"]))]))
  
  # summary statistics
  erp.basal.stats <- get_stats(erp.basal.dat)
  tnbc.basal.stats <- get_stats(tnbc.basal.dat)
  tnbc.nonbasal.stats <- get_stats(tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c("ERp.Basal\n",capture.output(erp.basal.stats), "\n",
                               "TNBC.Basal\n",capture.output(tnbc.basal.stats), "\n",
                               "TNBC.NonBasal\n",capture.output(tnbc.nonbasal.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  tnbc.basal.res <- wilcox.test(erp.basal.dat, tnbc.basal.dat)
  tnbc.nonbasal.res <- wilcox.test(erp.basal.dat, tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c(capture.output(tnbc.basal.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.nonbasal.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    data = list(TNBC_NonBasal=tnbc.nonbasal.dat,TNBC_Basal=tnbc.basal.dat,ERpHER2n_Basal=erp.basal.dat), 
    col = color.palette, 
    names = names(color.palette),
    ylab = "mRNA expression (log2)",
    main = metagene)
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