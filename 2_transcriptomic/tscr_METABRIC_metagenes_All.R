# Script: All part - Plotting and testing expression of metagenes in MB
# Author: Lennart Hohmann
# Date: 08.01.2024
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
pacman::p_load(readxl)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.0 <- "./data/SCANB/2_transcriptomic/raw/metagene_definitions.XLSX"
infile.1 <- "./data/METABRIC/0_GroupSamples/TNBC_sampleIDs.RData"
infile.2 <- "./data/Parameters/TNBC_color_palette.RData"
infile.3 <- "./data/METABRIC/2_transcriptomic/processed/All_LogScaled_gex.RData" # Scaled mRNA expression
infile.4 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.5 <- "./data/Parameters/color_palette.RData"
# output paths
outfile.1 <- "./data/METABRIC/2_transcriptomic/processed/Metagene_scores_All.RData"
plot.file <- paste0(output.path,cohort,"_metagenes_All.pdf")
txt.file <- paste0(output.path,cohort,"_metagenes_All.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# metagene definitions
metagene.def <- as.data.frame(read_excel(infile.0)) 
colnames(metagene.def) <- c("Module","Gene","Entrez_ID")

# load sampleIDs
sampleIDs.1 <- loadRData(infile.1)[c("ERpHER2n_Basal", "TNBC_NonBasal", "TNBC_Basal")]
sampleIDs.2 <- loadRData(infile.4)[c("ERpHER2n_LumB", "ERpHER2n_LumA")]
sampleIDs <- c(sampleIDs.1, sampleIDs.2)

# load palette
color.palette.1 <- loadRData(infile.2)[c("TNBC_NonBasal",
                                         "TNBC_Basal","ERpHER2n_Basal")]
color.palette.2 <- loadRData(infile.5)[c("LumB", "LumA")]
color.palette <- c(color.palette.1, color.palette.2)
names(color.palette)[names(color.palette) == "LumA"] <- "ERpHER2n_LumA"
names(color.palette)[names(color.palette) == "LumB"] <- "ERpHER2n_LumB"

# load gex
gex.df <- loadRData(infile.3)
gex.df <- gex.df[, colnames(gex.df) %in% unname(unlist(sampleIDs))]

# how many sample sit it
sapply(sampleIDs, function(sublist) sum(names(gex.df) %in% sublist))
#unlist(lapply(sampleIDs,length))

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
  basal.dat <- as.numeric(as.vector(
    metagene.gex[, colnames(metagene.gex) %in% unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  luma.dat <- as.numeric(as.vector(
    metagene.gex[, colnames(metagene.gex) %in% unname(unlist(sampleIDs["ERpHER2n_LumA"]))]))
  lumb.dat <- as.numeric(as.vector(
    metagene.gex[, colnames(metagene.gex) %in% unname(unlist(sampleIDs["ERpHER2n_LumB"]))]))
  tnbc.basal.dat <- as.numeric(as.vector(
    metagene.gex[, colnames(metagene.gex) %in% unname(unlist(sampleIDs["TNBC_Basal"]))]))
  tnbc.nonbasal.dat <- as.numeric(as.vector(
    metagene.gex[, colnames(metagene.gex) %in% unname(unlist(sampleIDs["TNBC_NonBasal"]))]))
  
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  luma.stats <- get_stats(luma.dat)
  lumb.stats <- get_stats(lumb.dat)
  tnbc.basal.stats <- get_stats(tnbc.basal.dat)
  tnbc.nonbasal.stats <- get_stats(tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c("ERpHER2n_Basal\n",capture.output(basal.stats), "\n",
                               "ERpHER2n_LumA\n",capture.output(luma.stats), "\n",
                               "ERpHER2n_LumB\n",capture.output(lumb.stats), "\n",
                               "TNBC_Basal\n",capture.output(tnbc.basal.stats), "\n",
                               "TNBC_NonBasal\n",capture.output(tnbc.nonbasal.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- wilcox.test(basal.dat, luma.dat)
  lumb.res <- wilcox.test(basal.dat, lumb.dat)
  tnbc.basal.res <- wilcox.test(basal.dat, tnbc.basal.dat)
  tnbc.nonbasal.res <- wilcox.test(basal.dat, tnbc.nonbasal.dat)
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.basal.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(tnbc.nonbasal.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    data = list(ERpHER2n_LumA=luma.dat,
                ERpHER2n_LumB=lumb.dat,
                ERpHER2n_Basal=basal.dat,
                TNBC_Basal=tnbc.basal.dat,
                TNBC_NonBasal=tnbc.nonbasal.dat), 
    col = color.palette[c("ERpHER2n_LumA","ERpHER2n_LumB",
                          "ERpHER2n_Basal","TNBC_Basal",
                          "TNBC_NonBasal")], 
    names = c("ERpHER2n_LumA","ERpHER2n_LumB",
              "ERpHER2n_Basal","TNBC_Basal",
              "TNBC_NonBasal"),
    ylab = "metagene score",
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
                main = plot.parameters[[i]]$main,
                las = 2)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)