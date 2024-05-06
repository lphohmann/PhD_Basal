# Script: PAM50 correlations in SCAN-B
# Author: Lennart Hohmann
# Date: 26.01.2024
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
pacman::p_load("vioplot")
#-------------------
# set/create output directories
# for plots
output.path <- "./output/1_clinical/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/1_clinical/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/cst.txt"
infile.4 <- "./data/Parameters/color_palette.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_ERpercentage.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- unlist(loadRData(infile.1)[c("ERpHER2n_LumA",
                                          "ERpHER2n_LumB",
                                          "ERpHER2n_Basal")])

color.palette <- loadRData(infile.4)[c("LumA","LumB","Basal")]

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
#anno <- anno[anno$Sample %in% sampleIDs, ]

# load data
erperc.dat <- read.table(infile.3,sep="\t",header=TRUE)
erperc.dat <- erperc.dat[erperc.dat$rba %in% anno$GEX.assay,]
erperc.dat <- erperc.dat[c("Specimen","INCA2_op_pad_erproc")]
erperc.dat$PAM50 <- anno$PAM50_NCN[match(erperc.dat$Specimen,anno$Sample)] 
#View(erperc.dat)

#table(erperc.dat$INCA2_op_pad_erproc)

# also add the tnbc data to this
# load full follow up cohort
# then pull

#######################################################################
# Plot
#######################################################################

luma.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="LumA"]
lumb.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="LumB"]
basal.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="Basal"]

#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)
vioplot(list(LumA=luma.dat,
             LumB=lumb.dat,
             Basal=basal.dat),
        col = color.palette,
        names = names(color.palette),
        ylab = "INCA2_op_pad_erproc",
        main = "ER+ cell count (%)")
#points(jitter(rep(1, length(luma.dat)), amount = 0.2), luma.dat, col = "red")
#points(jitter(rep(2, length(lumb.dat)), amount = 0.2), lumb.dat, col = "red")
points(jitter(rep(3, length(basal.dat)), amount = 0.2), basal.dat, col = rgb(0, 0, 0, alpha = 0.7))
axis(3,at=1:3,labels=c(2790,1337,431))
abline(h = 10, col = "red", lty = 2)
text(1, 10, paste0("≤10% Basal_n=",sum(basal.dat<=10,na.rm=TRUE)), adj = c(0, -0.5), col = "red")
abline(h = 1, col = "black", lty = 2)
text(1, 1, paste0("≤1% Basal_n=",sum(basal.dat<=1,na.rm=TRUE)), adj = c(0, -0.5), col = "black")

dev.off()
