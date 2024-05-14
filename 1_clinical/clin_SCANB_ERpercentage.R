# Script: ER+ cell counts (%) in SCAN-B whole follow up cohort
# Author: Lennart Hohmann
# Date: 06.05.2024
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
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
sampleIDs <- loadRData(infile.1)

#color.palette <- loadRData(infile.4)[c("LumA","LumB","Basal")]
color.palette <- c(LumA="#2176d5", LumB="#34c6eb", Basal="#ba0606", TNBCBasal_under1p = "#fff7bc",
  TNBCBasal_1to10p = "#fec44f",
  ERpBasal = "#d95f0e")
# add 3 colors
#color.palette

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
#anno <- anno[anno$Sample %in% sampleIDs, ]

# load data
erperc.dat <- read.table(infile.3,sep="\t",header=TRUE)
erperc.dat <- erperc.dat[erperc.dat$rba %in% anno$GEX.assay,]
erperc.dat <- erperc.dat[c("Specimen","INCA2_op_pad_erproc")]
erperc.dat$PAM50 <- anno$PAM50_NCN[match(erperc.dat$Specimen,anno$Sample)] 

#######################################################################
# Plot
#######################################################################

luma.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="LumA"]
lumb.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="LumB"]
basal.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$PAM50=="Basal"]

tnbc.basal.1 <- erperc.dat[erperc.dat$PAM50=="Basal" & 
                             erperc.dat$INCA2_op_pad_erproc < 1 ,] #$INCA2_op_pad_erproc

tnbc.basal.1.10 <- erperc.dat[erperc.dat$PAM50=="Basal" & 
                                erperc.dat$INCA2_op_pad_erproc >= 1 & 
                                erperc.dat$INCA2_op_pad_erproc < 10,] #$INCA2_op_pad_erproc
erp.basal <- erperc.dat[erperc.dat$PAM50=="Basal" & 
                          erperc.dat$INCA2_op_pad_erproc >= 10,]
#View(tnbc.basal.1.10)
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE)
par(las=2)
vioplot(list(LumA=luma.dat,
             LumB=lumb.dat,
             Basal=basal.dat,
             TNBCBasal_under1p = tnbc.basal.1$INCA2_op_pad_erproc,
             TNBCBasal_1to10p = tnbc.basal.1.10$INCA2_op_pad_erproc,
             ERpBasal = erp.basal$INCA2_op_pad_erproc),
        col = color.palette,
        names = names(color.palette),
        ylab = "INCA2_op_pad_erproc",
        main = "ER+ cell count (%)")
#points(jitter(rep(1, length(luma.dat)), amount = 0.2), luma.dat, col = "red")
#points(jitter(rep(2, length(lumb.dat)), amount = 0.2), lumb.dat, col = "red")
points(jitter(rep(3, length(basal.dat)), amount = 0.2), basal.dat, col = rgb(0, 0, 0, alpha = 0.7))
points(jitter(rep(4, length(tnbc.basal.1$INCA2_op_pad_erproc)), amount = 0.2), tnbc.basal.1$INCA2_op_pad_erproc, col = rgb(0, 0, 0, alpha = 0.7))

axis(3,at=1:6,labels=c(sum(!is.na(luma.dat)),
                       sum(!is.na(luma.dat)),
                       sum(!is.na(basal.dat)),
                       sum(!is.na(tnbc.basal.1$INCA2_op_pad_erproc)),
                       sum(!is.na(tnbc.basal.1.10$INCA2_op_pad_erproc)),
                       sum(!is.na(erp.basal$INCA2_op_pad_erproc))))
abline(h = 10, col = "red", lty = 2)
text(1, 10, "10%", adj = c(0, -0.5), col = "red")
abline(h = 1, col = "black", lty = 2)
text(1, 1, "1%", adj = c(0, -0.5), col = "black")

dev.off()
