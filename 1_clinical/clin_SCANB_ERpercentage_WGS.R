# Script: ER+ cell counts (%) in SCAN-B whole follow up cohort for basal wgs cases
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
infile.5 <- "./data/Parameters/TNBC_color_palette.RData"
infile.6 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.7 <- "./data/SCANB/2_transcriptomic/processed/All_LogScaled_gex.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_ERpercentage_WGS.pdf")
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
sampleIDs.2 <- loadRData(infile.6)[c("ERpHER2n_LumB", "ERpHER2n_LumA")]
sampleIDs <- c(sampleIDs.1, sampleIDs.2)

c1.ids <- c("S004527", "S001534", "S000006", "S003591", "S000299", "S005062")
#setdiff(sampleIDs,c1.ids)
c2.ids <- c("S000037","S000327","S001023","S001318","S003109",
            "S003433","S004629","S000939","S002552","S004454") 

# load palette
color.palette.1 <- loadRData(infile.5)[c("TNBC_NonBasal",
                                         "TNBC_Basal","ERpHER2n_Basal")]
color.palette.2 <- loadRData(infile.4)[c("LumB", "LumA")]
color.palette <- c(color.palette.1, color.palette.2)
names(color.palette)[names(color.palette) == "LumA"] <- "ERpHER2n_LumA"
names(color.palette)[names(color.palette) == "LumB"] <- "ERpHER2n_LumB"

# load gex
gex.df <- loadRData(infile.7)
gex.df <- gex.df[, unname(unlist(sampleIDs))]

## plot esr1

gene.gex <- gex.df["ESR1",] 

# subtype data
basal.dat.c1 <- as.numeric(as.vector(
  gene.gex[, c1.ids])) #adapt
basal.dat.c2 <- as.numeric(as.vector(
  gene.gex[, c2.ids])) #adapt



luma.dat <- as.numeric(as.vector(
  gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumA"]))]))
lumb.dat <- as.numeric(as.vector(
  gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumB"]))]))
tnbc.basal.dat <- as.numeric(as.vector(
  gene.gex[, unname(unlist(sampleIDs["TNBC_Basal"]))]))
tnbc.nonbasal.dat <- as.numeric(as.vector(
  gene.gex[, unname(unlist(sampleIDs["TNBC_NonBasal"]))]))
pdf(file = plot.file, onefile = TRUE)

par(mfrow = c(2, 2),las=2)
bp <- boxplot(list(ERpHER2n_LumA=luma.dat,
              ERpHER2n_LumB=lumb.dat,
              ERpHER2n_Basal_c1=basal.dat.c1,
              ERpHER2n_Basal_c2=basal.dat.c2,
              TNBC_Basal=tnbc.basal.dat,
              TNBC_NonBasal=tnbc.nonbasal.dat), 
  col = c("#2176d5","#34c6eb",
          "#f0fc03","#03fcb1", #adapt
          "#FF2400","#fec44f"), 
  names = c("ERpHER2n_LumA","ERpHER2n_LumB",
            "ERpHER2n_Basal_c1","ERpHER2n_Basal_c2","TNBC_Basal",
            "TNBC_NonBasal"),
  ylab = "mRNA expression (log2)",
  main = "ESR1")
axis(3,at=1:length(bp$n),labels=bp$n)

#######################################################################
# load data
#######################################################################

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
#anno <- anno[anno$Sample %in% sampleIDs, ]

# load data
erperc.dat <- read.table(infile.3,sep="\t",header=TRUE)
erperc.dat <- erperc.dat[erperc.dat$rba %in% anno$GEX.assay,]
erperc.dat <- erperc.dat[c("Specimen","INCA2_op_pad_erproc")]

#######################################################################
# Plot
#######################################################################

luma.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% unname(unlist(sampleIDs["ERpHER2n_LumA"]))]
lumb.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% unname(unlist(sampleIDs["ERpHER2n_LumB"]))]
basal.dat.c1 <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% c1.ids]
basal.dat.c2 <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% c2.ids]
tnbc.basal.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% unname(unlist(sampleIDs["TNBC_Basal"]))]
tnbc.nonbasal.dat <- erperc.dat$INCA2_op_pad_erproc[erperc.dat$Specimen %in% unname(unlist(sampleIDs["TNBC_NonBasal"]))]

#######################################################################

# save plots

bp <- boxplot(list(ERpHER2n_LumA=luma.dat,
                   ERpHER2n_LumB=lumb.dat,
                   ERpHER2n_Basal_c1=basal.dat.c1,
                   ERpHER2n_Basal_c2=basal.dat.c2,
                   TNBC_Basal=tnbc.basal.dat,
                   TNBC_NonBasal=tnbc.nonbasal.dat), 
              #ylim=c(-10,20),
              col = c("#2176d5","#34c6eb",
                      "#f0fc03","#03fcb1", #adapt
                      "#FF2400","#fec44f"), 
              names = c("ERpHER2n_LumA","ERpHER2n_LumB",
                        "ERpHER2n_Basal_c1","ERpHER2n_Basal_c2","TNBC_Basal",
                        "TNBC_NonBasal"),
              main = "ER+ cell count (%)")
#abline(h=0)
axis(3,at=1:length(bp$n),labels=bp$n)

dev.off()
