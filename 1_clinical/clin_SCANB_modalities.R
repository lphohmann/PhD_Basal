# Script: overview avaialble modalities SCAN-B
# Author: Lennart Hohmann
# Date: 20.01.2025
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
pacman::p_load(VennDiagram)
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
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_modalitiesVenn.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]

#View(anno)
anno <- anno[c("Sample","ER","PR","HER2","LN.spec",
     "NHG","Size.mm","TreatGroup","DRFi_days",
     "Age",  "OSbin","OS","RFIbin","RFI",
     "DRFIbin","DRFI","NCN.PAM50")]               

table(is.na(anno$Size.mm))

pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))

# Create the Venn diagrams and display them
par(mfrow = c(2, 2)) # Arrange the plots in a 2x2 grid

# FoxA1 Meth vs. FoxC1 Meth
venn1 <- venn.diagram(
  x = list(FoxA1_Meth = foxa1.meth, FoxC1_Meth = foxc1.meth),
  category.names = c("FoxA1 Meth", "FoxC1 Meth"),
  filename = NULL,
  main = "FoxA1 Meth vs FoxC1 Meth"
)
grid.newpage()
grid.draw(venn1)

dev.off()