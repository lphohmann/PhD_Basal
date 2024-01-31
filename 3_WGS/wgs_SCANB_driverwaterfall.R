# Script: Driver mutation waterfallplot in SCAN-B
# Author: Lennart Hohmann
# Date: 30.01.2024
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
pacman::p_load(readxl,
               reshape2,
               GenVisR)
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
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/SCANB/3_WGS/raw/Project2_Basal_like_Drivers_22Jan24_ForJohan.xlsx"
infile.4 <- "./data/BASIS/3_WGS/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx"
infile.5 <- "data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_waterfall_ERpHER2nBasal.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load SCANB data for ERpHER2nBasal samples
#######################################################################

# load sampleIDs
sampleIDs <- unname(unlist(loadRData(infile.1)[c("ERpHER2n_Basal")]))

# driver data
driv.indel <- read_excel(infile.3, sheet = "PindelDrivers")
driv.point <- read_excel(infile.3, sheet = "CavemanDrivers")
driv.rearr <- read_excel(infile.3, sheet = "BRASS_drivers")
#driv.amp <- NULL # need ID key, so ask Johan

wgs.sum <- read_excel(infile.2, sheet = "Summary")
sign.rearr <- read_excel(infile.2, sheet = "RearrangmentSigs")
sign.mut <- read_excel(infile.2, sheet = "SBSsigs")
wgs.hrd <- read_excel(infile.2, sheet = "HRDetect")

# load IDkey and correct sampleIDs -> ask Johan for key 

# driver genes from BASIS to pull amplification status
sub.indel.drivers <- as.data.frame(read_excel(infile.4, sheet = "Subs.indels"))
amp.drivers <- as.data.frame(read_excel(infile.4, sheet = "CopyNumber"))
driver.genes <- unique(c(sub.indel.drivers$Gene,amp.drivers$Gene))
driver.genes <- ifelse(driver.genes == "Chr8:(ZNF703/FGFR1)", "ZNF703", driver.genes)

# amplification status
#cna.genes <- loadRData(infile.5)
#driv.amp <- x

#######################################################################
# process to right format
#######################################################################

driv.rearr[driv.rearr == "_"] <- NA
# where gene1 != gene2 i duplicate the row 
driv.rearr <- melt(driv.rearr[c("sample","gene1","gene2","svclass")], id.vars = c("sample", "svclass"), measure.vars = c("gene1", "gene2"), variable.name = "gene_type", value.name = "gene")
driv.rearr <- driv.rearr[!is.na(driv.rearr$gene),]
driv.rearr <- driv.rearr[c("sample","gene","svclass")]
names(driv.rearr) <- c("Sample","VD_Gene","VC")
driv.df <- rbind(driv.indel[c("Sample","VD_Gene","VC")],
      driv.point[c("Sample","VD_Gene","VC")],
      driv.rearr)
names(driv.df) <- c("sample", "gene", "variant_class")

#######################################################################
# plot waterfall
#######################################################################

plot <- waterfall(driv.df, 
                  fileType = "Custom", 
                  variant_class_order = unique(driv.df$variant_class), 
                  mainGrid = TRUE,
                  main_geneLabSize = 15,
                  mainRecurCutoff = 0,
                  maxGenes = 10,
                  mainDropMut = TRUE, # drop unused types from legend
                  #rmvSilent = TRUE,
                  out= "grob",
                  mutBurdenLayer = layer,
                  plotMutBurden = FALSE) #

# append to list
plot.list <- append(plot.list,list(plot))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE, height = 10, width = 15)

for (i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
  
  #print(plot.list[[i]])
}

dev.off()
