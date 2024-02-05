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

# load IDkey and correct sampleIDs -> ask Johan for key 
id.key <- loadRData("./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData")
id.key <- id.key[c("Tumour","Specimen_id")]

# wgs QC 
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]
qc.samples.s <- id.key$Specimen_id[match(qc.samples,id.key$Tumour)]

# driver data
driv.indel <- read_excel(infile.3, sheet = "PindelDrivers")
driv.indel <- driv.indel[driv.indel$Sample %in% qc.samples,]
driv.indel$Sample <- id.key$Specimen_id[match(driv.indel$Sample,id.key$Tumour)]
driv.indel$VC <- paste0("indel_", driv.indel$VC)

driv.point <- read_excel(infile.3, sheet = "CavemanDrivers")
driv.point <- driv.point[driv.point$Sample %in% qc.samples,]
driv.point$Sample <- id.key$Specimen_id[match(driv.point$Sample,id.key$Tumour)]
driv.point$VC <- paste0("point_", driv.point$VC)

driv.rearr <- read_excel(infile.3, sheet = "BRASS_drivers")
driv.rearr <- driv.rearr[driv.rearr$sample %in% qc.samples,]
driv.rearr$sample <- id.key$Specimen_id[match(driv.rearr$sample,id.key$Tumour)]

#wgs.sum <- read_excel(infile.2, sheet = "Summary")
#sign.rearr <- read_excel(infile.2, sheet = "RearrangmentSigs")
#sign.mut <- read_excel(infile.2, sheet = "SBSsigs")
#wgs.hrd <- read_excel(infile.2, sheet = "HRDetect")

# driver genes from BASIS to pull amplification status (CN drivers)
#sub.indel.drivers <- as.data.frame(read_excel(infile.4, sheet = "Subs.indels"))
amp.drivers <- as.data.frame(read_excel(infile.4, sheet = "CopyNumber"))
driver.genes <- unique(c(amp.drivers$Gene,driv.indel$VD_Gene,driv.point$VD_Gene))
driver.genes <- ifelse(driver.genes == "Chr8:(ZNF703/FGFR1)", "ZNF703", driver.genes)

# amplification status
cna.genes <- loadRData(infile.5)
names(cna.genes) <- gsub("\\..*", "", names(cna.genes))
cna.genes <- cna.genes[names(cna.genes) %in% qc.samples.s]
cna.driver.genes <- lapply(cna.genes, function(x) {
  #x <- cna.genes[[1]]
  filt.x <- x[x$gene %in% driver.genes & x$Amp==1,]
  filt.x <- filt.x[,c("sample","gene")]
  filt.x$sample <- gsub("\\..*", "", filt.x$sample)
  return(filt.x)
})
driv.amp <- do.call(rbind, cna.driver.genes)
rownames(driv.amp) <- NULL
colnames(driv.amp) <- c("Sample","VD_Gene")
driv.amp$VC <- "CN_amplification"

#######################################################################
# process to right format
#######################################################################

driv.rearr[driv.rearr == "_"] <- NA
# where gene1 != gene2 i duplicate the row 
driv.rearr <- melt(driv.rearr[c("sample","gene1","gene2","svclass")], id.vars = c("sample", "svclass"), measure.vars = c("gene1", "gene2"), variable.name = "gene_type", value.name = "gene")
driv.rearr <- driv.rearr[!is.na(driv.rearr$gene),]
driv.rearr <- driv.rearr[c("sample","gene","svclass")]
names(driv.rearr) <- c("Sample","VD_Gene","VC")
driv.rearr$VC <- paste0("rearr_", driv.rearr$VC)
driv.df <- rbind(driv.indel[c("Sample","VD_Gene","VC")],
      driv.point[c("Sample","VD_Gene","VC")],
      driv.amp[c("Sample","VD_Gene","VC")],
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
                  maxGenes = 50,
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
