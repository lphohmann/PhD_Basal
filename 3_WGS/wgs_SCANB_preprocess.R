# Script: Preprocess SCANB mut data into correct format
# Author: Lennart Hohmann
# Date: 11.11.2024
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
pacman::p_load(readxl,
               reshape2)#,GenVisR) # crashes if loaded with pacman
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GenVisR")
library(GenVisR)
#-------------------
# set/create output directories
# for data
data.path <- "./data/SCANB/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/SCANB/3_WGS/raw/Project2_Basal_like_Drivers_22Jan24_ForJohan.xlsx"
infile.4 <- "./data/BASIS/3_WGS/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx"
infile.5 <- "./data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.7 <- "./data/BASIS/3_WGS/raw/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx"
infile.8 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
# output paths
outfile.1 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
outfile.2 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nAll.RData"
outfile.3 <- "./data/SCANB/3_WGS/processed/driver_HugoSymbols.RData"
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# process SCANB ERpHER2nBasal samples
#######################################################################

# load sampleIDs
sampleIDs <- unname(unlist(loadRData(infile.1)[c("ERpHER2n_Basal")]))

# load IDkey and correct sampleIDs -> ask Johan for key 
id.key <- loadRData(infile.6)
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
driv.rearr[driv.rearr == "_"] <- NA
# where gene1 != gene2 i duplicate the row 
driv.rearr <- reshape2::melt(driv.rearr[c("sample","gene1","gene2","svclass")], id.vars = c("sample", "svclass"), measure.vars = c("gene1", "gene2"), variable.name = "gene_type", value.name = "gene")
driv.rearr <- driv.rearr[!is.na(driv.rearr$gene),]
driv.rearr <- driv.rearr[c("sample","gene","svclass")]
names(driv.rearr) <- c("Sample","VD_Gene","VC")
driv.rearr$VC <- paste0("rearr_", driv.rearr$VC)

# driver genes from BASIS to pull amplification status (CN drivers)
amp.drivers <- as.data.frame(read_excel(infile.4, sheet = "CopyNumber"))
driver.genes <- unique(c(amp.drivers$Gene,driv.indel$VD_Gene,driv.point$VD_Gene,
                         driv.rearr$VD_Gene))
driver.genes <- ifelse(driver.genes == "Chr8:(ZNF703/FGFR1)", "ZNF703", driver.genes)
save(driver.genes, file=outfile.3)
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

# process to right format
driv.df <- rbind(driv.indel[c("Sample","VD_Gene","VC")],
                 driv.point[c("Sample","VD_Gene","VC")],
                 driv.amp[c("Sample","VD_Gene","VC")],
                 driv.rearr)
names(driv.df) <- c("sample", "gene", "variant_class")
save(driv.df, file= outfile.1)

#View(driv.df)

#######################################################################
# process BASIS data
#######################################################################

# driver data SCANB
driv.scanb <- driv.df
driv.scanb$PAM50 <- "Basal"

# driver genes from BASIS LumA/LumB
driv.basis <- as.data.frame(read_excel(infile.7, 
                                       sheet = "COMBINED_EVENTS"))
basis.anno <- loadRData(infile.8)
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),]
driv.basis <- driv.basis[driv.basis$Sample %in% basis.anno$sample_name,]
driv.basis$PAM50 <- basis.anno$PAM50_AIMS[match(driv.basis$Sample,basis.anno$sample_name)]
driv.basis$variant_class <- paste0(driv.basis$Mutation_Type,"_",driv.basis$Effect)
driv.basis <- driv.basis[c("Sample","Gene","variant_class","PAM50")]
names(driv.basis) <- names(driv.scanb)

# exclude CopyNumber_HD & Complex_frameshift because these were not assessed in SCANB
driv.basis <- driv.basis[driv.basis$variant_class %!in% 
                           c("Complex_frameshift","CopyNumber_HD"),]

# in one df
driv.dat <- as.data.frame(rbind(driv.scanb,driv.basis))
#head(driv.dat) #sample   gene    variant_class PAM50
driv.dat$gene <- ifelse(driv.dat$gene == "Chr8:(ZNF703/FGFR1)", "ZNF703", 
                   driv.dat$gene)

save(driv.dat, file= outfile.2)
