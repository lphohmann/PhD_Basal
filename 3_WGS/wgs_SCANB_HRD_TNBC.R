# Script: HRD in SCAN-B TNBC
# Author: Lennart Hohmann
# Date: 08.05.2024
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
pacman::p_load(readxl)
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
infile.1 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.5 <- "./data/Parameters/TNBC_color_palette.RData"
infile.7 <- "./data/SCANB/3_WGS/raw/2024_02_14_hrdetect_refsig_params_low_burden_sv_accounted_for.csv"
# output paths
plot.file <- paste0(output.path,cohort,"_HRD_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_HRD_TNBC.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load palette
color.palette <- loadRData(infile.5)

# wgs QC 
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]

# sampleid as rownames
hrd.dat <- read.table(infile.7, sep = ",", header = TRUE)
hrd.dat <- hrd.dat[hrd.dat$Tumour %in% qc.samples,]
hrd.dat$Lund.tumour.id <- gsub("\\..*", "", hrd.dat$Lund.tumour.id)

# HRD calling
hrd.dat$HRDetect <- ifelse(hrd.dat$Probability >= 0.7,"HRD-high","HRD-low")
#head(hrd.dat)
hrd.dat$Subtype <- "ERpHER2n_Basal"
hrd.dat <- hrd.dat[c("HRDetect","Subtype")]

#######################################################################
# load and process BASIS data
#######################################################################
tnbc.samples <- loadRData(infile.1)[c("TNBC_Basal","TNBC_NonBasal")]

tnbc.anno <- loadRData(infile.3)
tnbc.anno$Subtype <- sapply(tnbc.anno$External_ID_sample, function(sampleID) {
  for (sublist_name in names(tnbc.samples)) {
    if (sampleID %in% tnbc.samples[[sublist_name]]) {
      return(sublist_name)
    }
  }
  return(NA)
})
#View

hrd.tnbc <- tnbc.anno[tnbc.anno$External_ID_sample %in% unlist(tnbc.samples),
                       c("HRDetect.prob","Subtype")]
hrd.tnbc$HRDetect <- ifelse(hrd.tnbc$HRDetect.prob >= 0.7,"HRD-high","HRD-low")
hrd.tnbc$HRDetect.prob <- NULL
#row.names(hrd.basis) <- hrd.basis$sample_name
#hrd.basis$sample_name <- NULL
#names(hrd.basis) <- c("HRDetect","PAM50")

hrd.all <- rbind(hrd.dat,hrd.tnbc)
#table(hrd.all$Subtype,hrd.all$HRDetect)

#######################################################################

# total sample counts 
subtype.counts <- table(hrd.all$Subtype)

#######################################################################
# plot and stats
#######################################################################

#  to plot/test

txt.out <- append(txt.out, c("\nHRD Frequency\n",
                             "\n###########################################\n"))

# gene data
hrd.tbl <- as.data.frame.matrix(table(hrd.all$Subtype,hrd.all$HRDetect))
hrd.freqs <- (hrd.tbl$`HRD-high` / subtype.counts)*100

# statistics
txt.out <- append(txt.out, c("\n",capture.output(hrd.freqs), "\n",
                             "\n###########################################\n"))

# mann whitney u tests
bas.res <- fisher.test(hrd.tbl[c("ERpHER2n_Basal","TNBC_Basal"),])
nonbas.res <- fisher.test(hrd.tbl[c("ERpHER2n_Basal","TNBC_NonBasal"),])

txt.out <- append(txt.out, c(capture.output(bas.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(nonbas.res), "\n###########################################\n"))

# plot
plot.par <- list(
  height = hrd.freqs[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")], 
  names=names(color.palette), 
  col = color.palette,
  ylim=c(0,100),
  main="HRDetect",
  ylab="HRD-high frequency (%)")
plot.parameters <- append(plot.parameters, list(plot.par))

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
  axis(3,at=bp,labels=subtype.counts[names(color.palette)])
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)