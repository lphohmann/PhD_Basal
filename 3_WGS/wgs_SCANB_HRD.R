# Script: HRD in SCAN-B
# Author: Lennart Hohmann
# Date: 10.02.2024
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
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.4 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.5 <- "./data/Parameters/color_palette.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_HRD.pdf")
txt.file <- paste0(output.path,cohort,"_HRD.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load IDkey and correct sampleIDs -> ask Johan for key 
id.key <- loadRData(infile.4)
id.key <- id.key[c("Tumour","Specimen_id")]

# load palette
color.palette <- loadRData(infile.5)[c("LumA","LumB","Basal")]

# wgs QC 
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]
qc.samples.s <- id.key$Specimen_id[match(qc.samples,id.key$Tumour)]

# sampleid as rownames
hrd.dat <- as.data.frame(read_excel(infile.2, sheet = "HRDetect"))
hrd.dat <- hrd.dat[hrd.dat$Sample %in% qc.samples,]
hrd.dat$Sample <- id.key$Specimen_id[match(hrd.dat$Sample,id.key$Tumour)]
rownames(hrd.dat) <- hrd.dat$Sample
hrd.dat$Sample <- NULL

# HRD calling
hrd.dat$HRDetect <- ifelse(hrd.dat$Probability >= 0.7,"HRD-high","HRD-low")
#head(hrd.dat)
hrd.dat$PAM50 <- "Basal"
hrd.dat <- hrd.dat[c("HRDetect","PAM50")]

#######################################################################
# load and process BASIS data
#######################################################################

basis.anno <- loadRData(infile.3)
hrd.basis <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),
                         c("sample_name","HRDetect","PAM50_AIMS")]
hrd.basis <- hrd.basis[!is.na(hrd.basis$sample_name),]
row.names(hrd.basis) <- hrd.basis$sample_name
hrd.basis$sample_name <- NULL
names(hrd.basis) <- c("HRDetect","PAM50")

hrd.all <- rbind(hrd.dat,hrd.basis)
#table(hrd.all$PAM50,hrd.all$HRDetect)

#######################################################################

# total sample counts 
pam50.counts <- table(hrd.all$PAM50)

#######################################################################
# plot and stats
#######################################################################

# genes to plot/test

txt.out <- append(txt.out, c("\nHRD Frequency\n",
                             "\n###########################################\n"))

# gene data
hrd.tbl <- as.data.frame.matrix(table(hrd.all$PAM50,hrd.all$HRDetect))
hrd.freqs <- (hrd.tbl$`HRD-high` / pam50.counts)*100

# statistics
txt.out <- append(txt.out, c("\n",capture.output(hrd.freqs), "\n",
                             "\n###########################################\n"))

# mann whitney u tests
luma.res <- fisher.test(hrd.tbl[c("Basal","LumA"),])
lumb.res <- fisher.test(hrd.tbl[c("Basal","LumB"),])

txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))

# plot
plot.par <- list(
  height = hrd.freqs[c("LumA","LumB","Basal")], 
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
  axis(3,at=bp,labels=pam50.counts[names(color.palette)])
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)