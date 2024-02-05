# Script: Mutational & Rearrangement signatures in SCAN-B
# Author: Lennart Hohmann
# Date: 31.01.2024
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
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.4 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.5 <- "./data/Parameters/color_palette.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_signatures.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
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
sign.rearr <- as.data.frame(read_excel(infile.2, sheet = "RearrangmentSigs"))
rownames(sign.rearr) <- sign.rearr$Sample
sign.rearr$Sample <- NULL
sign.mut <- as.data.frame(read_excel(infile.2, sheet = "SBSsigs"))
rownames(sign.mut) <- sign.mut$Sample
sign.mut$Sample <- NULL

# add 6a and 6b
sign.rearr$RefSigR6 <- sign.rearr$RefSigR6a + sign.rearr$RefSigR6b
sign.rearr$RefSigR6a <- NULL
sign.rearr$RefSigR6b <- NULL

# convert SCANB to proportion per sample
sign.rearr <- as.data.frame(apply(sign.rearr,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.rearr[is.na(sign.rearr)] <- 0
sign.mut <- as.data.frame(apply(sign.mut,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.mut[is.na(sign.mut)] <- 0

# merge
sign.scanb <- as.data.frame(merge(t(sign.rearr), t(sign.mut), by = "row.names", all.x = TRUE))

# correct ids
sign.scanb$Sample <- id.key$Specimen_id[match(sign.scanb$Row.names,id.key$Tumour)]
sign.scanb$Row.names <- NULL
sign.scanb <- sign.scanb[sign.scanb$Sample %in% qc.samples.s,]

#######################################################################
# load and process BASIS data
#######################################################################

basis.anno <- loadRData(infile.3)
sign.basis <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),
                         c("PAM50_AIMS","RSproportions","MutSigProportions")]
# flatten
rs.df <- as.data.frame(sign.basis$RSproportions)
ms.df <- as.data.frame(sign.basis$MutSigProportions)
sign.basis <- cbind(sign.basis["PAM50_AIMS"], rs.df, ms.df)

# check which signatures are present in both cohorts
colnames(sign.scanb) <- gsub("RefSigR","RS", colnames(sign.scanb))
colnames(sign.scanb) <- gsub("SBS","S", colnames(sign.scanb))
common.sigs <- intersect(colnames(sign.scanb),colnames(sign.basis))

# exclude non-shared ones
sign.scanb <- sign.scanb[,common.sigs]
pam50.basis <- sign.basis$PAM50_AIMS
sign.basis <- sign.basis[,common.sigs]
sign.basis$PAM50 <- pam50.basis

# add subtype
sign.scanb$PAM50 <- "Basal"

# put together in an object ready for plotting 
sign.dat <- rbind(sign.scanb,sign.basis)
#View(sign.dat)

#######################################################################
# plot
#######################################################################

for (sig in common.sigs) {
  #sig="RS2"
  # sig data
  luma.dat <- sign.dat[sign.dat$PAM50=="LumA",sig]
  lumb.dat <- sign.dat[sign.dat$PAM50=="LumB",sig]
  basal.dat <- sign.dat[sign.dat$PAM50=="Basal",sig]
  
  # plot
  plot.par <- list(
    data = list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
    col = color.palette, 
    names = names(color.palette),
    ylab = "proportion",
    main = sig)
  # boxplot(plot_parameters$data, 
  #         col = plot_parameters$col,
  #         names = plot_parameters$names,
  #         ylab = plot_parameters$ylab,
  #         main = plot_parameters$main)
  # plot <- recordPlot()
  # plot.new()
  plot.parameters <- append(plot.parameters, list(plot.par))

}

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                main = plot.parameters[[i]]$main)
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)