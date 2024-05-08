# Script: Mutational & Rearrangement signatures in SCAN-B TNBC
# Author: Lennart Hohmann
# Date: 07.05.2024
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
infile.4 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.5 <- "./data/Parameters/TNBC_color_palette.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/2024_02_11_snv_sig_refsig_errperc20_correct_exposures_combined_502.csv"
infile.9 <- "./data/SCANB/3_WGS/raw/2024_02_14_SCANB_ERpos_rearr_sig_502.csv"
# output paths
plot.file <- paste0(output.path,cohort,"_signatures_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_signatures_TNBC.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load new data
#######################################################################

## check
#infile.c <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
#c <- loadRData(infile.c) # big anno obj with sig data
#c$SigFit.mutationalSignatures.prop
#c$RearrangementSignatures.prop



##

# load IDkey and correct sampleIDs 
id.key <- loadRData(infile.4)
id.key <- id.key[c("Tumour","Specimen_id")]

# load palette
color.palette <- loadRData(infile.5)

# wgs QC 
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]
qc.samples.s <- id.key$Specimen_id[match(qc.samples,id.key$Tumour)]

sign.mut <- read.table(infile.6, sep = ",", header = TRUE)
sign.mut <- sign.mut[c("Lund.tumour.id","SBS1","SBS2","SBS3","SBS5",
                     "SBS8","SBS13","SBS17","SBS18","SBS127", 
                     "unassigned","SBS6")]
names(sign.mut)[1] <- "Sample"
sign.mut$Sample <- gsub("\\..*", "", sign.mut$Sample)
rownames(sign.mut) <- sign.mut$Sample
sign.mut$Sample <- NULL

sign.rearr <- read.table(infile.9, sep = ",", header = TRUE)
sign.rearr <- sign.rearr[c("Lund.tumour.id","RefSigR2", "RefSigR4","RefSigR6a",
                         "RefSigR1","RefSigR5","RefSigR6b","RefSigR8","RefSigR3","RefSigR11",
                         "unassigned")]
names(sign.rearr)[1] <- "Sample"
sign.rearr$Sample <- gsub("\\..*", "", sign.rearr$Sample)
rownames(sign.rearr) <- sign.rearr$Sample
sign.rearr$Sample <- NULL

# add 6a and 6b
sign.rearr$RefSigR6 <- sign.rearr$RefSigR6a + sign.rearr$RefSigR6b
sign.rearr$RefSigR6a <- NULL
sign.rearr$RefSigR6b <- NULL

# which samples have <25 rearr snvs
sign.rearr.exclude <- sign.rearr
sign.rearr.exclude <- sign.rearr.exclude[qc.samples.s,]
sign.rearr.exclude$sum <-apply(sign.rearr.exclude,1,sum) #exclude and redo
sign.rearr.exclude <- row.names(sign.rearr.exclude[sign.rearr.exclude$sum < 25,])

# convert SCANB to proportion per sample
sign.rearr <- as.data.frame(apply(sign.rearr,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.rearr[is.na(sign.rearr)] <- 0
sign.rearr <- sign.rearr[row.names(sign.rearr) != "unassigned",]

sign.mut <- as.data.frame(apply(sign.mut,1,function(x) (x/sum(x, na.rm=TRUE))))
sign.mut[is.na(sign.mut)] <- 0
sign.mut <- sign.mut[row.names(sign.mut) != "unassigned",]

# merge
sign.scanb <- as.data.frame(merge(t(sign.rearr), t(sign.mut), by = "row.names", all.x = TRUE))
names(sign.scanb)[1] <- "Sample"

# select basal samples
sign.scanb <- sign.scanb[sign.scanb$Sample %in% qc.samples.s,]

#######################################################################
# load and process TNBC data
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

# cont here

sign.tnbc <- tnbc.anno[tnbc.anno$External_ID_sample %in% unlist(tnbc.samples),
                         c("Subtype","SigFit.mutationalSignatures.prop","RearrangementSignatures.prop")]
# flatten
rs.df <- as.data.frame(sign.tnbc$RearrangementSignatures.prop)
ms.df <- as.data.frame(sign.tnbc$SigFit.mutationalSignatures.prop)
sign.tnbc <- cbind(sign.tnbc["Subtype"], rs.df, ms.df)

str(sign.tnbc)
str(sign.scanb)

# check which signatures are present in both cohorts
colnames(sign.scanb) <- gsub("RefSigR","RS", colnames(sign.scanb))
colnames(sign.scanb) <- gsub("SBS","S", colnames(sign.scanb))
common.sigs <- intersect(colnames(sign.scanb),colnames(sign.tnbc))

# exclude non-shared ones
rownames.basal <- sign.scanb$Sample
sign.scanb <- sign.scanb[,common.sigs]
st.tnbc <- sign.tnbc$Subtype
sign.tnbc <- sign.tnbc[,common.sigs]
sign.tnbc$Subtype <- st.tnbc

# add subtype
sign.scanb$Subtype <- "ERpHER2n_Basal"

# put together in an object ready for plotting 
sign.dat <- rbind(sign.scanb,sign.tnbc)
#View(sign.dat)
#row.names(sign.dat)[1:length(rownames.basal)] <- rownames.basal

#######################################################################
# plot
#######################################################################
#sig="RS2"
for (sig in common.sigs) {
  
  txt.out <- append(txt.out, c("\n",sig,"\n",
                               "\n###########################################\n"))
  
  # sig data
  bas.dat <- sign.dat[sign.dat$Subtype=="TNBC_Basal",sig]
  nonbas.dat <- sign.dat[sign.dat$Subtype=="TNBC_NonBasal",sig]
  # have to filter samples for rearr signatures below <25
  if (grepl("RS", sig)) {
    basal.dat <- sign.dat[sign.dat$Subtype=="ERpHER2n_Basal",]
    basal.dat <- basal.dat[!row.names(basal.dat) %in% sign.rearr.exclude,sig]
  } else { basal.dat <- sign.dat[sign.dat$Subtype=="ERpHER2n_Basal",sig] }
  
  # statistics
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  bas.stats <- get_stats(bas.dat)
  nonbas.stats <- get_stats(nonbas.dat)
  
  txt.out <- append(txt.out, c("ERpHER2n_Basal\n",capture.output(basal.stats), "\n",
                               "TNBC_Basal\n",capture.output(bas.stats), "\n",
                               "TNBC_NonBasal\n",capture.output(nonbas.stats),
                               "\n###########################################\n"))
  
  # mann whitney u tests
  bas.res <- wilcox.test(basal.dat, bas.dat)
  nonbas.res <- wilcox.test(basal.dat, nonbas.dat)
  
  txt.out <- append(txt.out, c(capture.output(bas.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(nonbas.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    data = list(TNBC_NonBasal=nonbas.dat,TNBC_Basal=bas.dat,ERpHER2n_Basal=basal.dat), 
    col = color.palette, 
    names = names(color.palette),
    ylab = "proportion",
    main = sig)
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