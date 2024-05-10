# Script: Kataegis in SCAN-B TNBC
# Author: Lennart Hohmann
# Date: 10.05.2024
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
#pacman::p_load(readxl,
#               reshape2)
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
infile.2 <- "./data/SCANB/3_WGS/raw/SCANB_ERp_collected_kataegis.RData"
infile.3 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.4 <- "./data/SCANB/5_TNBC_NatMed/SCANB_TNBC_collected_kataegis.RData"
infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.7 <- "./data/Parameters/TNBC_color_palette.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_kataegis_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_kataegis_TNBC.txt")
#-------------------
# storing objects
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.7)
sample.ids <- loadRData(infile.1)
# scanb data
basal.ids <- sample.ids[["ERpHER2n_Basal"]]
id.key <- loadRData(infile.3)
id.key <- id.key[id.key$Specimen.External.ID %in% basal.ids,]
kat.scanb <- loadRData(infile.2)
names(kat.scanb) <- gsub("^(.*?)(?:_vs_.*|$)", "\\1", names(kat.scanb))
kat.scanb <- kat.scanb[id.key$Tumour]
names(kat.scanb) <- id.key$Specimen.External.ID[match(names(kat.scanb),id.key$Tumour)]

#tnbc data
kat.tnbc <- loadRData(infile.4)
tnbc.anno <- loadRData(infile.5)
tnbc.samples <- sample.ids[c("TNBC_Basal","TNBC_NonBasal")]
tnbc.anno$Subtype <- sapply(tnbc.anno$External_ID_sample, function(sampleID) {
  for (sublist_name in names(tnbc.samples)) {
    if (sampleID %in% tnbc.samples[[sublist_name]]) {
      return(sublist_name)
    }
  }
  return(NA)
})

tnbc.anno <- tnbc.anno[c("PD_ID","Subtype")]
kat.tnbc <- kat.tnbc[tnbc.anno$PD_ID]

############################################################################

# get kataegis counts per tumor 
kat.scanb.counts <- lapply(kat.scanb, function(x) {
  x <- x[which(x$confidence>=1),] # filter by confidence
  return(nrow(x)) # number of events
})
# Convert the list to a data frame
kat.scanb.counts <- data.frame(
  SampleID = names(kat.scanb.counts),
  Kat_events = unlist(kat.scanb.counts)
)
kat.scanb.counts$Subtype <-"ERpHER2n_Basal"

# get kataegis counts per tumor 
kat.tnbc.counts <- lapply(kat.tnbc, function(x) {
  x <- x[which(x$confidence>=1),] # filter by confidence
  return(nrow(x)) # number of events
})
# Convert the list to a data frame
kat.tnbc.counts <- data.frame(
  SampleID = names(kat.tnbc.counts), 
  Kat_events = unlist(kat.tnbc.counts)
)

kat.tnbc.counts$Subtype <- tnbc.anno$Subtype[match(kat.tnbc.counts$SampleID,tnbc.anno$PD_ID)]

kat.counts <- rbind(kat.tnbc.counts,kat.scanb.counts)

# binary column
kat.counts$Kat_binary <- ifelse(kat.counts$Kat_events == 0, 0, 1)

#######################################################################
#######################################################################
subtype.counts <- table(kat.counts$Subtype)
  
#  to plot/test
txt.out <- append(txt.out, c("\nKataegis Frequency\n",
                             "\n###########################################\n"))

# gene data
kat.tbl <- as.data.frame.matrix(table(kat.counts$Subtype,kat.counts$Kat_binary))
kat.freqs <- (kat.tbl$`1` / subtype.counts)*100

# statistics
txt.out <- append(txt.out, c("\n",capture.output(kat.freqs), "\n",
                             "\n###########################################\n"))

# mann whitney u tests
bas.res <- fisher.test(kat.tbl[c("ERpHER2n_Basal","TNBC_Basal"),])
nonbas.res <- fisher.test(kat.tbl[c("ERpHER2n_Basal","TNBC_NonBasal"),])

txt.out <- append(txt.out, c(capture.output(bas.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(nonbas.res), "\n###########################################\n"))

# plot
plot.par <- list(
  height = kat.freqs[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")], 
  names=names(color.palette), 
  col = color.palette,
  ylim=c(0,100),
  main="Kataegis",
  ylab="Kataegis frequency (%)")
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