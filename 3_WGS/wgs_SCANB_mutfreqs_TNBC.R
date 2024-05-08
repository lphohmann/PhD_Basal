# Script: Driver mutation frequecies in SCAN-B TNBC
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
pacman::p_load(readxl,
               reshape2)
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
infile.3 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
infile.4 <- "./data/SCANB/5_TNBC_NatMed/driver.df.scanb.complete.csv"
infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.7 <- "./data/Parameters/TNBC_color_palette.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_mutfreqs_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_mutfreqs_TNBC.txt")
#-------------------
# storing objects
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load IDkey and correct sampleIDs -> ask Johan for key 
id.key <- loadRData(infile.6)
id.key <- id.key[c("Tumour","Specimen_id")]

# load palette
color.palette <- loadRData(infile.7)

# wgs QC
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]
qc.samples.s <- id.key$Specimen_id[match(qc.samples,id.key$Tumour)]

# driver data SCANB
driv.scanb <- loadRData(infile.3)
driv.scanb$Subtype <- "ERpHER2n_Basal"

# driver genes from tnbc groups
tnbc.samples <- loadRData(infile.1)[c("TNBC_Basal","TNBC_NonBasal")]
driv.tnbc <- read.table(infile.4,sep=",",header=TRUE)
tnbc.anno <- loadRData(infile.5)
tnbc.anno$Subtype <- sapply(tnbc.anno$External_ID_sample, function(sampleID) {
  for (sublist_name in names(tnbc.samples)) {
    if (sampleID %in% tnbc.samples[[sublist_name]]) {
      return(sublist_name)
    }
  }
  return(NA)
})

# convert sample name
driv.tnbc$SampleID <- tnbc.anno$External_ID_sample[match(driv.tnbc$Sample,tnbc.anno$PD_ID)]
driv.tnbc$Subtype <- tnbc.anno$Subtype[match(driv.tnbc$SampleID,tnbc.anno$External_ID_sample)]

tnbc.anno <- tnbc.anno[tnbc.anno$External_ID_sample %in% unlist(tnbc.samples),]
driv.tnbc <- driv.tnbc[driv.tnbc$SampleID %in% unlist(tnbc.samples),]
#View(driv.tnbc)

# get into right format
# SampleID Gene Type+Effect Subtype
driv.tnbc$variant_class <- paste0(driv.tnbc$Type,"_",driv.tnbc$Effect)
driv.tnbc <- driv.tnbc[c("SampleID","Gene","variant_class","Subtype")]
names(driv.tnbc) <- names(driv.scanb)

#- cont here; match these
table(driv.tnbc$variant_class)
table(driv.scanb$variant_class)






# exclude CopyNumber_HD & Complex_frameshift because these were not assessed in SCANB
driv.basis <- driv.basis[driv.basis$variant_class %!in% 
                              c("Complex_frameshift","CopyNumber_HD"),]

#View(driv.basis)
#View(basis.anno)
#table(driv.basis$PAM50)

# in one df
#head(driv.scanb)
#head(driv.basis)
driv.dat <- as.data.frame(rbind(driv.scanb,driv.basis))
#View(driv.dat)

# total sample counts to calc. mut freqs
pam50.counts <- table(
  driv.dat[!duplicated(driv.dat[,c("sample")]),]$PAM50)

#######################################################################
# plot and stats
#######################################################################

# genes to plot/test
genes <- c("TP53","MYC","PIK3CA","BRCA2","ESR1")
for (gene in genes) {
  #print(gene)
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  # gene data
  gene.dat <- as.data.frame(table(driv.dat[driv.dat$gene==gene,]$PAM50))
  row.names(gene.dat) <- gene.dat$Var1
  gene.dat$Var1 <- NULL
  gene.dat$Not_mutated <- pam50.counts[row.names(gene.dat)] - gene.dat$Freq
  names(gene.dat) <- c("Mutated","Not_Mutated")
  gene.mut.freqs <- (gene.dat$Mutated / pam50.counts[row.names(gene.dat)])*100
  
  # statistics
  txt.out <- append(txt.out, c("\n",capture.output(gene.mut.freqs), "\n",
                               "\n###########################################\n"))
  
  # mann whitney u tests
  luma.res <- if(sum(is.na(gene.dat[c("Basal","LumA"),]))==0) {
    fisher.test(gene.dat[c("Basal","LumA"),]) } else {"NA"}
  lumb.res <- if(sum(is.na(gene.dat[c("Basal","LumB"),]))==0) {
    fisher.test(gene.dat[c("Basal","LumB"),]) } else {"NA"}
  
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    height = gene.mut.freqs[c("LumA","LumB","Basal")], 
    names=names(color.palette), 
    col = color.palette,
    ylim=c(0,100),
    main=gene,
    ylab="Mutation frequency (%)")
  plot.parameters <- append(plot.parameters, list(plot.par))
}

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
