# Script: Driver mutation frequecies in SCAN-B all
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
infile.3 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nAll.RData"
infile.4 <- "./data/SCANB/5_TNBC_NatMed/driver.df.scanb.complete.csv"
infile.5 <- "./data/SCANB/5_TNBC_NatMed/Updated_merged_annotations_n235_WGS_MethylationCohort.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.7 <- "./data/Parameters/TNBC_color_palette.RData"
infile.8 <- "./data/Parameters/color_palette.RData"

# output paths
plot.file <- paste0(output.path,cohort,"_mutfreqs_All.pdf")
txt.file <- paste0(output.path,cohort,"_mutfreqs_All.txt")
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
color.palette.1 <- loadRData(infile.7)[c("TNBC_NonBasal",
                                         "TNBC_Basal","ERpHER2n_Basal")]
color.palette.2 <- loadRData(infile.8)[c("LumB", "LumA")]
color.palette <- c(color.palette.1, color.palette.2)
names(color.palette)[names(color.palette) == "LumA"] <- "ERpHER2n_LumA"
names(color.palette)[names(color.palette) == "LumB"] <- "ERpHER2n_LumB"

# wgs QC
wgs.sum <- read_excel(infile.2, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
qc.samples <- wgs.sum$Tumour[-grep("FAIL", wgs.sum$`Final QC`)]
qc.samples.s <- id.key$Specimen_id[match(qc.samples,id.key$Tumour)]

# driver data SCANB
driv.scanb <- loadRData(infile.3)
driv.scanb$Subtype <- paste0("ERpHER2n_",driv.scanb$PAM50)
driv.scanb$PAM50 <- NULL

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

# get into right format to match scanb basal

driv.tnbc$variant_class <- paste0(driv.tnbc$Type,"_",driv.tnbc$Effect)
# not included in the basal drivers
driv.tnbc <- driv.tnbc[driv.tnbc$variant_class != "loss_NA", ]
# indel framesift, indel inframe
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class %in% c("Del_frameshift","Ins_frameshift"),
                                  "indel_frameshift",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class %in% c("Del_inframe","Ins_inframe"),
                                  "indel_inframe",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "amp_NA",
                                  "CN_amplification",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "Sub_missense",
                                  "point_missense",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "Sub_nonsense",
                                  "point_nonsense",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.deletion_NA",
                                  "rearr_deletion",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.inversion_NA",
                                  "rearr_inversion",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.tandem-duplication_NA",
                                  "rearr_tandem-duplication",driv.tnbc$variant_class)
driv.tnbc$variant_class <- ifelse(driv.tnbc$variant_class == "r.translocation_NA",
                                  "rearr_translocation",driv.tnbc$variant_class)

driv.tnbc <- driv.tnbc[c("SampleID","Gene","variant_class","Subtype")]
names(driv.tnbc) <- names(driv.scanb)

#table(driv.tnbc$variant_class)
#table(driv.scanb$variant_class)
#head(driv.tnbc)
#head(driv.scanb)

driv.dat <- as.data.frame(rbind(driv.scanb,driv.tnbc))
#View(driv.dat)

# total sample counts to calc. mut freqs, not all tnbc samples here, some have no driv muts?
#subtype.counts <- table(
#   driv.dat[!duplicated(driv.dat[,c("sample")]),]$Subtype)

subtype.counts <- c("ERpHER2n_Basal"=16,"ERpHER2n_LumA"=73,
                      "ERpHER2n_LumB"=105,"TNBC_Basal"=184,
                      "TNBC_NonBasal"=44)

#######################################################################
# plot and stats
#######################################################################

# genes to plot/test
genes <- c("TP53","MYC","PIK3CA","BRCA2","ESR1","RB1","PTEN")
for (gene in genes) {
  #print(gene)
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  
  # gene data
  gene.dat <- as.data.frame(table(driv.dat[driv.dat$gene==gene,]$Subtype))
  row.names(gene.dat) <- gene.dat$Var1
  gene.dat$Var1 <- NULL
  gene.dat$Not_mutated <- subtype.counts[row.names(gene.dat)] - gene.dat$Freq
  names(gene.dat) <- c("Mutated","Not_Mutated")
  gene.mut.freqs <- (gene.dat$Mutated / subtype.counts[row.names(gene.dat)])*100
  
  # statistics
  txt.out <- append(txt.out, c("\n",capture.output(gene.mut.freqs), "\n",
                               "\n###########################################\n"))
  
  # mann whitney u tests
  nonbas.res <- if(sum(is.na(gene.dat[c("ERpHER2n_Basal","TNBC_NonBasal"),]))==0) {
    fisher.test(gene.dat[c("ERpHER2n_Basal","TNBC_NonBasal"),]) } else {"NA"}
  bas.res <- if(sum(is.na(gene.dat[c("ERpHER2n_Basal","TNBC_Basal"),]))==0) {
    fisher.test(gene.dat[c("ERpHER2n_Basal","TNBC_Basal"),]) } else {"NA"}
  luma.res <- if(sum(is.na(gene.dat[c("ERpHER2n_Basal","ERpHER2n_LumA"),]))==0) {
    fisher.test(gene.dat[c("ERpHER2n_Basal","ERpHER2n_LumA"),]) } else {"NA"}
  lumb.res <- if(sum(is.na(gene.dat[c("ERpHER2n_Basal","ERpHER2n_LumB"),]))==0) {
    fisher.test(gene.dat[c("ERpHER2n_Basal","ERpHER2n_LumB"),]) } else {"NA"}
  txt.out <- append(txt.out, c(capture.output(nonbas.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(bas.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  
  # plot
  plot.par <- list(
    height = gene.mut.freqs[c("ERpHER2n_LumA","ERpHER2n_LumB","ERpHER2n_Basal","TNBC_Basal","TNBC_NonBasal")], 
    names=c("ERpHER2n_LumA","ERpHER2n_LumB","ERpHER2n_Basal","TNBC_Basal","TNBC_NonBasal"), 
    col = color.palette[c("ERpHER2n_LumA","ERpHER2n_LumB","ERpHER2n_Basal","TNBC_Basal","TNBC_NonBasal")],
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
  axis(3,at=bp,labels=subtype.counts[plot.parameters[[i]]$names]) # change to all samples?
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
