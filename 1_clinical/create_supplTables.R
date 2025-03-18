# Script: Create supplementary tables for Basal manuscript
# Author: Lennart Hohmann
# Date: 04.03.2025
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/1_clinical/src/clin_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(openxlsx, writexl,readxl)
#-------------------
# set/create output directories
output.path <- "./output/0_SupplementaryTables/"
dir.create(output.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/0_GroupSamples/TNBC_sampleIDs.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.4 <- "./data/SCANB/2_transcriptomic/processed/Metagene_scores_All.RData"
infile.5 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGSanno.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/ASCAT_segments_list_BASEprocessing.RData"
infile.7 <- "./data/SCANB/3_WGS/processed/processed_signatures.RData"
infile.8 <- "./data/SCANB/3_WGS/raw/2024_02_14_hrdetect_refsig_params_low_burden_sv_accounted_for.csv"
infile.9 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
infile.10 <- "./data/SCANB/3_WGS/raw/HRDproject_table_S1.xlsx"
infile.11 <- "./data/SCANB/4_CN/processed/CNA_genetest.RData"
infile.12 <- "./data/SCANB/4_CN/processed/CNA_genetest_TNBC.RData"
infile.13 <- "./data/SCANB/2_transcriptomic/processed/DE_result_FPKM_TNBC.RData"
infile.14 <- "./data/SCANB/2_transcriptomic/processed/DE_result_FPKM.RData"
infile.15 <- "./data/SCANB/1_clinical/raw/cst.txt"
infile.16 <- "./data/SCANB/1_clinical/raw/GSE278586_Annotations_Step1.txt"
infile.17 <- "./data/SCANB/3_WGS/raw/Project2_Basal_like_Drivers_22Jan24_ForJohan.xlsx"
infile.18 <- "./data/SCANB/3_WGS/raw/MergedAnnotations_ERp_Cohort_FailFiltered.RData"
infile.19 <- "./data/SCANB/3_WGS/raw/SCANB_ERpos_Project2.xlsx"
# output paths
outfile.1 <- paste0(output.path,"SCANB_SupplementaryTableS1.xlsx")
outfile.2 <- paste0(output.path,"SCANB_SupplementaryTableS2.xlsx")
outfile.3 <- paste0(output.path,"SCANB_SupplementaryTableS3.xlsx")
outfile.4 <- paste0(output.path,"SCANB_SupplementaryTableS4.xlsx")
outfile.5 <- paste0(output.path,"SCANB_SupplementaryTableS5.xlsx")

#######################################################################
# functions
#######################################################################

add_sheet <- function(wb, sheet_name, data) {
  addWorksheet(wb, sheet_name)   # Add a new worksheet with the specified name
  writeData(wb, sheet_name, data)  # Write the provided data frame to the worksheet
}

#######################################################################
# load data
#######################################################################

sample.ids <- append(loadRData(infile.1),
            loadRData(infile.2)[c("TNBC_Basal","TNBC_NonBasal")])

clin.dat <- loadRData(infile.3)
clin.dat <- clin.dat[clin.dat$Follow.up.cohort==TRUE,]

mg.dat <- loadRData(infile.4)

wgs.anno <- loadRData(infile.5)

ascat.segments <- do.call(rbind, loadRData(infile.6))
ascat.segments$sample <- sub("\\..*", "", ascat.segments$sample)
ascat.segments <- ascat.segments[ascat.segments$sample %in% wgs.anno$SampleID[wgs.anno$WGS==1],]

signature.dat <- loadRData(infile.7)

hrd.dat <- read.table(infile.8, sep = ",", header = TRUE)
hrd.dat$Lund.tumour.id <- sub("\\..*", "", hrd.dat$Lund.tumour.id)
hrd.dat <- hrd.dat[hrd.dat$Lund.tumour.id %in% wgs.anno$SampleID[wgs.anno$WGS==1],]
hrd.dat$HRDetect <- ifelse(hrd.dat$Probability >= 0.7,"HRD-high","HRD-low")

mut.dat <- loadRData(infile.9)
mb.dat <- read_excel(infile.10,sheet="Table_S1")
mb.dat <- mb.dat[which(mb.dat$Specimen %in% wgs.anno$SampleID[wgs.anno$WGS==1]),
                 c("Specimen","PAM50_NCN","Caveman.counts_CLPM.0_ASMD140_final",
                   "Pindel.counts_QUAL250_Repeats_10_final")]
CNA.genetest.erp <- loadRData(infile.11)
CNA.genetest.tn <- loadRData(infile.12)
CNA.genetest.all <- merge(CNA.genetest.erp, CNA.genetest.tn, by = "gene", all = TRUE)

DE.erp <- loadRData(infile.14)
DE.tn <- loadRData(infile.13)
DE.erp <- DE.erp[rownames(DE.erp) != "", ]
DE.tn <- DE.tn[rownames(DE.tn) != "", ]
DE.erp$gene <- rownames(DE.erp)
DE.tn$gene <- rownames(DE.tn)
DE <- merge(DE.erp, DE.tn, by = "gene", all = TRUE)
rownames(DE) <- DE$gene
#DE$gene <- NULL 

# ordered vectors of DEGs (ordered based on logFC for later GSEA)
fpkm.degs.luma.df <- DE[DE$Basal.LumA.padj <= 0.05 & abs(DE$Basal.LumA.logFC) > log2(2),]
fpkm.degs.luma <- rownames(fpkm.degs.luma.df)[order(fpkm.degs.luma.df$Basal.LumA.logFC,decreasing=TRUE)]
fpkm.degs.lumb.df <- DE[DE$Basal.LumB.padj <= 0.05 & abs(DE$Basal.LumB.logFC) > log2(2),]
fpkm.degs.lumb <- rownames(fpkm.degs.lumb.df)[order(fpkm.degs.lumb.df$Basal.LumB.logFC,decreasing=TRUE)]
core.degs <- intersect(fpkm.degs.luma,fpkm.degs.lumb)

gl.freqs <- loadRData("./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData")
gl.freqs.tn <- loadRData("./data/SCANB/4_CN/processed/CNA_GLFreqs_TNBC.RData")
gl.freqs.tn <- gl.freqs.tn[,!(names(gl.freqs.tn) %in% c("chr","start","end"))]

gl.freqs <- merge(gl.freqs,gl.freqs.tn,by="gene",all=TRUE)

genes.AG <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$LumA.Gain.padj <= 0.05, ]$gene, "gene"]

genes.AL <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$LumA.Loss.padj <= 0.05, ]$gene,  "gene"]

genes.BG <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$LumB.Gain.padj <= 0.05, ]$gene,  "gene"]

genes.BL <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$LumB.Loss.padj <= 0.05, ]$gene,  "gene"]

genes.bG <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$TNBC.Basal.Gain.padj <= 0.05, ]$gene,  "gene"]
genes.bL <- gl.freqs[gl.freqs$gene %in% 
                       CNA.genetest.all[CNA.genetest.all$TNBC.Basal.Loss.padj <= 0.05, ]$gene,  "gene"]
genes.nbG <- gl.freqs[gl.freqs$gene %in% 
                        CNA.genetest.all[CNA.genetest.all$TNBC.NonBasal.Gain.padj <= 0.05, ]$gene,  "gene"]
genes.nbL <- gl.freqs[gl.freqs$gene %in% 
                        CNA.genetest.all[CNA.genetest.all$TNBC.NonBasal.Loss.padj <= 0.05, ]$gene,  "gene"]

#######################################################################
# Table 1: Full data for all Basal cases + GSM ids
#######################################################################
erperc.dat <- read.table(infile.15,sep="\t",header=TRUE)
erperc.dat <- erperc.dat[erperc.dat$rba %in% clin.dat$GEX.assay,]
erperc.dat <- erperc.dat[c("Specimen","INCA2_op_pad_erproc")]
names(erperc.dat) <- c("Sample","INCA2_op_pad_erproc")
clin.dat <- merge(clin.dat,erperc.dat,by="Sample")
clin.dat <- clin.dat[c("Sample","ER","PR","HER2","NCN.PAM50","LN",
                         "NHG","Age","Size.mm","INCA2_op_pad_erproc",
                         "TreatGroup","InvCa.type","DRFI","DRFIbin",
                         "OS","OSbin")] 
clin.dat <- clin.dat[clin.dat$Sample %in% sample.ids$ERpHER2n_Basal,]
mg.dat <- as.data.frame(t(mg.dat))
colnames(mg.dat) <- paste("metagene.", colnames(mg.dat), sep = "")
mg.dat$Sample <- rownames(mg.dat)
clin.dat <- merge(clin.dat,mg.dat,by="Sample")
wgs.anno <- wgs.anno[c("SampleID","WGS")]
names(wgs.anno) <- c("Sample","WGS")
clin.dat <- merge(clin.dat,wgs.anno,by="Sample")

gsm.ids <- read.table(infile.16, sep = "\t",header = TRUE)
gsm.ids$Sample <- gsub(".*(S\\d+).*", "\\1", gsm.ids$Title)
gsm.ids <- gsm.ids[c("Sample","GSMid")]
gsm.ids <- gsm.ids[gsm.ids$Sample %in% wgs.anno$Sample[wgs.anno$WGS==1],]
clin.dat <- merge(clin.dat,gsm.ids,by="Sample", all=TRUE)

wb <- createWorkbook()
add_sheet(wb, "ERpHER2n-Basal_annotations", clin.dat)
# Save the workbook to an Excel file
saveWorkbook(wb, outfile.1, overwrite = TRUE)

#######################################################################
# Table 3: WGS data for all Basal cases
#######################################################################

hrd.dat <- hrd.dat[c("Lund.tumour.id","intercept","del.mh.prop","SNV3",
                  "SV3","SV5","hrd","SNV8","Probability","HRDetect")] 
#View(mut.dat)
id.key <- loadRData(infile.18)
wgs.sum <- read_excel(infile.19, sheet = "Summary")
wgs.sum$`Final QC` <- toupper(wgs.sum$`Final QC`)
wgs.sum <- wgs.sum[wgs.sum$`Final QC` %in% c("PASS","AMBER"),]
wgs.sum$Sample <- id.key$Specimen_id[match(wgs.sum$Tumour,id.key$Tumour)]
wgs.sum <- wgs.sum[c("Sample","Batch","Final QC","Tumour","Tumour duplicate rate",
                     "Tumour mean depth","Normal","Normal duplicate rate",
                     "Normal mean depth","Final_Ploidy","Final_Aberrant cell fraction",
                     "Caveman counts (CLPM=0,ASMD>=140)_final",
                     "Pindel counts (QUAL>=250,Repeats>10)_final",
                     "BRASS Counts Assembly score >0_final")]
#View(wgs.sum)
mut.1 <- read_excel(infile.17,sheet = "PindelCoding")
mut.2 <- read_excel(infile.17,sheet = "PindelDrivers")
mut.3 <- read_excel(infile.17,sheet = "CavemanCoding")
mut.4 <- read_excel(infile.17,sheet = "CavemanDrivers")
mut.5 <- read_excel(infile.17,sheet = "BRASS_Coding")
mut.6 <- read_excel(infile.17,sheet = "BRASS_drivers")
mut.list <- lapply(list(mut.1,mut.2,mut.3,mut.4,mut.5,mut.6), function(df) {
  names(df) <- toupper(names(df))
  df <- df[df$SAMPLE %in% wgs.sum$Tumour,]
  df$SAMPLE <- id.key$Specimen_id[match(df$SAMPLE,id.key$Tumour)]
  return(df)
})

wb <- createWorkbook()
add_sheet(wb, "WGS_summary", wgs.sum) # check
add_sheet(wb, "WGS_signatures", signature.dat) # check
add_sheet(wb, "HRDetect", hrd.dat) # check
mut.sheet.names <- c("PindelCoding","PindelDrivers","CavemanCoding","CavemanDrivers","BRASSCoding","BRASSDrivers")
for(i in 1:length(mut.list)) {
  add_sheet(wb, mut.sheet.names[i], mut.list[[i]])
}
#add_sheet(wb, "AllDriver_mutations", mut.dat)

#add_sheet(wb, "Total_mutation_counts", mb.dat)

# Save the workbook to an Excel file
saveWorkbook(wb, outfile.3, overwrite = TRUE)

# get complete somatic wgs data
infile.8 <- "./data/SCANB/3_WGS/raw/2024_02_14_hrdetect_refsig_params_low_burden_sv_accounted_for.csv"
infile.9 <- "./data/SCANB/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
hrd.dat


#######################################################################
# Table 2: DE results, core set
#######################################################################

wb <- createWorkbook()
add_sheet(wb, "DE_results", DE)
add_sheet(wb, "coreDEGs", core.degs)
# Save the workbook to an Excel file
saveWorkbook(wb, outfile.2, overwrite = TRUE)

#######################################################################
# Table 4: CNA results, freqs, genetests
#######################################################################
wb <- createWorkbook()

# Find the maximum length among all vectors
max_length <- max(length(genes.AG), length(genes.AL), length(genes.BG), length(genes.BL), 
                  length(genes.bG), length(genes.bL), length(genes.nbG), length(genes.nbL))

# Create the data frame, padding shorter vectors with NA
genes.vs.basal <- data.frame(
  genes.LumA.gain = c(genes.AG, rep(NA, max_length - length(genes.AG))),
  genes.LumA.loss = c(genes.AL, rep(NA, max_length - length(genes.AL))),
  genes.LumB.gain = c(genes.BG, rep(NA, max_length - length(genes.BG))),
  genes.LumB.loss = c(genes.BL, rep(NA, max_length - length(genes.BL))),
  genes.TNBCBasal.gain = c(genes.bG, rep(NA, max_length - length(genes.bG))),
  genes.TNBCBasal.loss = c(genes.bL, rep(NA, max_length - length(genes.bL))),
  genes.TNBCNonBasal.gain = c(genes.nbG, rep(NA, max_length - length(genes.nbG))),
  genes.TNBCNonBasal.loss = c(genes.nbL, rep(NA, max_length - length(genes.nbL)))
)

add_sheet(wb, "CNA_frequencies", gl.freqs)
add_sheet(wb, "CNA_genetests_vs_ERpHER2n-Basal", CNA.genetest.all)
add_sheet(wb, "Signif_genes_vs_ERpHER2n-Basal", genes.vs.basal)

# Save the workbook to an Excel file
saveWorkbook(wb, outfile.4, overwrite = TRUE)

#######################################################################
# Table 5: DNA methylation results ?
#######################################################################
