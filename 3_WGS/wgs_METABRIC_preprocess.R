# Script: Preprocess METABRIC mut data into correct format
# Author: Lennart Hohmann
# Date: 11.11.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "METABRIC"
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
#-------------------
# set/create output directories
# for data
data.path <- "./data/METABRIC/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.0 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.1 <- "./data/METABRIC/3_WGS/raw/data_mutations_extended.txt"
infile.2 <- "./data/METABRIC/4_CN/raw/data_CNA.txt"
infile.3 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
infile.4 <- "./data/SCANB/3_WGS/processed/driver_HugoSymbols.RData"
# output paths
outfile.1 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
outfile.2 <- "./data/METABRIC/3_WGS/processed/drivermutations_All.RData"
outfile.3 <- "./data/METABRIC/4_CN/processed/CNA_genelevel_all.RData"
outfile.4 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nAll.RData"
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# 
sample.IDs.erp <- loadRData(infile.0)

# load annotation data and select 
anno <- loadRData(infile.3)
anno <- anno[!is.na(anno$METABRIC_ID),]

# driver gene symbols
driver.genes <- loadRData(infile.4)

#######################################################################
# CN data
#######################################################################

gene.anno <- loadRData("./data/SCANB/4_CN/processed/CNA_genelevel_all.RData")[["ascat.df.all"]] #"gene"    "chr"     "start"   "end"     "S000006" /(-1,0,1)
gene.anno <- gene.anno[c("gene","chr","start","end")]
#head(gene.anno)

# load and prep gene expr. data
cn.data <- read.delim(infile.2,header = TRUE, sep = "\t", dec = ".")
colnames(cn.data)[2:ncol(cn.data)] <- gsub("\\.", "-", colnames(cn.data)[2:ncol(cn.data)])
cn.data$Entrez_Gene_Id <- NULL
cn.data <- cn.data[!is.na(cn.data$Hugo_Symbol), ]
cn.data <- cn.data[!duplicated(cn.data$Hugo_Symbol), ]
row.names(cn.data) <- cn.data$Hugo_Symbol
cn.data$Hugo_Symbol <- NULL
cn.data <- cn.data[, colnames(cn.data) %in% anno$METABRIC_ID]
cn.data[] <- lapply(cn.data, as.numeric)
cn.data <- na.omit(cn.data)

# save genes and samples that have amplifications for genes
results <- list()
for (sample in names(cn.data)) {
  matches <- which(cn.data[[sample]] == 2)  # Get rows where value is 2
  if (length(matches) > 0) {
    results[[sample]] <- data.frame(
      sample = sample,
      gene = rownames(cn.data)[matches],
      Value = cn.data[[sample]][matches]
    )
  }
}
amp.genes <- do.call(rbind, results)
amp.genes$variant_class <- "CN_amp"
amp.genes$PAM50 <- anno$PAM50[match(amp.genes$sample, anno$METABRIC_ID)]
amp.genes$Value <- NULL
row.names(amp.genes) <- NULL
#table(unlist(as.data.frame(cn.data))) # -2       -1        0        1        2 
# get same format as scanb
cn.data[cn.data < -1] <- -1
cn.data[cn.data > 1] <- 1
#table(unlist(as.data.frame(cn.data))) # -1        0        1
#head(cn.data)
cn.data$gene <- row.names(cn.data)

cn.data <- merge(cn.data, gene.anno, by="gene")
cn.data <- cn.data[, c(c("gene","chr","start","end"), 
                       setdiff(names(cn.data), 
                               c("gene","chr","start","end")))]
save(cn.data, file=outfile.3)
#sampleIDs.erp <- loadRData(infile.0)
#cn.data.erp <- cn.data[,names(cn.data) %in% c(unname(unlist(sampleIDs.erp))]

#######################################################################
# mutational data
#######################################################################

# mut data
mut.data <- as.data.frame(read.delim(infile.1, header = FALSE, sep = "\t", dec = "."))
colnames(mut.data) <- mut.data[2,]
mut.data <- mut.data[-c(1, 2), ]
mut.data <- mut.data[c("Hugo_Symbol", 
                       "Variant_Classification",
                       "Tumor_Sample_Barcode")]
mut.data <- mut.data[mut.data$Tumor_Sample_Barcode %in% anno$METABRIC_ID,]
names(mut.data) <- c("gene","variant_class","sample")
mut.data <- mut.data[c("sample","gene","variant_class")]
mut.data$PAM50 <- anno$PAM50[match(mut.data$sample, anno$METABRIC_ID)]
# filter genes to be same list as in scanb
mut.data <- rbind(mut.data,amp.genes)
mut.data <- mut.data[mut.data$gene %in% driver.genes,] # for mut burden change this step
#head(mut.data)
save(mut.data, file = outfile.2)

mut.data.basal <- mut.data[mut.data$sample %in% sample.IDs.erp[["ERpHER2n_Basal"]],]
#length(unique(mut.data.basal$sample)) #40

save(mut.data.basal, file = outfile.1)

mut.data.ERpAll <- mut.data[mut.data$sample %in% unname(unlist(sample.IDs.erp)),]
save(mut.data.ERpAll, file = outfile.4)
