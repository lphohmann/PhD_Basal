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
infile.5 <- "./data/METABRIC/0_GroupSamples/TNBC_sampleIDs.RData"
# output paths
outfile.1 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nBasal.RData"
outfile.2 <- "./data/METABRIC/3_WGS/processed/drivermutations_All.RData"
outfile.3 <- "./data/METABRIC/4_CN/processed/CNA_genelevel_all.RData"
outfile.4 <- "./data/METABRIC/3_WGS/processed/drivermutations_ERpHER2nAll.RData"
outfile.5 <- "./data/METABRIC/4_CN/processed/CNA_GLFreqs_ERpHER2nAll.RData"
outfile.6 <- "./data/METABRIC/4_CN/processed/CNA_GLFreqs_TNBC.RData"
outfile.7 <- "./data/METABRIC/3_WGS/processed/drivermutations_TNBC.RData"

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
sample.IDs.tnbc <- loadRData(infile.5)

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

# count sample n for mut freq calc.
head(amp.genes)
unique_ids <- amp.genes[!duplicated(amp.genes$sample), ]
unique_ids.erp <- unique_ids[unique_ids$sample%in% unname(unlist(sample.IDs.erp)),]
table(unique_ids.erp$PAM50) #Basal 45  LumA 614  LumB 397 
unique_ids.tnbc <- unique_ids[unique_ids$sample%in% unname(unlist(sample.IDs.tnbc)),]
table(unique_ids.tnbc$PAM50) # check in sample file, here no LumBs?

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

# calc gain loss freqs


# Function to calculate frequency of gain and loss for a given set of sample IDs
calculate_freqs <- function(data, sample_ids) {
  # Subset the data to include only columns for the provided sample IDs
  data_subset <- data[, sample_ids, drop = FALSE]
  # Calculate frequency of losses (-1) and gains (1) per gene
  freqloss <- apply(data_subset, 1, function(x) (length(which(x == -1)) / length(sample_ids)) * -100)
  freqgain <- apply(data_subset, 1, function(x) (length(which(x == 1)) / length(sample_ids)) * 100)
  return(list(freqloss = freqloss, freqgain = freqgain))
}

# data frame to store results
CNA.freqs <- cn.data[, c("gene", "chr", "start", "end")]
# Iterate over each subtype in sample.IDs.erp to calculate and add frequencies to CNA.freqs
for (subtype in names(sample.IDs.erp)) {
  sample_ids <- sample.IDs.erp[[subtype]]  # Get sample IDs for the current subtype
  # Calculate frequencies
  freqs <- calculate_freqs(cn.data[, 5:ncol(cn.data)], sample_ids)  
  # Add frequency columns to CNA.freqs with appropriate naming
  CNA.freqs[[paste0("freqloss.", subtype)]] <- freqs$freqloss
  CNA.freqs[[paste0("freqgain.", subtype)]] <- freqs$freqgain
}

# View the result
save(CNA.freqs, file=outfile.5)

# same for tnbc

# data frame to store results
CNA.freqs <- cn.data[, c("gene", "chr", "start", "end")]
# Iterate over each subtype in sample.IDs.tnbc to calculate and add frequencies to CNA.freqs
for (subtype in names(sample.IDs.tnbc)) {
  sample_ids <- sample.IDs.tnbc[[subtype]]  # Get sample IDs for the current subtype
  # Calculate frequencies
  freqs <- calculate_freqs(cn.data[, 5:ncol(cn.data)], sample_ids)  
  # Add frequency columns to CNA.freqs with appropriate naming
  CNA.freqs[[paste0("freqloss.", subtype)]] <- freqs$freqloss
  CNA.freqs[[paste0("freqgain.", subtype)]] <- freqs$freqgain
}

# View the result
save(CNA.freqs, file=outfile.6)

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
# filter amp genes and then rbind

mut.data <- rbind(mut.data,amp.genes)
mut.data <- mut.data[mut.data$gene %in% driver.genes,] #  for mut burden change this step
#head(mut.data)
save(mut.data, file = outfile.2)

mut.data.basal <- mut.data[mut.data$sample %in% sample.IDs.erp[["ERpHER2n_Basal"]],]
#length(unique(mut.data.basal$sample)) #40

save(mut.data.basal, file = outfile.1)

mut.data.ERpAll <- mut.data[mut.data$sample %in% unname(unlist(sample.IDs.erp)),]
save(mut.data.ERpAll, file = outfile.4)


mut.data.TNBC <- mut.data[mut.data$sample %in% unname(unlist(sample.IDs.tnbc)),]
save(mut.data.TNBC, file = outfile.7)
