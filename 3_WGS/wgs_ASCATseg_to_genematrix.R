# Script: Processing ASCAT data in SCAN-B
# Author: Lennart Hohmann
# Date: 18.01.2024
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
pacman::p_load(ggplot2,
               tidyverse,
               reshape2,
               #library(data.table)
               purrr,
               readxl,
               IRanges,
               GenomicFeatures,
               TxDb.Hsapiens.UCSC.hg38.knownGene,
               TxDb.Hsapiens.UCSC.hg19.knownGene,
               org.Hs.eg.db,
               Repitools)
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
infile.1 <- "./data/SCANB/3_WGS/raw/ASCAT_segments_list_BASEprocessing.RData"
# output paths
outfile.1 <- "./data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()
#-------------------
start.time <- Sys.time()

#######################################################################
# load data
#######################################################################

ascat.list <- loadRData(infile.1)

#######################################################################
# correct sample IDs
#######################################################################

# scanb
# # ID key file
# id.key <- as.data.frame(read_excel("./data/SCANB/3_genomic/raw/HER2_enriched_June23_ForJohan.xlsx", sheet = "Samples")) %>% 
#   dplyr::select(c("Sample","Tumour")) %>% dplyr::rename(sample=Sample)
# 
# # correct ids now
# res$sampleID <- id.key$sample[match(res$sample,id.key$Tumour)]
# res$sample <- NULL
# res$PAM50 <- rep("HER2E", length(sample))
# res$N_mut <- as.numeric(res$caveman_count) + as.numeric(res$pindel_count)
# scanb.muts <- res

################################################################################
# prepare the gene annotation data to which the segments are mapped
################################################################################

# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene) # meta = entrez ID)
ENTREZID2SYMBOL <- biomaRt::select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL

genes <- annoGR2DF(genes)
genes <- genes[!is.na(genes$SYMBOL), ]
genes$chr <- gsub("^chr", "", genes$chr)
genes[genes$chr=="X",]$chr <- 23
genes <- genes[genes$chr %in% 1:23, ]
genes <- genes[, !(names(genes) %in% c("gene_id", "strand", "width"))]
numeric_columns <- c('chr', 'start', 'end')
genes[, numeric_columns] <- lapply(genes[, numeric_columns], as.numeric)

################################################################################
# convert to gene matrix (1 row per gene), also stored in list
################################################################################

#colnames(ascat.list[[1]])
dat.cols <- c("CNA","LOH","cnnLOH","Amp","HomDel")
res.list <- list()

# loop over genes
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) { #nrow(genes)
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]
  
  # get the gainloss and amp state of that gene for all samples
  gene.statuses <- sapply(ascat.list, function(segment.df) {
    #segment.df <- ascat.list[[1]] # DEL LATER
    
    segment.df[segment.df$chr=="X",]$chr <- 23
    numeric_columns <- c('chr','startpos','endpos', dat.cols)
    segment.df[, numeric_columns] <- lapply(segment.df[, numeric_columns], as.numeric)
    
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    # check gene/segment overlap
    segments.query <- with(chr.segments, IRanges(startpos, endpos))
    gene.subject <- with(gene.dat, IRanges(start, end))
    
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0 # calculating overlaps
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      main.segment <- c("CNA"=NA,"LOH"=NA,"cnnLOH"=NA,"Amp"=NA,"HomDel"=NA)
    } else {
      hits <- findOverlaps(segments.query, gene.subject)
      overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
      percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
      chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
      chr.segments$percentOverlap <- percentOverlap
      # only select the segment with the highest overlap proportion
      main.segment <- unlist(chr.segments[which.max(chr.segments$percentOverlap),])
    }
    # gene data for that sample
    sample.gene.dat <- c("gene"=gene.dat$SYMBOL, "chr"=gene.dat$chr, "start"=gene.dat$start, 
      "end"=gene.dat$end, main.segment[dat.cols])
    return(sample.gene.dat)
  })
  
  #
  gene.res <- as.data.frame(t(gene.statuses))
  gene.res$sample <- rownames(gene.res)
  rownames(gene.res) <- NULL
  # store results
  res.list[[gene.dat$SYMBOL]] <- gene.res
  
  close(pb)
}


#######################################################################
# format results to get 1 df per sample
#######################################################################

# Create a list to store individual data frames for each sampleID
sample.res <- list()

# Loop through each sampleID and create a data frame containing corresponding rows
for (sampleID in names(ascat.list)) {
  sample.df <- do.call(rbind, lapply(res.list, function(df) {subset(df, sample == sampleID)}))
  rownames(sample.df) <- NULL
  sample.res[[as.character(sampleID)]] <- sample.df[c("sample", 
                                                      colnames(sample.df)[colnames(sample.df) != "sample"])]
}

# save
save(sample.res, file= outfile.1)

#x <- loadRData("./data/SCANB/3_WGS/processed/ASCAT_genelevel.RData")
#View(x[[1]])

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
# takes approx 9.749125 hours
