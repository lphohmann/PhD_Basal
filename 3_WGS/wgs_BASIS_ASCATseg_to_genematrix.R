# Script: Processing ASCAT data in SCAN-B
# Author: Lennart Hohmann
# Date: 18.01.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "BASIS"
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
data.path <- "./data/BASIS/3_WGS/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/BASIS/3_WGS/raw/ASCAT_CEL_Total_EasySegments.RData"
infile.2 <- "./data/BASIS/3_WGS/raw/ASCAT_InSilico_EasySegments.RData"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
# output paths
outfile.1 <- "./data/BASIS/3_WGS/processed/ASCAT_genelevel.RData"
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

# get segment data
cn.basis.segments.1 <- loadRData(infile.1)
cn.basis.segments.2 <- loadRData(infile.2)
cn.basis.segments <- rbind(cn.basis.segments.1,cn.basis.segments.2)
#View(cn.basis.segments)
#######################################################################
# correct sample IDs and data structure
#######################################################################

# get relevant sample IDs
basis.anno <- loadRData(infile.3)
basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & 
                           basis.anno$PAM50_AIMS %in% c("LumA","LumB"),
                         c("sample_name","PAM50_AIMS")]

# select relevant samples
cn.basis.segments$SampleID <- gsub("b","",gsub("a.*","",cn.basis.segments$SampleID))
cn.basis.segments <- cn.basis.segments[cn.basis.segments$SampleID %in% basis.anno$sample_name,]

# organize in list of dfs (1 per sample)
cn.basis.segments <- with(cn.basis.segments, split(
  cn.basis.segments, list(SampleID)))

################################################################################
# prepare the gene annotation data to which the segments are mapped
################################################################################

# get gene positions & convert to hgnc symbols
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) # meta = entrez ID)
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

#colnames(cn.basis.segments[[1]])
dat.cols <- c("num.mark","seg.mean","CNA","LOH","cnnLOH")
res.list <- list()

# loop over genes
pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) { #nrow(genes)) { 
  setTxtProgressBar(pb,i)
  
  gene.dat <- genes[i,]
  print(gene.dat$SYMBOL)
  
  # get the gainloss and amp state of that gene for all samples
  gene.statuses <- sapply(cn.basis.segments, function(segment.df) {
    #segment.df <- cn.basis.segments[[1]] # DEL LATER
    
    numeric_columns <- c('chr','start','stop', dat.cols)
    segment.df[, numeric_columns] <- lapply(segment.df[, numeric_columns], as.numeric)
    
    # get the cn status of the segment the gene is on
    # only relevant segments
    chr.segments <- segment.df[segment.df$chr==gene.dat$chr,]
    
    # check gene/segment overlap
    segments.query <- with(chr.segments, IRanges(start, stop))
    gene.subject <- with(gene.dat, IRanges(start, end))
    
    # check 
    chr.segments$overlap = countOverlaps(segments.query, gene.subject) != 0 # calculating overlaps
    if (sum(chr.segments$overlap)==0) { # no segment covers the gene region
      main.segment <- setNames(rep(NA, length(dat.cols)), dat.cols)
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

#save(res.list,file=outfile.1)
#######################################################################
# format results to get 1 df per sample
#######################################################################

# Create a list to store individual data frames for each sampleID
sample.res <- list()

# Loop through each sampleID and create a data frame containing corresponding rows
for (sampleID in names(cn.basis.segments)) {
  sample.df <- do.call(rbind, lapply(res.list, function(df) {subset(df, sample == sampleID)}))
  rownames(sample.df) <- NULL
  sample.res[[as.character(sampleID)]] <- sample.df[c("sample", 
                                                      colnames(sample.df)[colnames(sample.df) != "sample"])]
}

# save
save(sample.res, file= outfile.1)

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)
# takes approx 9.749125 hours
