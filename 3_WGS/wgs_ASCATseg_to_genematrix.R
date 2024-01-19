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
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv" # move ot basal project
# output paths
#outfile.1 <- ""
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

ascat.list <- loadRData("./data/SCANB/3_WGS/raw/ASCAT_segments_list_BASEprocessing.RData")

View(head(ascat.list))
View(head(ascat.list[[1]]))

#######################################################################
# correct sample IDs
#######################################################################



#######################################################################
# get chromosome lengths
#######################################################################
# why do i need that

# chr lengths
# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = infile.3, sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$chr <- as.numeric(gsub('^.{3}','',chr.lengths$chr))

#######################################################################
# convert to gene matrix, also stored in list
#######################################################################

# convert


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

# get segment data
#cn.scanb.segments <- loadRData(scanb.segments)
ascat.list
colnames(ascat.list[[1]])

# Create empty data frames for storing results
dat.cols <- c("CNA","LOH","cnnLOH","Amp","HomDel")
res.names <- c("gene", "chr", "start", "end", c(names(ascat.list)))
for (name in dat.cols) {
  assign(paste0(name, ".df"), data.frame(matrix(nrow = 0, ncol = length(res.names), 
                                                dimnames = list(NULL, res.names))))
}

# loop over genes
i=1 # DEL LATER


pb = txtProgressBar(min = 0, max = length(genes$SYMBOL), initial = 0, style = 3)
for (i in 1:nrow(genes)) {
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
  gene.res <- t(gene.statuses)
  View(gene.res)
  
  # CONT HERE
  # store results
  
  gl.df[i,] <- c(gene.dat$SYMBOL,
                 gene.dat$chr,
                 gene.dat$start,
                 gene.dat$end,
                 unname(unlist(gene.statuses["GainLoss",])))
  
  amp.df[i,] <- c(gene.dat$SYMBOL,
                  gene.dat$chr,
                  gene.dat$start,
                  gene.dat$end,
                  unname(unlist(gene.statuses["Amp",])))
  
  close(pb)
}

cn.list <- list("gainloss"=gl.df,"amp"=amp.df)

# add center position
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% mutate_at(c("chr","start","end"), as.numeric)
  centerPos.vec <- apply(df, 1, function(row) {
    gene.length <- as.numeric(row[["end"]]) - as.numeric(row[["start"]])
    centerPos <- as.numeric(row[["start"]]) + (gene.length/2)
    return(centerPos)
  })
  df$centerPos <- centerPos.vec
  df <- df %>% dplyr::relocate(centerPos, .after=chr)
  return(df)
})

# convert to genome position
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% 
    mutate_at(c("chr","start","end"), as.numeric) %>% 
    # add new chr genome position column
    mutate(genome = 0) %>% 
    # update the genome col to fill in the actual chr positions
    rows_update(chr.lengths[c("chr","genome")]) %>% 
    # add a column with the genome position of each probe
    mutate(Genome_pos = centerPos + genome) %>% 
    relocate(c(genome,Genome_pos), .after=centerPos) %>% 
    dplyr::select(-c(genome))
  return(df)
})

# remove start and end
cn.list <- lapply(cn.list, function(df) {
  df <- df %>% 
    dplyr::select(-c(start, end))
  return(df)
})

#save(cn.list,
#     file = scanb.gene.cna)

cn.scanb.list <- loadRData(scanb.gene.cna)
