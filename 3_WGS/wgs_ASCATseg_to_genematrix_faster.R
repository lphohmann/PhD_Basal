#!/usr/bin/env Rscript
# Script: Processing ASCAT segments to genelevel; Hg38
# Author: Lennart Hohmann
# Date: 18.01.2024
#-------------------
#nohup ./ASCATseg_to_genematrix_faster.R > output.log 2>&1 &
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("/Users/le7524ho/Documents/Johan_ASCATgenelevel")
#-------------------
# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(purrr, readxl, IRanges, GenomicFeatures, 
               TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, Repitools, 
               future.apply, data.table)
#-------------------
# Set up parallel processing
plan(multisession, workers = parallel::detectCores() - 1)
#-------------------
# set/create output directories
output.path <- "./output/"
dir.create(output.path)
#-------------------
# input paths
infile.1 <- "./data/Collected_ASCATsegmentlist_WGS_MOtrain_1076_samples.RData"
# output paths
outfile.1 <- paste0(output.path,"ASCAT_genelevel.RData")
#-------------------
start.time <- Sys.time()

#######################################################################
# functions
#######################################################################

# function: loads RData data file and allows to assign it directly to variable
loadRData <- function(file.path) {
  load(file.path)
  get(ls()[ls() != "file.path"])
}

#######################################################################
# load data
#######################################################################

# Load ASCAT segment data
ascat.list <- loadRData(infile.1)

# Convert all ASCAT list entries to data.table (faster operations)
ascat.list <- lapply(ascat.list, as.data.table)

# Prepare gene annotation
genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
ENTREZID2SYMBOL <- biomaRt::select(org.Hs.eg.db, 
                                   mcols(genes)$gene_id, 
                                   c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL

genes <- annoGR2DF(genes)
genes <- genes[!is.na(genes$SYMBOL), ]
genes <- genes[genes$chr %in% c(paste0("chr", 1:22),"chrX"), ]
genes$chr <- droplevels(genes$chr)
genes$chr <- gsub("^chr", "", genes$chr)
genes$chr[genes$chr == "X"] <- "23"
genes$chr <- as.numeric(genes$chr)
genes <- genes[, !names(genes) %in% c("gene_id", "strand", "width")]
genes[, c('chr', 'start', 'end')] <- lapply(genes[, c('chr', 'start', 'end')], as.numeric)
#View(genes)



# smaller set to test
# set.seed(123)
# genes <- do.call(rbind, lapply(1:23, function(ch)
#   genes[genes$chr == ch, ][sample(min(10, sum(genes$chr == ch))), ]))

################################################################################
# convert to gene matrix (1 row per gene), also stored in list
################################################################################

# create IRanges object that stores gene start and end positions
gene_ranges <- IRanges(start = genes$start, end = genes$end)

# cols extracted from ASCAT segment data
# unfortunately this doent work with future_lapply, some scope issue
#dat.cols <- c("nMajor", "nMinor", "nTot", "CNA", "LOH", "cnnLOH", "Amp", "HomDel")

#tests
#i <- which(genes$SYMBOL=="MIR190B")
#segment.df <- ascat.list[["S002940.l.d.a.lib.g.a.cn"]] #.cn

# Function to process a single gene by mapping it to ASCAT segments across samples
process_gene <- function(i) {
  # data for the current gene
  gene.dat <- genes[i, ]
  gene.subject <- gene_ranges[i]  # genomic range for this gene
  print(paste0(gene.dat$SYMBOL," --- Gene: ",i,"/",nrow(genes)))
  
  # Process all samples in parallel using future_lapply()
  gene.statuses <- future_lapply(ascat.list, function(segment.df) {
    
    #print(segment.df$sample[1])
    # Convert chromosome "X" to numeric value 23 for consistency
    segment.df[segment.df$chr == "X", "chr"] <- 23
    
    # Convert key columns to numeric type to ensure consistency in calculations
    segment.df$chr <- as.numeric(segment.df$chr) 
    
    # Filter segments that are on the same chromosome as the current gene
    chr.segments <- segment.df[segment.df$chr == gene.dat$chr, ]
    
    # Create IRanges object for the ASCAT segments
    segments.query <- IRanges(start = chr.segments$startpos, end = chr.segments$endpos)
    
    # Find overlapping segments between ASCAT segments and the gene's genomic region
    hits <- findOverlaps(segments.query, gene.subject)
    
    chr.segments$overlap <- 0  # Initialize with 0 (no overlap)
    chr.segments$overlap[queryHits(hits)] <- 1  # Set overlap to 1 for segments that overlap with the gene
    # If no segments overlap with the gene, return NA values
    if (length(hits) == 0) {
      return(c("gene" = gene.dat$SYMBOL, "chr" = gene.dat$chr, 
               "start" = gene.dat$start, "end" = gene.dat$end, 
               setNames(rep(NA, 
                            length(c("nMajor", "nMinor", "nTot", "CNA", "LOH", "cnnLOH", "Amp", "HomDel"))), 
                        c("nMajor", "nMinor", "nTot", "CNA", "LOH", "cnnLOH", "Amp", "HomDel"))))
    }
    
    # only keep overlapping sgemnts
    chr.segments <- chr.segments[chr.segments$overlap == TRUE,]
    
    # Compute the overlap proportion for each matching segment
    overlaps <- pintersect(segments.query[queryHits(hits)], gene.subject[subjectHits(hits)])
    percentOverlap <- width(overlaps) / width(gene.subject[subjectHits(hits)])
    
    # Store the overlap proportion in the segment data
    chr.segments$percentOverlap <- percentOverlap
    
    # Select the segment with the highest overlap proportion as the best match
    best_segment <- chr.segments[which.max(percentOverlap), ]
    
    # Extract relevant segment information for the gene
    sample.gene.dat <- unlist(c("gene" = gene.dat$SYMBOL, "chr" = gene.dat$chr, "start" = gene.dat$start, 
                         "end" = gene.dat$end, best_segment[,c("nMajor", "nMinor", "nTot", "CNA", "LOH", "cnnLOH", "Amp", "HomDel")]))
    
    return(sample.gene.dat)
  })
  
  # Convert the results to a data frame with samples as rows
  gene.res <- do.call(rbind, lapply(gene.statuses, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  gene.res$sample <- rownames(gene.res)  # Assign sample names
  rownames(gene.res) <- NULL  # Remove row names
  
  return(gene.res)  # Return the processed gene data
}
                           
# Process all genes in parallel
res.list <- future_lapply(seq_len(nrow(genes)), process_gene)

# Combine results into sample-wise data frames
sample.res <- lapply(names(ascat.list), function(sampleID) {
  do.call(rbind, lapply(res.list, function(df) df[df$sample == sampleID, ]))
})

# Save results
save(sample.res, file = outfile.1)

print(Sys.time() - start.time)
