# Script: Map probe coordinates to genes in BASIS for later comparions with SCANB-Basal
# Author: Lennart Hohmann
# Date: 12.02.2024
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
pacman::p_load(rtracklayer,
               GenomicRanges,
               rtracklayer,
               Repitools,
               GenomicFeatures,
               TxDb.Hsapiens.UCSC.hg19.knownGene,
               org.Hs.eg.db,
               plyranges)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/4_CN/"
dir.create(output.path)
# for data
data.path <- "./data/BASIS/4_CN/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/BASIS/4_CN/raw/LumA_CollectedFrequencyData.RData"
infile.2 <- "./data/BASIS/4_CN/raw/LumB_CollectedFrequencyData.RData"
#infile.3 <- 

# output paths
outfile.1 <- paste0(data.path,"ASCAT_genelevel.RData.RData")
#plot.file <- paste0(output.path,cohort,"_i.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

cn.luma <- loadRData("./data/BASIS/4_CN/raw/LumA_CollectedFrequencyData.RData")
cn.luma <- do.call("cbind", list(cn.luma$fData,"LumA_Gain"=cn.luma$CN_Gain,"LumA_Loss"=cn.luma$CN_Loss))

cn.lumb <- loadRData("./data/BASIS/4_CN/raw/LumB_CollectedFrequencyData.RData")
cn.lumb <- do.call("cbind", list(cn.lumb$fData,"LumB_Gain"=cn.lumb$CN_Gain,"LumB_Loss"=cn.lumb$CN_Loss))

cn.basis <- merge(cn.luma,cn.lumb,by=c("reporterId","chromosome","centerPosition"))
#View(cn.basis)

#######################################################################
# map probes to genes
#######################################################################

probe.positions <- cn.basis[1:3]
probe.positions$chromosome <- paste0("chr",probe.positions$chromosome)
probe.positions$chromosome[probe.positions$chromosome=="chr23"] <- "chrX"

# create granges objects
probes <- GRanges(seqnames = probe.positions$chromosome,
                  ranges = IRanges(probe.positions$centerPosition),
                  reporterId = probe.positions$reporterId)

# get hg19 gnes
genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) # meta = entrez ID

# convert to hgnc symbols
ENTREZID2SYMBOL <- biomaRt::select(org.Hs.eg.db, mcols(genes)$gene_id, c("ENTREZID", "SYMBOL"))
stopifnot(identical(ENTREZID2SYMBOL$ENTREZID, mcols(genes)$gene_id))
mcols(genes)$SYMBOL <- ENTREZID2SYMBOL$SYMBOL
#View(genes)

genes <- filter(genes,seqnames %in% unique(probe.positions$chromosome))
gene.anno <- as.data.frame(genes)
probe.anno <- as.data.frame(probes)

# which 
overlap.res <- findOverlaps(genes,probes)
# queryHits(): indexes of the gene coordinates that overlap the corresponding 
# subjectHits(): indexes of the probes
# line up the query column identifier (gene) that overlaps each probe
f1 <- factor(subjectHits(overlap.res), levels=seq_len(subjectLength(overlap.res)))
# use of factor() with exactly as many levels as there are subjects ensures that the splitAsList() command returns a 1:1 mapping between the subjects (probes) and the genes in the corresponding CharacterList
overlap.list <- splitAsList(mcols(genes)[["SYMBOL"]][queryHits(overlap.res)], f1) # split the column of gene symbols into lists corresponding to the regions of overlap
mcols(probes) <- overlap.list

key.df <- merge(as.data.frame(probes),probe.anno,
                by=c("seqnames","start","end","width","strand")) 

names(key.df)[names(key.df) == "X"] <- "Gene_symbol"
names(key.df)[names(key.df) == "start"] <- "Position"
names(key.df)[names(key.df) == "seqnames"] <- "Chr"

key.df <- key.df[c("reporterId","Gene_symbol")]

# add to matrix
probe.positions <- merge(probe.positions,key.df,by="reporterId") 

head(probe.positions)

# replace character(0) with NA for later filtering of probes without annotation
probe.positions$Gene_symbol <- lapply(probe.positions$Gene_symbol, 
                                      function(x) {
                                        if (identical(x, character(0))) {
                                          return(NA)
                                        } else {return(x)}
                                      })

################################################################################
# prep data for saving
################################################################################
# key file for gene to probe mapping - make one for BASIS or just use table14? 
#map.key <- loadRData("data//4_CN/processed/CN_mapped_probes.RData")


# add "chr"
#cnfile$Chr <- sub("^", "chr", cnfile$Chr) 

# merge with gene positions
cn.basis.mapped <- merge(cn.basis,probe.positions[c("reporterId", "Gene_symbol")],by="reporterId") 

# filter NA rows
cn.basis.mapped <- cn.basis.mapped[which(!is.na(cn.basis.mapped$Gene_symbol)),]

cn.basis.mapped <- cn.basis.mapped[c("Gene_symbol","LumA_Gain",
                                     "LumA_Loss","LumB_Gain","LumB_Loss")]

# save
save(cnfile, 
     file = )


