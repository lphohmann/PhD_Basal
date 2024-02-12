# Script: Plotting genome wide CNA Frequencies in SCAN-B 
# Author: Lennart Hohmann
# Date: 11.02.2024
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
pacman::p_load(
               IRanges,
               GenomicFeatures,
               TxDb.Hsapiens.UCSC.hg38.knownGene,
               TxDb.Hsapiens.UCSC.hg19.knownGene,
               org.Hs.eg.db,
               Repitools)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/4_CN/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/4_CN/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/3_WGS/processed/ASCAT_genelevel.RData"

#infile.3 <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv" # move ot basal project
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

# load Basal ids
basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load ASCAT gene data
ascat.dat <- loadRData(infile.2)
names(ascat.dat) <- gsub("\\..*", "", names(ascat.dat))
ascat.dat <- ascat.dat[names(ascat.dat) %in% basal.ids]



#######################################################################
# add center position for later plotting; or do in different script?
#######################################################################
# why do i need that

# chr lengths
# get chr lengths to get genome positions of probes (excluding chr X)
chr.lengths <- as.data.frame(read.table(file = infile.3, sep = '\t', header = FALSE))[1:22,] %>% 
  dplyr::rename(chr=V1,length=V2) %>% 
  mutate(genome = cumsum(as.numeric(length))) %>% 
  mutate(genome = lag(genome,default = 0)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))
chr.lengths$chr <- as.numeric(gsub('^.{3}','',chr.lengths$chr))


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


final.names <- c("gene","chr","start","end","CNA","LOH","cnnLOH","Amp","HomDel","sample")

final.df <- data.frame(matrix(nrow = 0, ncol = length(final.names), 
                              dimnames = list(NULL, final.names)))
View(rbind(final.df,gene.res))

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
