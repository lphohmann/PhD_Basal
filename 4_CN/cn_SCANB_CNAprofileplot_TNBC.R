# Script: TNBC Plotting genome wide CNA Frequencies in SCAN-B 
# Author: Lennart Hohmann
# Date: 01.04.2024
#TODO: check chr23, exclude?
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
#infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
#infile.2 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
#infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.4 <- "./data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
infile.5 <- "./data/Parameters/TNBC_color_palette.RData"
infile.6 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.8 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_TNBC.RData"
infile.9 <- "./data/SCANB/4_CN/processed/CNA_genetest_TNBC.RData"
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAprofile_TNBC.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
#plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.5)[c("TNBC_NonBasal","TNBC_Basal","ERpHER2n_Basal")]

# load Basal ids
#basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load
gl.freqs.basal <- loadRData(infile.6)
gl.freqs.basal <- gl.freqs.basal[c("gene","chr","start","end","freqloss.Basal","freqgain.Basal")]
gl.freqs.tnbc <- loadRData(infile.8)

head(gl.freqs.basal)
head(gl.freqs.tnbc)

# only include common genes in ASCAT
common.genes <- intersect(gl.freqs.basal$gene,gl.freqs.tnbc$gene)
gl.freqs.basal <- gl.freqs.basal[gl.freqs.basal$gene %in% common.genes, ]
gl.freqs.tnbc <- gl.freqs.tnbc[gl.freqs.tnbc$gene %in% common.genes, ]
identical(gl.freqs.basal$gene,gl.freqs.tnbc$gene)
# order
gl.freqs.basal <- gl.freqs.basal[order(gl.freqs.basal$gene), ]
gl.freqs.tnbc <- gl.freqs.tnbc[order(gl.freqs.tnbc$gene), ]
identical(gl.freqs.tnbc$gene,gl.freqs.basal$gene)

# 
gl.freqs.all <- cbind(gl.freqs.basal,gl.freqs.tnbc[5:ncol(gl.freqs.tnbc)])

#basis.anno <- loadRData(infile.3)
#basis.anno <- basis.anno[basis.anno$ClinicalGroup == "ERposHER2neg" & basis.anno$PAM50_AIMS %in% c("LumA","LumB"),c("sample_name","PAM50_AIMS")]

#######################################################################
# add genome center position for genes
#######################################################################

# Read data from file
chr.lengths <- as.data.frame(read.table(file = infile.4, sep = '\t', header = FALSE))[1:23, ]

# Rename columns
names(chr.lengths) <- c("Chr", "length")

# Create 'genome' column
chr.lengths$genome <- cumsum(as.numeric(chr.lengths$length)) # lag by 1 position (cause I have to add the length of the previous chr to the probe positions (0 for chr1 probes))

# Shift 'genome' column by one row
chr.lengths$genome <- c(0, chr.lengths$genome[-nrow(chr.lengths)])

chr.lengths[chr.lengths$Chr=="chrX",]$Chr <- "chr23"
chr.lengths$Chr <- as.numeric(gsub('^.{3}','',chr.lengths$Chr))

# add center pos
gl.freqs.all[c("chr","start","end")] <- lapply(gl.freqs.all[c("chr","start","end")], as.numeric)
gl.freqs.all$centerPos <- (gl.freqs.all$start + gl.freqs.all$end) / 2

# convert to genome position
gl.freqs.all$genome <- chr.lengths$genome[match(gl.freqs.all$chr,chr.lengths$Chr)]

# Add a column with the genome position of each probe
gl.freqs.all$Genome_pos <- gl.freqs.all$centerPos + gl.freqs.all$genome

###############################################################################
# signif gene data
###############################################################################

# prep dat
gene.test.df <- loadRData(infile.9)

# subset signif genes
genes.TB.gain <- gl.freqs.all[gl.freqs.all$gene %in% 
                       gene.test.df[gene.test.df$TNBC.Basal.Gain.padj <= 0.05, ]$gene, ]

genes.TB.loss <- gl.freqs.all[gl.freqs.all$gene %in% 
                       gene.test.df[gene.test.df$TNBC.Basal.Loss.padj <= 0.05, ]$gene, ]

genes.TNB.gain <- gl.freqs.all[gl.freqs.all$gene %in% 
                       gene.test.df[gene.test.df$TNBC.NonBasal.Gain.padj <= 0.05, ]$gene, ]

genes.TNB.loss <- gl.freqs.all[gl.freqs.all$gene %in% 
                       gene.test.df[gene.test.df$TNBC.NonBasal.Loss.padj <= 0.05, ]$gene, ]

# Pick higher frequency between freq.gain.luma and freq.gain.her2e
genes.TB.gain$y <- with(genes.TB.gain,ifelse(freqgain.tnbc.Basal > freqgain.Basal,
                                             freqgain.tnbc.Basal, freqgain.Basal))
                   
genes.TB.loss$y <- with(genes.TB.loss,ifelse(freqloss.tnbc.Basal < freqloss.Basal,
                                             freqloss.tnbc.Basal, freqloss.Basal))

genes.TNB.gain$y <- with(genes.TNB.gain,ifelse(freqgain.tnbc.NonBasal > freqgain.Basal,
                                               freqgain.tnbc.NonBasal, freqgain.Basal))

genes.TNB.loss$y <- with(genes.TNB.loss,ifelse(freqloss.tnbc.NonBasal < freqloss.Basal,
                                               freqloss.tnbc.NonBasal, freqloss.Basal))

gl.freqs <- gl.freqs.all
###############################################################################
# only plot autosomes
###############################################################################
#str(chr.lengths)
gl.freqs <- gl.freqs[gl.freqs$chr != 23,]
genes.TB.gain <- genes.TB.gain[genes.TB.gain$chr != 23,]
genes.TB.loss <- genes.TB.loss[genes.TB.loss$chr != 23,]
genes.TNB.gain <- genes.TNB.gain[genes.TNB.gain$chr != 23,]
genes.TNB.loss <- genes.TNB.loss[genes.TNB.loss$chr != 23,]
chr.lengths <- chr.lengths[chr.lengths$Chr != 23,]

###############################################################################
# plot: all profiles with points 
###############################################################################

pdf(file = plot.file, height = 21.0, width = 72.0)

plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.tnbc.Basal, 
    color = "TNBC_Basal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.tnbc.Basal, 
    color = "TNBC_Basal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.tnbc.NonBasal, 
    color = "TNBC_NonBasal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.tnbc.NonBasal,  
    color = "TNBC_NonBasal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqgain.Basal,  
    color = "ERpHER2n_Basal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$freqloss.Basal,  
    color = "ERpHER2n_Basal"),size=4) + 
  scale_colour_manual(name="Subtype", values = color.palette) + 
  geom_point(aes(x = genes.TB.gain$Genome_pos, y = genes.TB.gain$y), size=20, colour="black") +
  geom_point(aes(x = genes.TB.loss$Genome_pos, y = genes.TB.loss$y), size=20, colour="black") +
  geom_point(aes(x = genes.TNB.gain$Genome_pos, y = genes.TNB.gain$y), size=12, colour="blue") +
  geom_point(aes(x = genes.TNB.loss$Genome_pos, y = genes.TNB.loss$y), size=12, colour="blue") +
  geom_vline(xintercept = chr.lengths$genome,
             linetype="dashed",size=1) + 
  scale_x_continuous(name="Genome position (chromosome)",
                     breaks=chr.lengths$genome, 
                     labels=as.character(1:22), #as.character(1:23),
                     limits = c(0,max(chr.lengths$genome)+50000000), #+50000000
                     expand = c(0, 0)) +
  scale_y_continuous(name="Alteration frequency (%)",
                     breaks=c(seq(-100,100,25)),
                     labels=c(100,75,50,25,0,25,50,75,100),
                     expand = c(0, 0),
                     limits = c(-100,100)) +
  theme_bw() +
  theme(text=element_text(size=30),
        legend.title = element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        legend.position = c(0.97, 0.95),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black",linewidth=2),
        axis.ticks = element_line(colour = "black", linewidth = 2),
        axis.ticks.length=unit(0.5, "cm")) + #legend.position = "none") +
  annotate(x=min(gl.freqs$Genome_pos)+30000000,
           y=c(-50,50), label=c("Loss","Gain"), 
           geom="text", angle=90, hjust=0.5, 
           size=9, colour=c("black","black")) 
print(plot)
dev.off()

###############################################################################
# plot: all profiles with points specific to Basal
###############################################################################

# subset signif genes
genes.all.gain <- gl.freqs[gl.freqs$gene %in% 
                             intersect(genes.TB.gain$gene,genes.TNB.gain$gene), ]
genes.all.loss <- gl.freqs[gl.freqs$gene %in% 
                             intersect(genes.TB.loss$gene,genes.TNB.loss$gene), ]