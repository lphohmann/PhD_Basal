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
infile.2 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.3 <- "./data/BASIS/1_clinical/raw/Summarized_Annotations_BASIS.RData"
infile.4 <- "./data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv"
infile.5 <- "./data/Parameters/color_palette.RData"


#infile.3 <- "data/BASIS/4_CN/raw/GRCh38_EBV.chrom.sizes.tsv" # move ot basal project
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAprofile.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.5)[c("LumA","LumB","Basal")]

# load Basal ids
basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load
gl.freqs <- loadRData(infile.2)

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
gl.freqs[c("chr","start","end")] <- lapply(gl.freqs[c("chr","start","end")], as.numeric)
gl.freqs$centerPos <- (gl.freqs$start + gl.freqs$end) / 2

# convert to genome position
gl.freqs$genome <- chr.lengths$genome[match(gl.freqs$chr,chr.lengths$Chr)]

# Add a column with the genome position of each probe
gl.freqs$Genome_pos <- gl.freqs$centerPos + gl.freqs$genome

###############################################################################
# plot: all profiles without points
###############################################################################

pdf(file = plot.file, height = 21.0, width = 72.0)

plot <- ggplot() +  
  ggtitle("Genome-wide frequency of gain/loss CN alterations") +
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$LumA_Gain, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$LumA_Loss, 
    color = "LumA"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$LumB_Gain, 
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$LumB_Loss,  
    color = "LumB"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$Basal_Loss,  
    color = "Basal"),size=4) + 
  geom_line(aes(
    x = gl.freqs$Genome_pos, 
    y = gl.freqs$Basal_Gain,  
    color = "Basal"),size=4) + 
  scale_colour_manual(name="Subtype", values = color.palette) + 
  geom_vline(xintercept = chr.lengths$genome[-length(chr.lengths$genome)],
             linetype="dashed",size=1) + 
  scale_x_continuous(name="Genome position (chromosome)",
                     breaks=chr.lengths$genome, 
                     labels=as.character(1:23),
                     limits = c(0,max(chr.lengths$genome)), #+50000000
                     expand = c(0, 0)) +
  scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
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
# ###############################################################################
# # signif gene data
# ###############################################################################
# 
# # prep dat
# gene.test.df <- loadRData("data/COMBINED/4_CN/processed/CNA_genelevel.RData")
# genes.AG <- gene.freqs %>% 
#   filter(gene %in% gene.test.df[gene.test.df$LumA.Gain.padj<=0.05,]$gene) %>% 
#   mutate(y = ifelse(
#     freq.gain.luma>freq.gain.her2e,freq.gain.luma,freq.gain.her2e)) #pick what is higher luma or her2e
# genes.BG <- gene.freqs %>% 
#   filter(gene %in% gene.test.df[gene.test.df$LumB.Gain.padj<=0.05,]$gene) %>% 
#   mutate(y = ifelse(
#     freq.gain.lumb>freq.gain.her2e,freq.gain.lumb,freq.gain.her2e))
# 
# genes.AL <- gene.freqs %>% 
#   filter(gene %in% gene.test.df[gene.test.df$LumA.Loss.padj<=0.05,]$gene) %>% 
#   mutate(y = ifelse(
#     freq.loss.luma<freq.loss.her2e,freq.loss.luma,freq.loss.her2e))
# genes.BL <- gene.freqs %>% 
#   filter(gene %in% gene.test.df[gene.test.df$LumB.Loss.padj<=0.05,]$gene) %>% 
#   mutate(y = ifelse(
#     freq.loss.lumb<freq.loss.her2e,freq.loss.lumb,freq.loss.her2e))
# 
# 
# ###############################################################################
# # plot: all profiles with points
# ###############################################################################
# 
# #pdf(file = paste(output.path,cohort,"_GLsubtypeprofiles.pdf", sep =""), height = 21.0, width = 72.0)
# 
# plot <- ggplot() +  
#   ggtitle("Genome-wide frequency of gain/loss CN alterations") +
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.gain.luma)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.gain.luma),]$freq.gain.luma, 
#     color = "LUMA"),size=4) + 
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.loss.luma)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.loss.luma),]$freq.loss.luma, 
#     color = "LUMA"),size=4) + 
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.gain.lumb)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.gain.lumb),]$freq.gain.lumb, 
#     color = "LUMB"),size=4) + 
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.loss.lumb)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.loss.lumb),]$freq.loss.lumb, 
#     color = "LUMB"),size=4) + 
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.gain.her2e)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.gain.her2e),]$freq.gain.her2e, 
#     color = "HER2E"),size=4) + 
#   geom_line(aes(
#     x = cn.data[which(!is.na(cn.data$freq.loss.her2e)),]$Genome_pos, 
#     y = cn.data[!is.na(cn.data$freq.loss.her2e),]$freq.loss.her2e, 
#     color = "HER2E"),size=4) + 
#   scale_colour_manual(name="Subtype", values = c("HER2E"="#d334eb", "LUMA"="#2176d5", "LUMB"="#34c6eb")) + 
#   geom_point(aes(x = genes.AG$Genome_pos, y = genes.AG$y), size=12) +
#   geom_point(aes(x = genes.AL$Genome_pos, y = genes.AL$y), size=12) +
#   geom_point(aes(x = genes.BG$Genome_pos, y = genes.BG$y), size=12, colour="red") +
#   geom_point(aes(x = genes.BL$Genome_pos, y = genes.BL$y), size=12, colour="red") +
#   geom_vline(xintercept = chr.lengths$genome[-length(chr.lengths$genome)],
#              linetype="dashed",size=1) +
#   scale_x_continuous(name="Genome position (chromosome)",
#                      breaks=chr.lengths$chrbreaks, 
#                      labels=as.character(1:22),
#                      limits = c(0,max(chr.lengths$genome)), #+50000000
#                      expand = c(0, 0)) +
#   scale_y_continuous(name="Alteration frequency (%)", # \n Loss          Gain", # \n Loss & G
#                      breaks=c(seq(-100,100,25)),
#                      labels=c(100,75,50,25,0,25,50,75,100),
#                      expand = c(0, 0),
#                      limits = c(-100,100)) +
#   theme_bw() +
#   theme(text=element_text(size=30),
#         legend.title = element_blank(),
#         axis.title.y = element_text(vjust = 0.5),
#         legend.position = c(0.97, 0.95),
#         panel.border = element_blank(), 
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line = element_line(colour = "black",linewidth=2),
#         axis.ticks = element_line(colour = "black", linewidth = 2),
#         axis.ticks.length=unit(0.5, "cm")) + #legend.position = "none") +
#   annotate(x=min(cn.data$Genome_pos)+30000000,
#            y=c(-50,50), label=c("Loss","Gain"), 
#            geom="text", angle=90, hjust=0.5, 
#            size=9, colour=c("black","black")) 
# print(plot)
# 
