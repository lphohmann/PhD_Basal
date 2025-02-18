# Script: Plotting CNA boxplots in SCAN-B/BASIS
# Author: Lennart Hohmann
# Date: 15.02.2024
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
#pacman::p_load()
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
infile.5 <- "./data/Parameters/color_palette.RData"
infile.6 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.7 <- "./data/SCANB/4_CN/processed/CNA_genetest.RData"
infile.8 <- "./data/SCANB/4_CN/processed/CNA_genelevel_all.RData"
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAboxplots.pdf")
txt.file <- paste0(output.path,cohort,"_CNAboxplots.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.5)[c("LumA","LumB","Basal")]

# load Basal ids
#basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load
gl.freqs <- loadRData(infile.6)

###############################################################################
# signif gene data
###############################################################################

# prep dat
gene.test.df <- loadRData(infile.7)

#length(gene.test.df$gene)
#length(genes.AG$gene)+length(genes.AL$gene)

# subset signif genes
genes.AG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumA.Gain.padj <= 0.05, ]$gene, ]

genes.AL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumA.Loss.padj <= 0.05, ]$gene, ]

genes.BG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumB.Gain.padj <= 0.05, ]$gene, ]

genes.BL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$LumB.Loss.padj <= 0.05, ]$gene, ]

###############################################################################
# only plot autosomes
###############################################################################
#str(chr.lengths)
gl.freqs <- gl.freqs[gl.freqs$chr != 23,]
genes.AG <- genes.AG[genes.AG$chr != 23,]
genes.AL <- genes.AL[genes.AL$chr != 23,]
genes.BG <- genes.BG[genes.BG$chr != 23,]
genes.BL <- genes.BL[genes.BL$chr != 23,]
#chr.lengths <- chr.lengths[chr.lengths$Chr != 23,]

###############################################################################
# plot: all profiles with points specific to Basal
###############################################################################

# # subset signif genes
genes.all.gain <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.AG$gene,genes.BG$gene), ]
genes.all.loss <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.AL$gene,genes.BL$gene), ]
# 
# # Pick higher frequency between freq.gain.luma and freq.gain.her2e
# genes.all.gain$y <- apply(genes.all.gain[, c("freqgain.Basal", "freqgain.LumA", "freqgain.LumB")], 1, max)
# genes.all.loss$y <- apply(genes.all.loss[, c("freqloss.Basal", "freqloss.LumA", "freqloss.LumB")], 1, min)
# 
# ####### ADAPT 

################################################################################
# CN boxplot freq differences of significant genes
################################################################################

# 2 data objects with luma and lumb
genes.AG.diff <- abs(genes.AG$freqgain.Basal - genes.AG$freqgain.LumA)
genes.AL.diff <- abs(genes.AL$freqloss.Basal - genes.AL$freqloss.LumA)
luma.dat <- c(genes.AG.diff,genes.AL.diff)
genes.BG.diff <- abs(genes.BG$freqgain.Basal - genes.BG$freqgain.LumB)
genes.BL.diff <- abs(genes.BL$freqloss.Basal - genes.BL$freqloss.LumB)
lumb.dat <- c(genes.BG.diff,genes.BL.diff)

txt.out <- append(txt.out, c("\nGenes with significant difference in CNA Frequency compared to Basal-like\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("LumA Gain: n=",
                             paste0(nrow(genes.AG)," (",nrow(genes.AG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumA Loss: n=",
                             paste0(nrow(genes.AL)," (",nrow(genes.AL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumB Gain: n=",
                             paste0(nrow(genes.BG)," (",nrow(genes.BG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "LumB Loss: n=",
                             paste0(nrow(genes.BL)," (",nrow(genes.BL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Gain: n=",
                             paste0(nrow(genes.all.gain), " (",nrow(genes.all.gain)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Loss: n=",
                             paste0(nrow(genes.all.loss), " (",nrow(genes.all.loss)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Abs. freq diff to Basal, median LumA=",
                             median(luma.dat),"\n",
                             "Abs. freq diff to Basal, median LumB=",
                             median(lumb.dat),"\n",
                             "\n###########################################\n"))

#median(luma.dat)
#median(lumb.dat)

plot.par <- list(
  data = list(LumA=luma.dat,LumB=lumb.dat), 
  col = color.palette[c("LumA","LumB")], 
  names = names(color.palette[c("LumA","LumB")]),
  ylab = "CNAFreq diff to Basal",
  main = "Abs. diff in CNA Freq to Basal")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################
# CN boxplot of % genome altered
################################################################################

cna.df <- loadRData(infile.8)[[1]]
row.names(cna.df) <- cna.df$gene
sample.ids <- loadRData(infile.8)[[2]]

cna.df <- cna.df[cna.df$chr != 23,]

# sample # %altered
basal.dat <- as.vector(apply(cna.df[sample.ids$Basal],2,function(x) {
  (sum(x!=0)/length(x))*100}))
luma.dat <-  as.vector(apply(cna.df[sample.ids$LumA],2,function(x) {
  (sum(x!=0)/length(x))*100}))
lumb.dat <-  as.vector(apply(cna.df[sample.ids$LumB],2,function(x) {
  (sum(x!=0)/length(x))*100}))

##
txt.out <- append(txt.out, c("\nGenome altered (%) compared to Basal\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("Basal: mean=", round(mean(basal.dat),2),"\n",
                             "LumA: mean=", round(mean(luma.dat),2),"\n",
                             "LumB: mean=", round(mean(lumb.dat),2),"\n",
                             "\n###########################################\n"))

# mann whitney u tests
luma.res <- wilcox.test(basal.dat, luma.dat) 
lumb.res <- wilcox.test(basal.dat, lumb.dat)

txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))

plot.par <- list(
  data = list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
  col = color.palette, 
  names = names(color.palette),
  ylab = "Genome altered (%)",
  main = "Genome altered")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################

# num singif gener per chromosme
a.dat <- rbind(genes.AG,genes.AL)
a.dat$chr <- as.numeric(a.dat$chr)
a.tbl <- table(a.dat$chr)

barplot(height=a.tbl, # num vec with gener per chromosome 
        names=names(a.tbl), 
        main="LumA: Signif genes by chromosme",
        ylab="Signif altered genes")

b.dat <- rbind(genes.BG,genes.BL)
b.dat$chr <- as.numeric(b.dat$chr)
b.tbl <- table(b.dat$chr)

barplot(height=b.tbl, # num vec with gener per chromosome 
        names=names(b.tbl), 
        main="LumB: Signif genes by chromosme",
        ylab="Signif altered genes")


################################################################################
################################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                main = plot.parameters[[i]]$main,
                ylim = c(0,100))
  axis(3,at=1:length(bp$n),labels=bp$n)
}
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)
