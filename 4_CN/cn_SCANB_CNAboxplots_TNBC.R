# Script: Plotting CNA boxplots in SCAN-B/BASIS tnbc
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
infile.5 <- "./data/Parameters/TNBC_color_palette.RData"
infile.6 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_all.RData"
infile.7 <- "./data/SCANB/4_CN/processed/CNA_genetest_TNBC.RData"
infile.8 <- "./data/SCANB/4_CN/processed/CNA_genelevel_all.RData"
infile.9 <- "./data/SCANB/4_CN/processed/CNA_GLFreqs_TNBC.RData"
infile.10 <- "./data/SCANB/4_CN/processed/CNA_genelevel_TNBC.RData"
# output paths
#outfile.1 <- ""
plot.file <- paste0(output.path,cohort,"_CNAboxplots_TNBC.pdf")
txt.file <- paste0(output.path,cohort,"_CNAboxplots_TNBC.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.5)#[c("LumA","LumB","Basal")]

# load Basal ids
#basal.ids <- unname(unlist(loadRData(infile.1)["ERpHER2n_Basal"]))

# load
gl.freqs.basal <- loadRData(infile.6)
gl.freqs.basal <- gl.freqs.basal[c("gene","chr","start","end","freqloss.Basal","freqgain.Basal")]


gl.freqs.tnbc <- loadRData(infile.9)

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
gl.freqs <- cbind(gl.freqs.basal,gl.freqs.tnbc[5:ncol(gl.freqs.tnbc)])

###############################################################################
# signif gene data
###############################################################################

# prep dat
gene.test.df <- loadRData(infile.7)
# subset signif genes
genes.bG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$TNBC.Basal.Gain.padj <= 0.05, ]$gene, ]

genes.bL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$TNBC.Basal.Loss.padj <= 0.05, ]$gene, ]

genes.nbG <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$TNBC.NonBasal.Gain.padj <= 0.05, ]$gene, ]

genes.nbL <- gl.freqs[gl.freqs$gene %in% 
                       gene.test.df[gene.test.df$TNBC.NonBasal.Loss.padj <= 0.05, ]$gene, ]

###############################################################################
# plot: all profiles with points specific to Basal
###############################################################################

# # subset signif genes
genes.all.gain <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.bG$gene,genes.nbG$gene), ]
genes.all.loss <- gl.freqs[gl.freqs$gene %in%
                             intersect(genes.bL$gene,genes.nbL$gene), ]
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
genes.bG.diff <- abs(genes.bG$freqgain.Basal - genes.bG$freqgain.tnbc.Basal)
genes.bL.diff <- abs(genes.bL$freqloss.Basal - genes.bL$freqloss.tnbc.Basal)
b.dat <- c(genes.bG.diff,genes.bL.diff)
genes.nbG.diff <- abs(genes.nbG$freqgain.Basal - genes.nbG$freqgain.tnbc.NonBasal)
genes.nbL.diff <- abs(genes.nbL$freqloss.Basal - genes.nbL$freqloss.tnbc.NonBasal)
nb.dat <- c(genes.nbG.diff,genes.nbL.diff)

txt.out <- append(txt.out, c("\nGenes with significant difference in CNA Frequency compared to Basal-like\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("TN_B Gain: n=",
                             paste0(nrow(genes.bG)," (",nrow(genes.bG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "TN_B Loss: n=",
                             paste0(nrow(genes.bL)," (",nrow(genes.bL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "TN_NB Gain: n=",
                             paste0(nrow(genes.nbG)," (",nrow(genes.nbG)/nrow(gl.freqs)*100,"%)"),"\n",
                             "TN_NB Loss: n=",
                             paste0(nrow(genes.nbL)," (",nrow(genes.nbL)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Gain: n=",
                             paste0(nrow(genes.all.gain), " (",nrow(genes.all.gain)/nrow(gl.freqs)*100,"%)"),"\n",
                             "Both intersected Loss: n=",
                             paste0(nrow(genes.all.loss), " (",nrow(genes.all.loss)/nrow(gl.freqs)*100,"%)"),"\n",
                             "\n###########################################\n"))


plot.par <- list(
  data = list(TNBC_Basal=b.dat,TNBC_NonBasal=nb.dat), 
  col = color.palette[c("TNBC_Basal","TNBC_NonBasal")], 
  names = c("TNBC_Basal","TNBC_NonBasal"),
  ylab = "CNAFreq diff to Basal",
  main = "Abs. diff in CNA Freq to Basal")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################
# CN boxplot of % genome altered
################################################################################

tnbc.cna.df <- loadRData(infile.10)[["ascat.df.tnbc"]]
tnbc.sample.ids <- loadRData(infile.10)[["subtype.samples"]]
cna.df <- loadRData(infile.8)[[1]]
row.names(cna.df) <- cna.df$gene
sample.ids <- loadRData(infile.8)[[2]]

# sample # %altered
basal.dat <- as.vector(apply(cna.df[sample.ids$Basal],2,function(x) {
  (sum(x!=0)/length(x))*100}))
b.dat <-  as.vector(apply(tnbc.cna.df[tnbc.sample.ids$TNBC.Basal],2,function(x) {
  (sum(x!=0)/length(x))*100}))
nb.dat <-  as.vector(apply(tnbc.cna.df[tnbc.sample.ids$TNBC.NonBasal],2,function(x) {
  (sum(x!=0)/length(x))*100}))

##
txt.out <- append(txt.out, c("\nGenome altered (%) compared to Basal\n",
                             "\n###########################################\n"))

txt.out <- append(txt.out, c("Basal: mean=", round(mean(basal.dat),2),"\n",
                             "TN_B: mean=", round(mean(b.dat),2),"\n",
                             "TN_NB: mean=", round(mean(nb.dat),2),"\n",
                             "\n###########################################\n"))

# mann whitney u tests
b.res <- wilcox.test(basal.dat, b.dat) 
nb.res <- wilcox.test(basal.dat, nb.dat)

txt.out <- append(txt.out, c(capture.output(b.res), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(nb.res), "\n###########################################\n"))

plot.par <- list(
  data = list(Basal=basal.dat,TNBC_Basal=b.dat,TNBC_NonBasal=nb.dat), 
  col = color.palette[c("ERpHER2n_Basal","TNBC_Basal","TNBC_NonBasal")], 
  names = c("ERpHER2n_Basal","TNBC_Basal","TNBC_NonBasal"),
  ylab = "Genome altered (%)",
  main = "Genome altered")

plot.parameters <- append(plot.parameters, list(plot.par))

################################################################################

# # num singif gener per chromosme
# a.dat <- rbind(genes.AG,genes.AL)
# a.dat$chr <- as.numeric(a.dat$chr)
# a.tbl <- table(a.dat$chr)
# 
# barplot(height=a.tbl, # num vec with gener per chromosome 
#         names=names(a.tbl), 
#         main="LumA: Signif genes by chromosme",
#         ylab="Signif altered genes")
# 
# b.dat <- rbind(genes.BG,genes.BL)
# b.dat$chr <- as.numeric(b.dat$chr)
# b.tbl <- table(b.dat$chr)
# 
# barplot(height=b.tbl, # num vec with gener per chromosome 
#         names=names(b.tbl), 
#         main="LumB: Signif genes by chromosme",
#         ylab="Signif altered genes")


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
