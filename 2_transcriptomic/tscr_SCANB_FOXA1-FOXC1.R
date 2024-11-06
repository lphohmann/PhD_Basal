# Script: Follow-up observed heterogeneity in FOXA1/FOXC1 expression 
# Author: Lennart Hohmann
# Date: 06.11.2024
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
source("./scripts/2_transcriptomic/src/tscr_functions.R")
if (!require("pacman")) install.packages("pacman")
#pacman::p_load()
#-------------------
# set/create output directories
# for plots
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.4 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGS_sampleIDs.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_FOXA1-FOXC1.pdf")
txt.file <- paste0(output.path,cohort,"_FOXA1-FOXC1.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# load palette
color.palette <- loadRData(infile.2)[c("LumA","LumB","Basal")]

# load gex
gex.df <- loadRData(infile.3)
gex.df <- gex.df[, unname(unlist(sampleIDs))]

#######################################################################
# mRNA expression boxplots
#######################################################################

gene.vec <- c("FOXA1", "FOXC1") 
#setdiff(gene.vec, rownames(gex.df))

for (gene in gene.vec) {
  txt.out <- append(txt.out, c("\n",gene,"\n",
                               "\n###########################################\n"))
  gene.gex <- gex.df[gene,] 
  # subtype data
  basal.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_Basal"]))]))
  luma.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumA"]))]))
  lumb.dat <- as.numeric(as.vector(
    gene.gex[, unname(unlist(sampleIDs["ERpHER2n_LumB"]))]))
  # summary statistics
  basal.stats <- get_stats(basal.dat)
  luma.stats <- get_stats(luma.dat)
  lumb.stats <- get_stats(lumb.dat)
  txt.out <- append(txt.out, c("Basal\n",capture.output(basal.stats), "\n",
                               "LumA\n",capture.output(luma.stats), "\n",
                               "LumB\n",capture.output(lumb.stats),
                               "\n###########################################\n"))
  # mann whitney u tests
  luma.res <- wilcox.test(basal.dat, luma.dat)
  lumb.res <- wilcox.test(basal.dat, lumb.dat)
  txt.out <- append(txt.out, c(capture.output(luma.res), "\n###########################################\n"))
  txt.out <- append(txt.out, c(capture.output(lumb.res), "\n###########################################\n"))
  # plot
  plot.par <- list(
    data = list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
    col = color.palette, 
    names = names(color.palette),
    ylab = "mRNA expression (log2)",
    main = gene)
  plot.parameters <- append(plot.parameters, list(plot.par))
}

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))
for (i in 1:length(plot.parameters)) {
  bp <- boxplot(plot.parameters[[i]]$data,
                col = plot.parameters[[i]]$col,
                names = plot.parameters[[i]]$names,
                ylab = plot.parameters[[i]]$ylab,
                main = plot.parameters[[i]]$main)
  axis(3,at=1:length(bp$n),labels=bp$n)
}

#######################################################################
# mRNA expression scatterplot
#######################################################################

plot.dat <- as.data.frame(t(gex.df[c("FOXA1","FOXC1"),]))
plot.dat$sampleID <- row.names(plot.dat)
plot.dat$PAM50 <- NA
plot.dat$PAM50 <- ifelse(plot.dat$sampleID %in% sampleIDs$ERpHER2n_Basal, "Basal",
                         ifelse(plot.dat$sampleID %in% sampleIDs$ERpHER2n_LumA, "LumA",
                                ifelse(plot.dat$sampleID %in% sampleIDs$ERpHER2n_LumB, "LumB", NA)))
# bring basal to foreground
plot(plot.dat$FOXA1, plot.dat$FOXC1, type = "n", 
     xlab = "FOXA1", ylab = "FOXC1")
points(plot.dat$FOXA1[which(plot.dat$PAM50 == "LumA")],
       plot.dat$FOXC1[which(plot.dat$PAM50 == "LumA")],
       col = color.palette["LumA"], pch = 16)
points(plot.dat$FOXA1[which(plot.dat$PAM50 == "LumB")],
       plot.dat$FOXC1[which(plot.dat$PAM50 == "LumB")],
       col = color.palette["LumB"], pch = 16)
# Plot Basal last to keep it on top
points(plot.dat$FOXA1[which(plot.dat$PAM50 == "Basal")],
       plot.dat$FOXC1[which(plot.dat$PAM50 == "Basal")],
       col = color.palette["Basal"], pch = 16, cex= 1.3)

#######################################################################
# mRNA expression by HRD status in Basal-like cases
#######################################################################

hrd.dat <- loadRData(infile.4)
hrd.dat$sampleID <- hrd.dat$Lund.tumour.id
hrd.dat <- merge(plot.dat, hrd.dat, by = "sampleID")
hrd.dat <- hrd.dat[c("sampleID","FOXA1","FOXC1","HRDetect")]
#str(hrd.dat)

# Boxplot for FOXA1 split by HRDetect status
bp.foxa1 <- boxplot(FOXA1 ~ HRDetect, data = hrd.dat,
        col = c("lightblue", "lightcoral"),  # Custom colors for each group
        xlab = "HRD Status", ylab = "FOXA1 mRNA expression")
axis(3,at=1:length(bp.foxa1$n),labels=bp.foxa1$n)
# Boxplot for FOXC1 split by HRDetect status
bp.foxc1 <- boxplot(FOXC1 ~ HRDetect, data = hrd.dat,
        col = c("lightblue", "lightcoral"),  # Custom colors for each group
        xlab = "HRD Status", ylab = "FOXC1 Expression")
axis(3,at=1:length(bp.foxc1$n),labels=bp.foxc1$n)

#######################################################################

#plot.list <- append(plot.list, list(plot))
# for(i in 1:length(plot.list)) {
#   print(i)
#   print(plot.list[[i]])
# }
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)