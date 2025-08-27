# Script: Danenberg immune part in metabric continued
# Author: Lennart Hohmann
# Date: 01.06.2024
#-------------------
# empty environment 
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "METABRIC"
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
data.path <- "./data/METABRIC/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.0 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.1 <- "./data/METABRIC/Danenberg/MBTMEIMCPublic/IMCClinical.csv"
infile.2 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTcounts.RData"
infile.3 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTprop.RData"
infile.4 <- "./data/METABRIC/2_transcriptomic/processed/SingleCells_PTmedians.RData"
infile.5 <- "./data/METABRIC/2_transcriptomic/raw/data_mRNA_median_all_sample_Zscores.txt"

# output paths
plot.file <- paste0(output.path,cohort,"_singlecellimmune.pdf")
txt.file <- paste0(output.path,cohort,"_singlecellimmune.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

anno <- read.table(infile.1,
                   sep = ",", header=TRUE)
lum.anno <- anno[which(anno$ERStatus=="pos"& anno$ERBB2_pos==FALSE),]
table(lum.anno$PAM50)
lum.anno$metabric_id[which(lum.anno$PAM50 == "Basal")]
counts.dat <- loadRData(infile.2)
prop.dat <- loadRData(infile.3)
median.dat <- loadRData(infile.4)
gex.dat <- read.table(infile.5,sep = "\t", header = TRUE)
names(gex.dat)[3:ncol(gex.dat)] <- gsub("\\.", "-", 
                                        colnames(gex.dat)[3:ncol(gex.dat)])

#######################################################################
# check mRNA expression vs pt proportions
#######################################################################

pt.plot <- function(pt.dat,gex.dat,pt,gene, plot.line=NULL) {
  pt.dat <- t(pt.dat[pt.dat$Phenotype == pt,
                       2:ncol(pt.dat)])
  pt.dat <- merge(t(gex.dat[gex.dat$Hugo_Symbol == gene,]), pt.dat, by = "row.names")
  names(pt.dat) <- c("sampleID","mRNA_gex", "PT_proportion")
  plot(pt.dat$mRNA_gex, pt.dat$PT_proportion, 
       main = pt, 
       xlab = paste0(gene," mRNA expr."), ylab = paste0(pt,"prop."), 
       pch = 19)
  if(!is.null(plot.line)) {
    abline(lm(PT_proportion ~ mRNA_gex, data = pt.dat), col = "red", lwd = 2)
    }
}
prop.dat$Phenotype

pdf(file = plot.file, onefile = TRUE, width = 7, height = 7) 
par(mfrow = c(3, 3))

pt.plot(prop.dat,gex.dat,pt="CD8^{+} T cells",gene="CD8A")
pt.plot(prop.dat,gex.dat,pt="CD4^{+} T cells",gene="CD4")
pt.plot(prop.dat,gex.dat,pt="CD15^{+}",gene="FUT4")
pt.plot(prop.dat,gex.dat,pt="HER2^{+}",gene="ERBB2")
pt.plot(prop.dat,gex.dat,pt="CD57^{+}",gene="B3GAT1")
pt.plot(prop.dat,gex.dat,pt="B cells",gene="MS4A1") #cd20
pt.plot(prop.dat,gex.dat,pt="Ki67^{+}",gene="MKI67") 
pt.plot(prop.dat,gex.dat,pt="CD38^{+} lymphocytes",gene="CD38") 

#######################################################################
# pam50 phenotype proportion boxplots
#######################################################################

color.palette <- loadRData("./data/Parameters/color_palette.RData")[c("LumA","LumB","Basal","Her2")]

for(pt in unique(prop.dat$Phenotype)) {
  pt.dat <- as.vector(prop.dat[which(prop.dat$Phenotype == pt),2:ncol(prop.dat)])
  luma.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal A")]])
  lumb.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal B")]])
  #her2e.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "HER2")]])
  basal.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Basal")]])
  
  bp <- boxplot(list(LumA=luma.dat,LumB=lumb.dat,
               Basal=basal.dat),#,Her2=her2e.dat),
          col = color.palette,
          names = names(color.palette)[1:3],
          ylab = paste0(pt," proportion"),
          main = pt)
  axis(3,at=1:length(bp$n),labels=bp$n)
}

#######################################################################
# pam50 phenotype count boxplots
#######################################################################

for(pt in unique(counts.dat$Phenotype)) {
  pt.dat <- as.vector(counts.dat[which(counts.dat$Phenotype == pt),2:ncol(counts.dat)])
  luma.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal A")]])
  lumb.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal B")]])
  #her2e.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "HER2")]])
  basal.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Basal")]])
  
  bp <- boxplot(list(LumA=luma.dat,LumB=lumb.dat,
                     Basal=basal.dat),#,Her2=her2e.dat),
                col = color.palette,
                names = names(color.palette)[1:3],
                ylab = paste0(pt," count"),
                main = pt)
  axis(3,at=1:length(bp$n),labels=bp$n)
}

#######################################################################
# combine counts
#######################################################################

# combine cd4 dc8 mac bcells
sum.rows <- c("CD4^{+} T cells",
              "CD8^{+} T cells",
              "Macrophages",
              "B cells",
              "CD4^{+} T cells & APCs",
              "Macrophages & granulocytes")
add.counts <- colSums(counts.dat[which(counts.dat$Phenotype %in% sum.rows),
                   2:ncol(counts.dat)],na.rm = TRUE)
add.prop <- colSums(prop.dat[which(prop.dat$Phenotype %in% sum.rows),
                                 2:ncol(prop.dat)],na.rm = TRUE)
combined.dat <- rbind(c("Comb.prop CD4,CD8,Mac,Bcell",add.prop), 
                      c("Comb.count CD4,CD8,Mac,Bcell",add.counts))
colnames(combined.dat)[1] <- "Phenotype"
combined.dat <- as.data.frame(combined.dat)

combined.dat[,-1] <- lapply(combined.dat[,-1], as.numeric)

pt = unique(combined.dat$Phenotype)[1]
for(pt in unique(combined.dat$Phenotype)) {
  pt.dat <- as.vector(combined.dat[which(combined.dat$Phenotype == pt),2:ncol(combined.dat)])
  luma.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal A")]])
  lumb.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Luminal B")]])
  #her2e.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "HER2")]])
  basal.dat <- unlist(pt.dat[lum.anno$metabric_id[which(lum.anno$PAM50 == "Basal")]])
  
  
  # Run Kruskal–Wallis test
  test.df <- data.frame(
    value = c(luma.dat, lumb.dat, basal.dat),
    group = factor(c(rep("LumA", length(luma.dat)),
                     rep("LumB", length(lumb.dat)),
                     rep("Basal", length(basal.dat))))
  )
  kw <- kruskal.test(value ~ group, data = test.df)
  print(paste(pt, "KW p =", signif(kw$p.value, 3)))
  
  bp <- boxplot(list(LumA=luma.dat,LumB=lumb.dat,
                     Basal=basal.dat),#,Her2=her2e.dat),
                col = color.palette,
                names = names(color.palette)[1:3],
                ylab = paste0(pt),
                main = pt)
  axis(3,at=1:length(bp$n),labels=bp$n)
}

#######################################################################
# phenotype marker median boxplots
# 1 maker in all pts
#######################################################################

par(mfrow = c(2, 1))

# final product: 1 df per marker with rows being samples and 1 column per phenotype

# extract marker col from all pts
markers <- colnames(median.dat[[1]])[3:ncol(median.dat[[1]])]
for (marker in markers) {
  #marker = markers[1]
  # create empty marker df
  # res.df <- setNames(
  #   data.frame(metabric_id = character(), marker = numeric()), 
  #   c("metabric_id", marker))
  res.ls <- list()
  for (pt in names(median.dat)) {
    pt.dat <- as.data.frame(median.dat[[pt]])
    marker.pt.df <- data.frame(metabric_id = pt.dat$metabric_id,
                               marker = pt.dat[[marker]])
    #row.names(marker.pt.df) <- pt.dat$metabric_id
    names(marker.pt.df)[2] <- pt
    #View(marker.pt.df)
    res.ls <- append(res.ls,list(marker.pt.df))
  }
  
  # assign and name dataframe for that marker
  final.dat <- Reduce(function(x, y) merge(x, y, by = "metabric_id", all = TRUE), res.ls)
  final.dat[,2:ncol(final.dat)] <- lapply(final.dat[,2:ncol(final.dat)], as.numeric)
  
  # plot
  boxplot(final.dat[2:ncol(final.dat)],las=2,ylab=marker,main=marker)
}

#######################################################################
# boxplots with
#######################################################################





# save plots
dev.off()