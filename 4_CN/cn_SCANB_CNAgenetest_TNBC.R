# Script: TNBC Test significnace of CNA freqs for each gene, create file
# Author: Lennart Hohmann
# Date: 30.03.2024
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
infile.1 <- "./data/SCANB/4_CN/processed/CNA_genelevel_all.RData"
infile.3 <- "./data/SCANB/4_CN/processed/CNA_genelevel_TNBC.RData"
# output paths
outfile.1 <- paste0(data.path,"CNA_genetest_TNBC.RData")
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

ascat.df.all <- loadRData(infile.1)[["ascat.df.all"]]
ascat.df.all[,2:ncol(ascat.df.all)] <- lapply(ascat.df.all[,2:ncol(ascat.df.all)], as.numeric)
basal.ids <- loadRData(infile.1)[["subtype.samples"]][["Basal"]]
ascat.df.basal <- ascat.df.all[c("gene","chr","start","end",basal.ids)]

#load tnbc data
ascat.df.tnbc <- loadRData(infile.3)[["ascat.df.tnbc"]]

ascat.df.tnbc[,2:ncol(ascat.df.tnbc)] <- lapply(ascat.df.tnbc[,2:ncol(ascat.df.tnbc)], as.numeric)
tnbc.basal.ids <- loadRData(infile.3)[["subtype.samples"]][["TNBC.Basal"]]
tnbc.nonbasal.ids <- loadRData(infile.3)[["subtype.samples"]][["TNBC.NonBasal"]]

# only include common genes in ASCAT
common.genes <- intersect(ascat.df.tnbc$gene,ascat.df.basal$gene)
ascat.df.tnbc <- ascat.df.tnbc[ascat.df.tnbc$gene %in% common.genes, ]
ascat.df.basal <- ascat.df.basal[ascat.df.basal$gene %in% common.genes, ]
identical(ascat.df.tnbc$gene,ascat.df.basal$gene)
# order
ascat.df.tnbc <- ascat.df.tnbc[order(ascat.df.tnbc$gene), ]
ascat.df.basal <- ascat.df.basal[order(ascat.df.basal$gene), ]
identical(ascat.df.tnbc$gene,ascat.df.basal$gene)

# 
ascat.df.all <- cbind(ascat.df.basal,ascat.df.tnbc[5:ncol(ascat.df.tnbc)])

#

#######################################################################
# test genes Loss and Gain freqs
#######################################################################

TNBC.Basal.Loss.pval <- c()
TNBC.Basal.Gain.pval <- c()
TNBC.NonBasal.Loss.pval <- c()
TNBC.NonBasal.Gain.pval <- c()

pb = txtProgressBar(min = 0, max = length(ascat.df.all$gene), initial = 0, style = 3)
for(i in 1:length(ascat.df.all$gene)) {
  setTxtProgressBar(pb,i)
  # gene metadata
  gene <- ascat.df.all$gene[i]
  gene.anno <- ascat.df.all[ascat.df.all$gene==gene,
                            c("gene","chr","start","end")]
  
  # # gene CNA data with PAM50 sample annotation
  # gene.dat <- as.data.frame(
  #   t(ascat.df.all[ascat.df.all$gene==gene,
  #                  !names(ascat.df.all) %in% c("gene","chr","start","end")]))
  # gene.dat$sampleID <- rownames(gene.dat)
  # gene.dat$PAM50 <- ifelse(gene.dat$sampleID %in% basal.ids, "Basal",
  #                         ifelse(gene.dat$sampleID %in% luma.ids, "LumA",
  #                                ifelse(gene.dat$sampleID %in% lumb.ids, "LumB", NA)))
  # gene.dat$sampleID <- NULL
  # names(gene.dat)[1] <- "gene.CNA"
  
  #View(gene.dat)
  basal.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,basal.ids])
  tnbc.basal.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,tnbc.basal.ids])
  tnbc.nonbasal.dat <- unlist(ascat.df.all[ascat.df.all$gene==gene,tnbc.nonbasal.ids])

  comp.list <- list("TNBC.Basal"=tnbc.basal.dat,"TNBC.NonBasal"=tnbc.nonbasal.dat)
  
  for(j in 1:length(comp.list)) {
    comp.dat <- comp.list[[j]] 
    comp.group <- names(comp.list)[j]
    
    gain.tbl <- data.frame(basal=c(sum(basal.dat==1),sum(basal.dat!=1)), 
                           comp=c(sum(comp.dat==1),sum(comp.dat!=1)), 
                           row.names = c("gain","no_gain"))
    
    loss.tbl <- data.frame(basal=c(sum(basal.dat==-1),sum(basal.dat!=-1)), 
                           comp=c(sum(comp.dat==-1),sum(comp.dat!=-1)), 
                           row.names = c("loss","no_loss"))
    
    # test gain
    gain.pval <- fisher.test(gain.tbl)$p.value
    # test loss
    loss.pval <- fisher.test(loss.tbl)$p.value
    
    # save in vectors
    if (comp.group=="TNBC.Basal") {
      TNBC.Basal.Gain.pval[i] <- gain.pval
      TNBC.Basal.Loss.pval[i] <- loss.pval
    } else if (comp.group=="TNBC.NonBasal") {
      TNBC.NonBasal.Gain.pval[i] <- gain.pval
      TNBC.NonBasal.Loss.pval[i] <- loss.pval
    }
    
  }
  close(pb)
}



gene.test.df <- data.frame("gene"=ascat.df.all$gene, 
                           "TNBC.Basal.Gain.pval"=TNBC.Basal.Gain.pval,
                           "TNBC.Basal.Loss.pval"=TNBC.Basal.Loss.pval,
                           "TNBC.NonBasal.Gain.pval"=TNBC.NonBasal.Gain.pval,
                           "TNBC.NonBasal.Loss.pval"=TNBC.NonBasal.Loss.pval)


#save(gene.test.df, file= signif.genes)

# adjust pval
gene.test.df$TNBC.Basal.Gain.padj <- p.adjust(gene.test.df$TNBC.Basal.Gain.pval, method = "fdr") #fdr
gene.test.df$TNBC.NonBasal.Gain.padj <- p.adjust(gene.test.df$TNBC.NonBasal.Gain.pval, method = "fdr")
gene.test.df$TNBC.Basal.Loss.padj <- p.adjust(gene.test.df$TNBC.Basal.Loss.pval, method = "fdr")
gene.test.df$TNBC.NonBasal.Loss.padj <- p.adjust(gene.test.df$TNBC.NonBasal.Loss.pval, method = "fdr")
#View(gene.test.df)
length(gene.test.df[gene.test.df$TNBC.Basal.Gain.padj<=0.05,]$gene)
length(gene.test.df[gene.test.df$TNBC.NonBasal.Gain.padj<=0.05,]$gene)
length(gene.test.df[gene.test.df$TNBC.Basal.Loss.padj<=0.05,]$gene)
length(gene.test.df[gene.test.df$TNBC.NonBasal.Loss.padj<=0.05,]$gene)

save(gene.test.df, file= outfile.1)
