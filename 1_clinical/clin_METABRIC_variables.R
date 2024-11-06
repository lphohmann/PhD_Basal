# Script: Create table with clinicopathological variables comparisons between subtypes in the METABRIC cohort  
# Author: Lennart Hohmann
# Date: 17.09.2024
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
source("./scripts/1_clinical/src/clin_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, matrixStats, 
               DescTools, openxlsx, writexl)
#-------------------
# set/create output directories
# for plots
output.path <- "./output/1_clinical/"
dir.create(output.path)
# for data
data.path <- "./data/METABRIC/1_clinical/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/METABRIC/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/METABRIC/1_clinical/raw/Merged_annotations.RData"
# output paths
outfile.1 <- paste0(output.path,cohort,"_clinpath_table1.xlsx")
#outfile.2 <- paste0(output.path,cohort,"_clinpath_table2.xlsx")
#plot.file <- paste0(output.path,cohort,"_SA.pdf")
txt.file <- paste0(output.path,cohort,"_variables.txt")
#-------------------
# storing objects 
#plot.list <- list() # object to store plots
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs.all <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB","ERpHER2n_HER2E")]
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# select relevant variables
anno <- loadRData(infile.2) 
anno <- anno[which(anno$METABRIC_ID %in% unname(unlist(sampleIDs))),]
# Create HER2_Low column based on conditions for HER2_IHC_status and HER2_SNP6_state
anno$HER2_Low <- ifelse(anno$HER2_IHC_status == 1, 1,
                        ifelse(anno$HER2_IHC_status == 2 & anno$HER2_SNP6_state != "GAIN", 1, NA))

# Create LN column based on lymph_nodes_positive
anno$LN <- ifelse(anno$lymph_nodes_positive > 0, 1, 0)

# Rename columns
colnames(anno)[colnames(anno) == "METABRIC_ID"] <- "sampleID"
colnames(anno)[colnames(anno) == "Grade"] <- "NHG"
colnames(anno)[colnames(anno) == "TumSize"] <- "Size"
colnames(anno)[colnames(anno) == "IC10_Rueda"] <- "IC10"
colnames(anno)[colnames(anno) == "PR.Expr"] <- "PR"

# Select relevant columns
anno <- anno[, c("sampleID", "PAM50", "Age", "Size", "NHG", "LN", "HER2_Low", "PR", "IC10")]

# replace na values with 0
anno$HER2_Low[is.na(anno$HER2_Low)] <- 0

# order IC10 
anno$IC10 <- factor(anno$IC10, levels=c("1","2","3","4ER+","6","7","8","9","10"))

# data type correction
anno$PR <- as.factor(anno$PR)
anno$PAM50 <- as.factor(anno$PAM50)
anno$NHG <- as.factor(anno$NHG)
anno$LN <- as.factor(anno$LN)
anno$HER2_Low <- as.factor(anno$HER2_Low)

# sample names
Basal.ids <- unname(unlist(sampleIDs["ERpHER2n_Basal"]))
#Her2.ids <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.ids <- unname(unlist(sampleIDs["ERpHER2n_LumA"]))
LumB.ids <- unname(unlist(sampleIDs["ERpHER2n_LumB"]))
rownames(anno) <- anno$sampleID
anno <- anno[, !(names(anno) %in% "sampleID")]

# storing format: excel document with sheets for each variable
# keep column names consistent so later I just need to copy & paste tables together
col.names <- c("Variable","Basal(ref)","Basal.%","LumA","LumA.%","LumB","LumB.%")

########################################################################
# 1. data vars (not review)
# table 1: N, Age, size, PR, NHG, LN, % review, H2low (review)

#######################################################################
# PAM50 count
#######################################################################

# matrix
pam50.count <- as.data.frame(matrix(ncol=length(col.names)))
names(pam50.count) <- col.names
pam50.count$Variable <- "N"

t <- table(anno$PAM50)

# add count data
pam50.count$`Basal(ref)`[1] <- unname(t["Basal"])
pam50.count$LumA[1] <- unname(t["LumA"])
pam50.count$LumB[1] <- unname(t["LumB"])

# add % data
pam50.count$`Basal.%`[1] <- round((unname(t["Basal"])/sum(t))*100,2)
pam50.count$`LumA.%`[1] <- round((unname(t["LumA"])/sum(t))*100,2)
pam50.count$`LumB.%`[1] <- round((unname(t["LumB"])/sum(t))*100,2)

pam50.count

#######################################################################
# Age comparison
#######################################################################
# note: I store the sd for continuous vars in the % column

age.data <- subset(anno, !is.na(Age), select = c(PAM50, Age))

#table(anno$PAM50,is.na(anno$Age))

# matrix
age.df <- as.data.frame(matrix(ncol=length(col.names)))
names(age.df) <- col.names
age.df$Variable <- "Age"

res.luma <- wilcox.test(age.data[Basal.ids,"Age"],age.data[LumA.ids,"Age"])
res.lumb <- wilcox.test(age.data[Basal.ids,"Age"],age.data[LumB.ids,"Age"])
#res.luma <- t.test(age.data[Basal.ids,"Age"],age.data[LumA.ids,"Age"])
#res.lumb <- t.test(age.data[Basal.ids,"Age"],age.data[LumB.ids,"Age"])

txt.out <- append(txt.out, c("\nAge\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

means <- tapply(age.data$Age, age.data$PAM50, mean)
sds <- tapply(age.data$Age, age.data$PAM50, sd)

age.df$`Basal(ref)` <- unname(means["Basal"])
age.df$LumA <- unname(means["LumA"])
age.df$LumB <- unname(means["LumB"])

age.df$`Basal.%` <- unname(sds["Basal"])
age.df$`LumA.%` <- unname(sds["LumA"])
age.df$`LumB.%` <- unname(sds["LumB"])

age.df

#######################################################################
# 4. Size
#######################################################################

size.data <- subset(anno, !is.na(Size), select = c(PAM50, Size))

#table(anno$PAM50,is.na(anno$Size))

# matrix
size.df <- as.data.frame(matrix(ncol=length(col.names)))
names(size.df) <- col.names
size.df$Variable <- "Size"

res.luma <- wilcox.test(size.data[Basal.ids,"Size"],size.data[LumA.ids,"Size"])
res.lumb <- wilcox.test(size.data[Basal.ids,"Size"],size.data[LumB.ids,"Size"])

txt.out <- append(txt.out, c("\nSize\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

means <- tapply(size.data$Size, size.data$PAM50, mean)
sds <- tapply(size.data$Size, size.data$PAM50, sd)

size.df$`Basal(ref)` <- unname(means["Basal"])
size.df$LumA <- unname(means["LumA"])
size.df$LumB <- unname(means["LumB"])

size.df$`Basal.%` <- unname(sds["Basal"])
size.df$`LumA.%` <- unname(sds["LumA"])
size.df$`LumB.%` <- unname(sds["LumB"])
size.df

#######################################################################
# 4. NHG (2x3 contingency table)
#######################################################################

nhg.data <- subset(anno, !is.na(NHG), select = c(PAM50, NHG))
#table(anno$PAM50,is.na(anno$NHG))

ct <- table(nhg.data$PAM50,nhg.data$NHG)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nNHG\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("NHG1","NHG2","NHG3")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

nhg.df <- as.data.frame(cbind(variable,
      basal.counts,basal.perc,
      luma.counts,luma.perc,
      lumb.counts,lumb.perc))

names(nhg.df) <- col.names

nhg.df

#######################################################################
# 4. LN (2x2 contingency table)
#######################################################################

ln.data <- subset(anno, !is.na(LN), select = c(PAM50, LN))
#table(anno$PAM50,is.na(anno$LN))

ct <- table(ln.data$PAM50,ln.data$LN)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nLN\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("LN.N0","LN.N+")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

ln.df <- as.data.frame(cbind(variable,
                             basal.counts,basal.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

names(ln.df) <- col.names
ln.df

#######################################################################
# HER2low frequency (2x2 contingency table)
#######################################################################
  
h2low.data <- subset(anno, !is.na(HER2_Low), select = c(PAM50, HER2_Low))
#table(anno$PAM50,is.na(anno$HER2_Low))

ct <- table(h2low.data$PAM50,h2low.data$HER2_Low)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nHER2low\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("HER2low.No","HER2low.Yes")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

h2low.df <- as.data.frame(cbind(variable,
                                basal.counts,basal.perc,
                                luma.counts,luma.perc,
                                lumb.counts,lumb.perc))

names(h2low.df) <- col.names
h2low.df


#######################################################################
# PR status (2x2 contingency table)
#######################################################################

pr.data <- subset(anno, !is.na(PR), select = c(PAM50, PR))

ct <- table(pr.data$PAM50,pr.data$PR)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nPR\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("PR.neg","PR.pos")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

pr.df <- as.data.frame(cbind(variable,
                             basal.counts,basal.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

names(pr.df) <- col.names
pr.df

#######################################################################
# save output
#######################################################################
  
# rel4+part of rs1
sheet.list.table1 <- list("N" = pam50.count, 
                         "Age" = age.df,
                         "PR" = pr.df,
                         "Size" = size.df,
                         "NHG" = nhg.df,
                         "LN" = ln.df,
                         "HER2low" = h2low.df)

# rbind
table1 <- as.data.frame(do.call("rbind", sheet.list.table1))

write_xlsx(table1, outfile.1)

# save text
writeLines(txt.out, txt.file)
