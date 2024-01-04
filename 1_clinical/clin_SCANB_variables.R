# Script: Create table with clinicopathological variables comparisons between subtypes in the SCANB cohort (rs1, rel4) 
# Author: Lennart Hohmann
# Date: 01.01.2024
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
data.path <- "./data/SCANB/1_clinical/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
outfile.1 <- paste0(output.path,cohort,"_clinpath_table1.xlsx")
outfile.2 <- paste0(output.path,cohort,"_clinpath_table2.xlsx")
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
sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]

# select relevant variables
anno <- loadRData(infile.2) %>% 
  # sample selection
  filter(Follow.up.cohort==TRUE) %>%
  filter(Sample %in% unname(unlist(sampleIDs))) %>% 
  #
  mutate(LN_Bosch = case_when(LNstatus_Bosch == "N0" ~ 0, 
                              LNstatus_Bosch == "N+" ~ 1)) %>% 
  dplyr::select(
    Sample, 
    # table 1: N, Age, size, PR, NHG, LN, % review, H2low
    NCN.PAM50, Age, PR, Size.mm, NHG, LN, HER2_Low,
    # table 2: N, Age, size, PR, NHG, LN, H2low
    Bosch_RS1, TumSize_Bosch, NHG_Bosch, LN_Bosch) %>% 
  dplyr::rename(
    sampleID=Sample,
    PAM50=NCN.PAM50,
    Size=Size.mm,
    Size_Bosch = TumSize_Bosch)

# data type correction
anno$Bosch_RS1[is.na(anno$Bosch_RS1)] <- 0
anno$NHG_Bosch <- as.factor(anno$NHG_Bosch)
anno$LN_Bosch <- as.factor(anno$LN_Bosch)
anno$PR[anno$PR == ""] <- NA
anno$PAM50 <- as.factor(anno$PAM50)
anno$NHG[anno$NHG == ""] <- NA
anno$NHG <- as.factor(anno$NHG)
anno$LN <- as.factor(anno$LN)
anno$HER2_Low <- as.factor(anno$HER2_Low)
anno$Age <- as.numeric(anno$Age)
anno$Size <- as.numeric(anno$Size)

# sample names
Basal.ids <- unname(unlist(sampleIDs["ERpHER2n_Basal"]))
#Her2.ids <- anno %>% filter(PAM50=="Her2") %>% pull(sampleID)
LumA.ids <- unname(unlist(sampleIDs["ERpHER2n_LumA"]))
LumB.ids <- unname(unlist(sampleIDs["ERpHER2n_LumB"]))
anno <- anno %>% column_to_rownames(var = "sampleID")

# storing format: excel document with sheets for each variable
# keep column names consitent so later I just need to copy & paste tables together
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

age.data <- anno %>% 
  filter(!is.na(Age)) %>% 
  dplyr::select(PAM50,Age)

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

size.data <- anno %>% 
  filter(!is.na(Size)) %>% 
  dplyr::select(PAM50,Size)

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

nhg.data <- anno %>% filter(!is.na(NHG)) %>% 
  dplyr::select(PAM50,NHG)
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

ln.data <- anno %>% filter(!is.na(LN)) %>% 
  dplyr::select(PAM50,LN)
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
  
h2low.data <- anno %>% filter(!is.na(HER2_Low)) %>% 
  dplyr::select(PAM50,HER2_Low)
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

pr.data <- anno %>% filter(!is.na(PR)) %>% 
  dplyr::select(PAM50,PR)

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
# % reviewed 
#######################################################################

ct <- table(anno$PAM50,anno$Bosch_RS1)
variable <- c("N.non_reviewed","N.reviewed")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)


rev.count <- as.data.frame(cbind(variable,
                             basal.counts,basal.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

names(rev.count) <- col.names
rev.count

#######################################################################
# filter anno to only include review data
anno <- anno %>% filter(Bosch_RS1 == 1)

#######################################################################
# PAM50 count in review
#######################################################################

# matrix
pam50.count.rs1 <- as.data.frame(matrix(ncol=length(col.names)))
names(pam50.count.rs1) <- col.names
pam50.count.rs1$Variable <- "N.RS1"

t <- table(anno$PAM50)

# add count data
pam50.count.rs1$`Basal(ref)`[1] <- unname(t["Basal"])
pam50.count.rs1$LumA[1] <- unname(t["LumA"])
pam50.count.rs1$LumB[1] <- unname(t["LumB"])

# add % data
pam50.count.rs1$`Basal.%`[1] <- round((unname(t["Basal"])/sum(t))*100,2)
pam50.count.rs1$`LumA.%`[1] <- round((unname(t["LumA"])/sum(t))*100,2)
pam50.count.rs1$`LumB.%`[1] <- round((unname(t["LumB"])/sum(t))*100,2)

pam50.count.rs1

#######################################################################
# LN Bosch
#######################################################################

ln.data <- anno %>% filter(!is.na(LN_Bosch)) %>% 
  dplyr::select(PAM50,LN_Bosch)
#table(anno$PAM50,is.na(anno$LN_Bosch))

ct <- table(ln.data$PAM50,ln.data$LN_Bosch)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nLN.RS1\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("LN.N0.RS1","LN.N+.RS1")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

ln.df.rs1 <- as.data.frame(cbind(variable,
                             basal.counts,basal.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

names(ln.df.rs1) <- col.names
ln.df.rs1

#######################################################################
# Size Bosch 
#######################################################################

size.data <- anno %>% 
  filter(!is.na(Size_Bosch)) %>% 
  dplyr::select(PAM50,Size_Bosch)

#table(anno$PAM50,is.na(anno$Size_Bosch))

# matrix
size.df.rs1 <- as.data.frame(matrix(ncol=length(col.names)))
names(size.df.rs1) <- col.names
size.df.rs1$Variable <- "Size.RS1"

res.luma <- wilcox.test(size.data[Basal.ids,"Size_Bosch"],
                        size.data[LumA.ids,"Size_Bosch"])
res.lumb <- wilcox.test(size.data[Basal.ids,"Size_Bosch"],
                        size.data[LumB.ids,"Size_Bosch"])

txt.out <- append(txt.out, c("\nSize.RS1\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

means <- tapply(size.data$Size_Bosch, size.data$PAM50, mean)
sds <- tapply(size.data$Size_Bosch, size.data$PAM50, sd)

size.df.rs1$`Basal(ref)` <- unname(means["Basal"])
size.df.rs1$LumA <- unname(means["LumA"])
size.df.rs1$LumB <- unname(means["LumB"])

size.df.rs1$`Basal.%` <- unname(sds["Basal"])
size.df.rs1$`LumA.%` <- unname(sds["LumA"])
size.df.rs1$`LumB.%` <- unname(sds["LumB"])
size.df.rs1

#######################################################################
# NHG Bosch
#######################################################################

nhg.data <- anno %>% filter(!is.na(NHG_Bosch)) %>% 
  dplyr::select(PAM50,NHG_Bosch)
#table(anno$PAM50,is.na(anno$NHG_Bosch))

ct <- table(nhg.data$PAM50,nhg.data$NHG_Bosch)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nNHG.RS1\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("NHG1.RS1","NHG2.RS1","NHG3.RS1")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

nhg.df.rs1 <- as.data.frame(cbind(variable,
                              basal.counts,basal.perc,
                              luma.counts,luma.perc,
                              lumb.counts,lumb.perc))

names(nhg.df.rs1) <- col.names

nhg.df.rs1

#######################################################################
# Age for review
#######################################################################

age.data <- anno %>% 
  filter(!is.na(Age)) %>% 
  dplyr::select(PAM50,Age)

#table(anno$PAM50,is.na(anno$Age))

# matrix
age.df.rs1 <- as.data.frame(matrix(ncol=length(col.names)))
names(age.df.rs1) <- col.names
age.df.rs1$Variable <- "Age.RS1"

res.luma <- wilcox.test(age.data[Basal.ids,"Age"],age.data[LumA.ids,"Age"])
res.lumb <- wilcox.test(age.data[Basal.ids,"Age"],age.data[LumB.ids,"Age"])

txt.out <- append(txt.out, c("\nAge.RS1\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

means <- tapply(age.data$Age, age.data$PAM50, mean)
sds <- tapply(age.data$Age, age.data$PAM50, sd)

age.df.rs1$`Basal(ref)` <- unname(means["Basal"])
age.df.rs1$LumA <- unname(means["LumA"])
age.df.rs1$LumB <- unname(means["LumB"])

age.df.rs1$`Basal.%` <- unname(sds["Basal"])
age.df.rs1$`LumA.%` <- unname(sds["LumA"])
age.df.rs1$`LumB.%` <- unname(sds["LumB"])

age.df.rs1

#######################################################################
# PR for review
#######################################################################

pr.data <- anno %>% filter(!is.na(PR)) %>% 
  dplyr::select(PAM50,PR)

ct <- table(pr.data$PAM50,pr.data$PR)
res.luma <- fisher.test(ct[c("Basal","LumA"),])
res.lumb <- fisher.test(ct[c("Basal","LumB"),])

txt.out <- append(txt.out, c("\nPR.RS1\n",
                             "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.luma), "\n###########################################\n"))
txt.out <- append(txt.out, c(capture.output(res.lumb), "\n###########################################\n"))

variable <- c("PR.neg.RS1","PR.pos.RS1")
basal.counts <- unname(ct["Basal",])
luma.counts <- unname(ct["LumA",])
lumb.counts <- unname(ct["LumB",])
basal.perc <- round((basal.counts/sum(ct["Basal",]))*100,2)
luma.perc <- round((luma.counts/sum(ct["LumA",]))*100,2)
lumb.perc <- round((lumb.counts/sum(ct["LumB",]))*100,2)

pr.df.rs1 <- as.data.frame(cbind(variable,
                             basal.counts,basal.perc,
                             luma.counts,luma.perc,
                             lumb.counts,lumb.perc))

names(pr.df.rs1) <- col.names
pr.df.rs1

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
                         "%reviewed" = rev.count,
                         "HER2low" = h2low.df)

# rbind
table1 <- as.data.frame(do.call("rbind", sheet.list.table1))

write_xlsx(table1, outfile.1)

# review
sheet.list.table2 <- list("N" = pam50.count.rs1, 
                          "Age" = age.df.rs1,
                          "PR" = pr.df.rs1,
                          "Size" = size.df.rs1,
                          "NHG" = nhg.df.rs1,
                          "LN" = ln.df.rs1,
                          "HER2low" = h2low.df)

# rbind
table2 <- as.data.frame(do.call("rbind", sheet.list.table2))

write_xlsx(table2, outfile.2)

# save text
writeLines(txt.out, txt.file)