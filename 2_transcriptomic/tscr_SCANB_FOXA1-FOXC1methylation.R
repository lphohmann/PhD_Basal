# Script: Follow-up observed heterogeneity in FOXA1/FOXC1 promoter region methylation patterns 
# Author: Lennart Hohmann
# Date: 19.02.2025
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
output.path <- "./output/2_transcriptomic/"
dir.create(output.path)
# for data
data.path <- "./data/SCANB/2_transcriptomic/processed/"
dir.create(data.path)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2nBasal_WGS_sampleIDs.RData"
infile.2 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/pam50ann_correlations_summarizedByMean_REL4.RData"
infile.4 <- "./data/SCANB/2_transcriptomic/processed/ERp_LogScaled_gex.RData"
infile.5 <- "./data/SCANB/2_transcriptomic/processed/Metagene_scores_All.RData"
infile.6 <- "./data/SCANB/3_WGS/raw/2024_02_14_hrdetect_refsig_params_low_burden_sv_accounted_for.csv"
# output paths
plot.file <- paste0(output.path,cohort,"_FOXA1-FOXC1methylation.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# clusters
foxa1.A <- c("S004527", "S001534", "S000006",
             "S003591", "S000299", "S005062", 
             "S004454", "S004629")
foxa1.B <- c("S000037","S000327","S001023",
             "S001318","S003109","S003433",
             "S000939","S002552")
foxc1.A <- c("S003109","S002552","S000037","S001023",
             "S000939","S000327","S003433","S004454",
             "S003591","S005062","S000006","S000299")
foxc1.B <- c("S001318","S004629","S001534")
meth.clust <- list("FOXA1"=list("A"=foxa1.A,"B"=foxa1.B),
                   "FOXC1"=list("A"=foxc1.A,"B"=foxc1.B))

# load sampleIDs of wgs basal cases
sampleIDs <- loadRData(infile.1)[["Lund.tumour.id"]]

# annotation
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% sampleIDs, ]

# load pam50 centroid corr 
corr.data <- loadRData(infile.3)
corr.data <- corr.data[corr.data$Assay %in% anno$GEX.assay,]
corr.data <- corr.data[!duplicated(corr.data$Assay),]
corr.data$Sample <- anno$Sample[match(corr.data$Assay,anno$GEX.assay)]

gex.data <- loadRData(infile.4)
gex.data <- gex.data[sampleIDs]

# metagene scores
mg.scores <- loadRData(infile.5)
mg.scores <- mg.scores[sampleIDs]

#comb.dat <- as.data.frame(t(gex.data))
#comb.dat$meanBasal <- corr.data$meanBasal[match(row.names(comb.dat),corr.data$Sample)]

# hrd dat
hrd.dat <- read.table(infile.6, sep = ",", header = TRUE)
hrd.dat$Lund.tumour.id <- gsub("\\..*", "", hrd.dat$Lund.tumour.id)
hrd.dat <- hrd.dat[hrd.dat$Lund.tumour.id %in% sampleIDs,]
hrd.dat$HRDetect <- ifelse(hrd.dat$Probability >= 0.7,"HRD-high","HRD-low")

#######################################################################
# - Cluster overlap
#######################################################################

intersect(foxa1.A,foxc1.A) #0
intersect(foxa1.A,foxc1.B) #5
intersect(foxa1.B,foxc1.B) #3
intersect(foxa1.B,foxc1.A) #7 # all meth foxc1 are foxa1 unmeth

#######################################################################
# Plots
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))

# loop genes
for (fox in c("FOXA1","FOXC1")) {
  fox.dat <- meth.clust[[fox]]
  
  #######################################################################
  # 1. Clusters vs Basal centroid
  #######################################################################
  
  bp <- boxplot(list(corr.data[which(corr.data$Sample %in% fox.dat$A), ]$meanBasal, 
                     corr.data[which(corr.data$Sample %in% fox.dat$B), ]$meanBasal),
          ylab="Basal-like centroid correlation",
          xlab="shore/CGI CpGs methylation cluster",
          names = c("A","B"),
          col = c("#fc5603", "#edb374"),#adjustcolor(c("#fc5603", "#edb374"), alpha.f = 0.8),  # make boxes transparent
          #border = "gray40",
          #boxwex = 0.3,   # thinner box width
          main=paste0(fox,": Basal centroid vs MethClust"))
  axis(3,at=1:length(bp$n),labels=bp$n)
  
  stripchart(
    list(corr.data[which(corr.data$Sample %in% fox.dat$A), ]$meanBasal, 
         corr.data[which(corr.data$Sample %in% fox.dat$B), ]$meanBasal),
    method = "jitter",
    vertical = TRUE,pch = 16,cex=2,col = 1, add = TRUE)
  
  # bp <- plot(x=comb.dat$FOXA1, 
  #           y=comb.dat$meanBasal,
  #           xlab="FOXA1 mRNA expression",
  #           ylab="Basal-like centroid correlation",
  #           main="FOXA1 expr. vs Basal centroid")
  
  #######################################################################
  # 2. Clusters vs HRD
  #######################################################################
  
  bp <- boxplot(list(hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$A), ]$Probability, 
                     hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$B), ]$Probability),
                ylab="HRD probability",
                xlab="shore/CGI CpGs methylation cluster",
                names = c("A","B"),
                col=c("#fc5603","#edb374"),
                main=paste0(fox,": HRD prob vs MethClust"))
  axis(3,at=1:length(bp$n),labels=bp$n)
  
  stripchart(
    list(hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$A), ]$Probability, 
         hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$B), ]$Probability),
    method = "jitter",vertical = TRUE,pch = 16,cex=2,col = 1, add = TRUE)
  
  # Count occurrences of "high" and "low" for each group
  A_counts <- table(hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$A), ]$HRDetect)
  B_counts <- table(hrd.dat[which(hrd.dat$Lund.tumour.id %in% fox.dat$B), ]$HRDetect)
  
  # Convert counts to percentages
  A_percent <- A_counts / sum(A_counts) * 100
  B_percent <- B_counts / sum(B_counts) * 100
  
  # Combine percentages into a matrix
  percent_matrix <- rbind(A_percent, B_percent=B_percent[seq(A_percent)])
  percent_matrix[is.na(percent_matrix)] <- 0
  # Create stacked barplot (percentage format)
  bp <- barplot(t(percent_matrix), 
                beside=TRUE, # Stacked bars
                col=c("#fa43a1","#fccfe6"), 
                names.arg=c("A", "B"), 
                ylab="Percentage", 
                ylim= c(0,100),
                xlab="shore/CGI CpGs methylation cluster", 
                main=paste0(fox, ": HRD prob vs MethClust"),
                legend.text = c("HRD-High","HRD-Low"))
  
  lab <- c(A_counts,B_counts[seq(A_counts)])
  lab[is.na(lab)] <- 0
  
  axis(3,at=1:4,labels=c(A_counts,B_counts[seq(A_counts)]))
  
  #######################################################################
  # 3. Clusters vs mRNA expression (FOXA1,FOXC1,ESR1)
  #######################################################################
  
  for(gene in c("FOXA1","FOXC1","ESR1")) {
    bp <- boxplot(list(as.numeric(gex.data[gene, fox.dat$A]), 
                       as.numeric(gex.data[gene, fox.dat$B])),
                  ylab="mRNA expression",
                  xlab="shore/CGI CpGs methylation cluster",
                  names = c("A","B"),
                  col=c("#fc5603","#edb374"),
                  main=paste0(fox,": ",gene," vs MethClust"))
    axis(3,at=1:length(bp$n),labels=bp$n)
    stripchart(
      list(as.numeric(gex.data[gene, fox.dat$A]), 
           as.numeric(gex.data[gene, fox.dat$B])),
      method = "jitter",vertical = TRUE,pch = 16,cex=2,col = 1, add = TRUE)
  }
  
  #######################################################################
  # 4. Clusters vs metagene rank scores
  #######################################################################
  
  for(metagene in row.names(mg.scores)) {
    bp <- boxplot(list(as.numeric(mg.scores[metagene, fox.dat$A]), 
                       as.numeric(mg.scores[metagene, fox.dat$B])),
                  ylab="rank score",
                  xlab="shore/CGI CpGs methylation cluster",
                  names = c("A","B"),
                  col=c("#fc5603","#edb374"),
                  main=paste0(fox,": ",metagene," vs MethClust"))
    axis(3,at=1:length(bp$n),labels=bp$n)
    stripchart(
      list(as.numeric(mg.scores[metagene, fox.dat$A]), 
           as.numeric(mg.scores[metagene, fox.dat$B])),
      method = "jitter",vertical = TRUE,pch = 16,cex=2,col = 1, add = TRUE)
  }
}

#######################################################################
#######################################################################
dev.off()