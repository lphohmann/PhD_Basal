# Script: PAM50 correlations in SCAN-B
# Author: Lennart Hohmann
# Date: 09.12.2024
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
# output paths
plot.file <- paste0(output.path,cohort,"_FOXA1vsBasalCentroid.pdf")
#txt.file <- paste0(output.path,cohort,"_i.txt")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
plot.parameters <- list() # object to store parameters to plot base R plots again later
#txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

# load sampleIDs
sampleIDs <- loadRData(infile.1)[["Lund.tumour.id"]]

c1.ids <- c("S004527", "S001534", "S000006", "S003591", "S000299", "S005062")
#setdiff(sampleIDs,c1.ids)
c2.ids <- c("S000037","S000327","S001023","S001318","S003109",
            "S003433","S004629","S000939","S002552","S004454") 

# annotation 
anno <- loadRData(infile.2)
anno <- anno[anno$Follow.up.cohort == TRUE,]
anno <- anno[anno$Sample %in% c(c2.ids,c1.ids), ]

# load data
corr.data <- loadRData(infile.3)
corr.data <- corr.data[corr.data$Assay %in% anno$GEX.assay,]
corr.data <- corr.data[!duplicated(corr.data$Assay),]
corr.data$Sample <- anno$Sample[match(corr.data$Assay,anno$GEX.assay)]


gex.data <- loadRData(infile.4)
gex.data <- gex.data["FOXA1",c(c1.ids,c2.ids)]

comb.dat <- as.data.frame(t(gex.data))
comb.dat$meanBasal <- corr.data$meanBasal[match(row.names(comb.dat),corr.data$Sample)]
#######################################################################
# Plot
#######################################################################
# save plots
pdf(file = plot.file, onefile = TRUE) 
par(mfrow = c(2, 2))

# LumA vs LumB
#plot(x=corr.data$meanLumA,y=corr.data$meanLumB, pch = 16)
#boxplot(corr.data[c("meanBasal","meanHer2","meanLumA","meanLumB","meanNormal")])

#######################################################################

bp <- boxplot(list(corr.data[which(corr.data$Sample %in% c1.ids), ]$meanBasal, 
                   corr.data[which(corr.data$Sample %in% c2.ids), ]$meanBasal),
        ylab="Basal-like centroid correlation",
        xlab="shore/CGI CpGs methylation cluster",
        col=c("#ba0606","#2176d5"),
        main="Basal centroid vs shoreCpG methyl.")
axis(3,at=1:length(bp$n),labels=bp$n)

bp <- boxplot(list(as.numeric(gex.data[, c1.ids]), 
                   as.numeric(gex.data[, c2.ids])),
              ylab="FOXA1 mRNA expression",
              xlab="shore/CGI CpGs methylation cluster",
              col=c("#ba0606","#2176d5"),
              main="FOXA1 expr. by shoreCpG methyl.")
axis(3,at=1:length(bp$n),labels=bp$n)

bp <- plot(x=comb.dat$FOXA1, 
          y=comb.dat$meanBasal,
          xlab="FOXA1 mRNA expression",
          ylab="Basal-like centroid correlation",
          main="FOXA1 expr. vs Basal centroid")
dev.off()

# also plot esr1
