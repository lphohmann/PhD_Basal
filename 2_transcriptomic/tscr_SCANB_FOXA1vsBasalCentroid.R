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


#######################################################################
#######################################################################
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
library(VennDiagram)

foxa1.meth <- c("S004527", "S001534", "S000006", 
                "S003591", "S000299", "S005062")
foxa1.unmeth <- c("S000037","S000327","S001023",
                  "S001318","S003109","S003433",
                  "S004629","S000939","S002552","S004454") 
foxc1.meth <- c("S003109","S002552","S000037",
                "S001023","S000939","S000327","S003433")
foxc1.unmeth <- c("S004454","S003591","S005062","S000006",
                  "S000299","S001318","S004629","S001534")

intersect(foxa1.meth,foxc1.meth) #0
intersect(foxa1.meth,foxc1.unmeth) #5
intersect(foxa1.unmeth,foxc1.unmeth) #3
intersect(foxa1.unmeth,foxc1.meth) #7 # all meth foxc1 are foxa1 unmeth

# Create the Venn diagrams and display them
par(mfrow = c(2, 2)) # Arrange the plots in a 2x2 grid

# FoxA1 Meth vs. FoxC1 Meth
venn1 <- venn.diagram(
  x = list(FoxA1_Meth = foxa1.meth, FoxC1_Meth = foxc1.meth),
  category.names = c("FoxA1 Meth", "FoxC1 Meth"),
  filename = NULL,
  main = "FoxA1 Meth vs FoxC1 Meth"
)
grid.newpage()
grid.draw(venn1)

# FoxA1 Meth vs. FoxC1 Unmeth
venn2 <- venn.diagram(
  x = list(FoxA1_Meth = foxa1.meth, FoxC1_Unmeth = foxc1.unmeth),
  category.names = c("FoxA1 Meth", "FoxC1 Unmeth"),
  filename = NULL,
  main = "FoxA1 Meth vs FoxC1 Unmeth"
)
grid.newpage()
grid.draw(venn2)

# FoxA1 Unmeth vs. FoxC1 Unmeth
venn3 <- venn.diagram(
  x = list(FoxA1_Unmeth = foxa1.unmeth, FoxC1_Unmeth = foxc1.unmeth),
  category.names = c("FoxA1 Unmeth", "FoxC1 Unmeth"),
  filename = NULL,
  main = "FoxA1 Unmeth vs FoxC1 Unmeth"
)
grid.newpage()
grid.draw(venn3)

# FoxA1 Unmeth vs. FoxC1 Meth
venn4 <- venn.diagram(
  x = list(FoxA1_Unmeth = foxa1.unmeth, FoxC1_Meth = foxc1.meth),
  category.names = c("FoxA1 Unmeth", "FoxC1 Meth"),
  filename = NULL,
  main = "FoxA1 Unmeth vs FoxC1 Meth"
)
grid.newpage()
grid.draw(venn4)

#######################################################################
# Immune response vs methylation cluster for FOXA1?

# get IR score for foxa1 clusters
mg.scores <- loadRData("./data/SCANB/2_transcriptomic/processed/Metagene_scores_All.RData")
ir.scores <- mg.scores["IR",] 
foxa1.meth.ir <- as.numeric(ir.scores[foxa1.meth])
foxa1.unmeth.ir <- as.numeric(ir.scores[foxa1.unmeth])

# Create a list of the two groups
data_list <- list(
  "FoxA1 Meth" = foxa1.meth.ir,
  "FoxA1 Unmeth" = foxa1.unmeth.ir
)

# Create the boxplot
bp <- boxplot(
  data_list,
  main = "Comparison of IR Scores for FoxA1 Groups",
  xlab = "Group",
  ylab = "IR Score",
  col = c("lightblue", "lightgreen"),
  border = "black"
)
axis(3,at=1:length(bp$n),labels=bp$n)


dev.off()