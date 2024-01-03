# general code snippets

################################################################################
# starting block all scripts
################################################################################

# Script: Aim
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
if (!require("pacman")) install.packages("pacman")
pacman::p_load(package1, package2)
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
infile.1 <- "./data/SCANB/4_CN/processed/Segment_CN_states.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_name.pdf")
txt.file <- paste0(output.path,cohort,"_name.txt")
outfile.1 <- paste0(data.path,"name.RData")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

################################################################################
# save plots
################################################################################

plot.list <- append(plot.list, list(plot))
pdf(file = plot.file, onefile = TRUE, height = 5, width = 5)
for(i in 1:length(plot.list)) {
  print(i)
  print(plot.list[[i]])
}
dev.off()

################################################################################
# save text
#################################################################################

txt.out <- append(txt.out, c(gene," statistics: Her2.dat vs. LumA.dat",
                    capture.output(res)))
writeLines(txt.out, txt.file)

################################################################################
# take execution time
################################################################################

start.time <- Sys.time()
# CODE
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

################################################################################
# compile files with pattern into one list
################################################################################

# ascat data file paths
temp <- list.files(
  path="./data/SCANB/4_CN/raw/to_lennart/",
  pattern="*.RData", full.names = TRUE, recursive = TRUE)
# load the files
ascat.files <- lapply(temp, loadRData)

################################################################################
# convert IDs
################################################################################

ascat.ids$Sample <- id.key$SENT.TUMOR[match(ascat.ids$Ascat.id, id.key$SENT.TUMOR.aliquot)] 

################################################################################
# loop progress bar
################################################################################

pb = txtProgressBar(min = 0, max = 10, initial = 0, style = 3)
for(i in 1:10) { 
  setTxtProgressBar(pb,i)
  close(pb)
}

################################################################################
# load excel file sheet
################################################################################

GSEA.HER2E.core.set <- openxlsx::read.xlsx("./output/supplementary_data/HER2n_pwenrichment_WikiPathway_2023_Human.xlsx",
                                           sheet="CoreDEGs")

################################################################################
# save as excel file with multiple sheets
################################################################################

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb=wb, sheetName="Core_gene_sets")
openxlsx::writeDataTable(wb=wb,sheet="Core_gene_sets",keepNA = TRUE,x=coreset.df)
openxlsx::addWorksheet(wb=wb, sheetName="GSEA.HER2E.core.set")
openxlsx::writeDataTable(wb=wb,sheet="GSEA.HER2E.core.set",x=GSEA.HER2E.core.set)
openxlsx::saveWorkbook(wb=wb,file="./output/supplementary_data/Manuscript_files/Table_S2.xlsx",
                       overwrite=TRUE)