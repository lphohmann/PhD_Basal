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

# or
pdf(file = plot.file, onefile = TRUE)
for(i in 1:length(plot.list)) {
  grid::grid.newpage()
  grid::grid.draw(plot.list[[i]])
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

################################################################################
# save ggplots on one page
# save as one page
grob.list <- lapply(plot.list, function(x) {
  p <- ggplotGrob(x)
  return(p)
})
combined.plot <- arrangeGrob(grobs = grob.list, ncol = 2)
ggsave(plot.file.2, combined.plot, 
       width = 21, height = 29.7, units = "cm")

################################################################################
# save n plots per page

# Create a list of ggplotGrob objects
grob_list <- lapply(plot.list, function(x) {
  p <- ggplotGrob(x)
  return(p)
})

# Set the maximum number of plots per page and orientation
plots_per_page <- 6
ncol <- 2
nrow <- 3

# Calculate the number of pages needed
num_pages <- ceiling(length(grob_list) / plots_per_page)

# Create a multi-page PDF
pdf(plot.file.2, width = 8.27, height = 11.69)

# Loop through pages and save each page
for (page in 1:num_pages) {
  start_index <- (page - 1) * plots_per_page + 1
  end_index <- min(page * plots_per_page, length(grob_list))
  
  grob_list_page <- grob_list[start_index:end_index]
  
  # Arrange the ggplotGrob objects in a grid for each page
  grid_arranged_page <- grid.arrange(grobs = grob_list_page, 
                                     ncol = ncol, nrow = nrow)
  
  # Print the arranged grob to the current PDF page
  print(grid_arranged_page)
}

# Close the PDF device
dev.off()

################################################################################
################################################################################
# remove duplicates 

gex <- gex[!duplicated(gex$Hgnc_ID), ] # keep 1st row of each symbol



# convet gene symbols
library(biomaRt)

ensemble_to_hgnc <- function(ensemble.ids) {
  # Connect to the Ensembl database
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # Retrieve Hugo gene names using biomaRt
  res <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
               filters = "ensembl_gene_id",
               values = ensemble.ids,
               mart = mart,
               uniqueRows = TRUE)
  return(res)
}



### base r current plot saving

for(i in 1:3){
  # plot
  plot <- boxplot(list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
                col = color.palette, names = names(color.palette),
                ylab = "mRNA expression (log2)",
                main=gene)
  plot <- recordPlot() #records current plot
  plot.new()

  plot.list <- append(plot.list, list(plot))
}

# save plots
pdf(file = plot.file, onefile = TRUE) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()



#### base r saving mutiple plots per page by saving paramters for later plotting

## IMPORTANT: the list with data has to be in same order as color palette names!
plot.parameters <- list()
# plot
plot.par <- list(
  data = list(LumA=luma.dat,LumB=lumb.dat,Basal=basal.dat), 
  col = color.palette, 
  names = names(color.palette),
  ylab = "mRNA expression (log2)",
  main = gene)
# boxplot(plot_parameters$data, 
#         col = plot_parameters$col,
#         names = plot_parameters$names,
#         ylab = plot_parameters$ylab,
#         main = plot_parameters$main)
# plot <- recordPlot()
# plot.new()
plot.parameters <- append(plot.parameters, list(plot.par))
#plot.list <- append(plot.list, list(plot))



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
par(mfrow = c(1, 1))
dev.off()

# save text
writeLines(txt.out, txt.file)


#######################################################################

# get normal IDs in SCANB
id.key <- read_excel(
  "data/SCANB/3_genomic/raw/JVCimpression2022-11-09_ERpos_WGS_Batch1-3.xlsx") %>% 
  dplyr::select(c(SENT.TUMOR,SENT.TUMOR.aliquot,TUMOR.alias)) %>% 
  mutate(SENT.TUMOR.aliquot = gsub("\\.","_",SENT.TUMOR.aliquot))
id.df <- sapply(list(names(kat.scanb)), function(x) {sub("_vs_.*", "", x)}) %>% 
  as.data.frame() %>% 
  dplyr::rename(Old.id=1)
# 1. convert epb IDs to normal sample IDs
id.df$New.id <- id.key$SENT.TUMOR[match(
  id.df$Old.id, id.key$SENT.TUMOR.aliquot)] 
# 2. convert the other IDs to normal sample IDs
id.df$New.id.2 <- id.key$SENT.TUMOR[match(
  id.df$Old.id, id.key$TUMOR.alias)] # 2. convert the other IDs to normal sample IDs
id.df$New.id <- ifelse(is.na(id.df$New.id), id.df$New.id.2, id.df$New.id)
id.df$New.id.2 <- NULL
# name
names(kat.scanb) <- id.df$New.id

# get normal IDs in BASIS
names(kat.basis) <- sapply(list(names(kat.basis)), 
                           function(x) {sub("PD(.*?)[A-Za-z].*", "PD\\1", x)})

#######################################################################

# modify package functions

# modify
trace(ggforest, edit = TRUE)
# stop using updated function
untrace(ggforest)



# when defining function with required arguments
save(..., list = character(),
     file = stop("'file' must be specified"),
     ascii = FALSE, version = NULL)



#######
x1<-c(1,2,3)
x2<-c(4,2,3)
x3<-c(6,2,3)

x <- c(x1,x2,x3)
y <- list(x1,x2,x3)



dat_to_excel <- function(sheet.ls, file_name) {
  if (!require("writexl")) install.packages("writexl", dependencies = TRUE)
  library(writexl)
  # convert each list element into a df
  df_list <- lapply(names(sheet.ls), function(group_name) {
    data.frame(value = sheet.ls[[group_name]], group = group_name)
  })
  # rbind all 
  final_df <- do.call(rbind, df_list)
  # Write to an Excel file
  write_xlsx(final_df, file_name)
  cat("Data has been successfully written to", file_name, "\n")
}




save_dat <- function(dat.ls, name.vec) {
  # Check if the length of the list and the vector of names are the same
  if (length(dat.ls) != length(name.vec)) {
    stop("The number of names must match the number of data vectors.")
  }
  # Create a named list
  sheet.ls <- setNames(dat.ls, name.vec)
  # Return the named list
  return(sheet.ls)
}

# Define the groups to filter by
groups <- c("Her2", "LumA", "LumB")

# Use lapply to create data vectors and sample IDs for each group, with names
results <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.numeric(as.vector(dat[indices])),
         sample_ids = indices)
    }), 
  groups)

View(results)



# final

#### save source data ####
source_dat_path <- "./output/source_data/R_objects/"  
groups <- c("Her2", "LumA", "LumB")
results <- setNames(
  lapply(groups, function(group) {
    # Get the sample indices for the current group
    indices <- anno$Sample[anno$NCN.PAM50 == group]
    # Create the data vector and sample IDs
    list(data_vector = as.vector(dat[indices]),
         sample_ids = indices)
  }), 
  groups)
saveRDS(results, paste0(source_dat_path,"OncoDX_score.RData")) # Save 
#### exported source data ####


x <- TNBC_sampleIDs.metabric[c("TNBC_Basal","TNBC_NonBasal")]
lapply(x,length)
sum(unlist(lapply(x,length)))
(unlist(lapply(x,length)) / sum(unlist(lapply(x,length)))) * 100
round((unlist(lapply(x,length)) / sum(unlist(lapply(x,length)))) * 100)
