# Function definitions for transcriptomic data analyses

################################################################################
# functions
################################################################################

# convert ensemble  gene ids to hugo
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


get_stats <- function(vec) {
  mean <- mean(vec)
  median <- median(vec)
  sd <- sd(vec)
  return(c("mean"=mean, "median"=median, "sd"=sd))
}

