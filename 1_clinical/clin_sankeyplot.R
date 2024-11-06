# Script: SCANB PAM50 Sankey plot, puts out to desktop
# Author: Lennart Hohmann
# Date: 27.08.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/1_clinical/src/clin_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               ggplot2,
               remotes)
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
#-------------------
# input paths
infile.1 <- "./data/SCANB/0_GroupSamples/ERpHER2n_sampleIDs.RData"
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
plot.file <- "../../Desktop/sankpeyplot_pam50.pdf"
#-------------------
# storing objects 
plot.list <- list() # object to store plots

#######################################################################
# load data
#######################################################################

color.palette <- loadRData(infile.2)
color.palette["unclassified"] <- "#bdbdbd"

dat <- loadRData(infile.3) %>% 
  filter(Follow.up.cohort==TRUE) %>% 
  filter(NCN.PAM50 != "Normal") %>% 
  dplyr::select(c("Sample","NCN.PAM50","ClinGroup"))
dat$ClinGroup[grepl("ERpHER2n", dat$ClinGroup)] <- "ERpHER2n"
dat$ClinGroup[dat$ClinGroup==""] <- "not available"

plot.dat <- dat %>% 
  make_long('NCN.PAM50', 'ClinGroup')  %>% 
  mutate(node = factor(node, levels = c("unclassified","Basal","Her2",
                                        "LumB","LumA","not available",
                                        "TNBC","HER2pERn","HER2pERp","ERpHER2n")))

ggplot(plot.dat,aes(x = x,
           next_x = next_x,
           node = node,
           next_node = next_node,
           fill = factor(node),
           label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = 1) +
  geom_sankey_label(size = 2, color = 1, fill = "white") + 
  scale_fill_manual(values = color.palette) +
  theme_sankey() +
  theme(legend.position = "none") +
  labs(x = NULL)

ggsave(plot.file)
