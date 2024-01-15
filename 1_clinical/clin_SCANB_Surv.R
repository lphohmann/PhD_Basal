# Script: Survival analyses in the SCANB cohort based on IDFS (SCANB rs1), OS (SCANB rel4) 
# Author: Lennart Hohmann
# Date: 01.01.2024
#-------------------
# empty environment
rm(list=ls())
# set working directory to the project directory
setwd("~/PhD_Workspace/Project_Basal/")
# cohort
cohort <- "SCANB"
# plot also HER2E?
plot.HER2E <- "no" #no
#-------------------
# packages
source("./scripts/src/general_functions.R")
source("./scripts/1_clinical/src/clin_functions.R")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, 
               ggplot2,
               ggfortify,
               survival,
               survminer,
               grid,
               gridExtra,
               ggplotify,
               remotes)
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
infile.2 <- "./data/Parameters/color_palette.RData"
infile.3 <- "./data/SCANB/1_clinical/raw/Summarized_SCAN_B_rel4_NPJbreastCancer_with_ExternalReview_Bosch_data.RData"
# output paths
plot.file <- paste0(output.path,cohort,"_SA.pdf")
plot.file.2 <- paste0(output.path,cohort,"_SA_combined.pdf")
txt.file <- paste0(output.path,cohort,"_SA.txt")
#outfile.1 <- paste0(data.path,"SA.RData")
#-------------------
# storing objects 
plot.list <- list() # object to store plots
txt.out <- c() # object to store text output, if the output is not in string format use capture.output()

#######################################################################
# load data
#######################################################################

if (plot.HER2E == "yes") {
  sampleIDs <- loadRData(infile.1)
  color.palette <- loadRData(infile.2)
} else { 
  # load sampleIDs
  sampleIDs <- loadRData(infile.1)[c("ERpHER2n_Basal", "ERpHER2n_LumA", "ERpHER2n_LumB")]
  # load palette
  color.palette <- loadRData(infile.2)[c("LumA","LumB","Basal")]
  }

names(color.palette) <- paste0("PAM50=", names(color.palette))

# review 1 data
svdata.rs1 <- loadRData(infile.3) %>% 
  # sample selection
  filter(Follow.up.cohort==TRUE) %>%
  filter(Sample %in% unname(unlist(sampleIDs))) %>% 
  # variable selection and correction
  mutate(Treatment = case_when(Chemo_Bosch & ET_Bosch ~ "CE", !Chemo_Bosch & ET_Bosch ~ "E")) %>% 
  mutate(Bosch_RS1 = ifelse(!is.na(Treatment), 1, NA)) %>% 
  filter(!is.na(Bosch_RS1)) %>%
  dplyr::select(c("Sample","Treatment","Age","NCN.PAM50",
                  "IDFS_bin_Bosch","IDFS_Bosch",
                  "TumSize_Bosch","NHG_Bosch","LNstatus_Bosch")) %>%
  dplyr::rename(IDFS = IDFS_Bosch, 
                IDFSbin = IDFS_bin_Bosch,
                TumSize = TumSize_Bosch, 
                PAM50 = NCN.PAM50,
                Grade = NHG_Bosch, 
                LN = LNstatus_Bosch) %>%
  mutate(across(c(IDFS,IDFSbin,Age,TumSize), as.numeric)) %>%
  mutate(Grade = factor(Grade, levels = c("1","2","3"))) %>%
  mutate(PAM50 = relevel(factor(PAM50), ref = "LumA")) %>%
  mutate(LN = factor(LN, levels = c("N0","N+")))

# rel4 data
svdata.rel4 <- loadRData(infile.3) %>% 
  # sample selection
  filter(Follow.up.cohort==TRUE) %>%
  filter(Sample %in% unname(unlist(sampleIDs))) %>% 
  # variable selection and correction
  mutate(Treatment = case_when(TreatGroup=="ChemoEndo" ~ "CE",
                               TreatGroup=="Endo" ~ "E")) %>%
  dplyr::select(c("Sample","Treatment","Age","NCN.PAM50",
                  "OS","OSbin","DRFI","DRFIbin","Size.mm","NHG","LN")) %>%
  dplyr::rename(
    TumSize = Size.mm,
    PAM50 = NCN.PAM50,
    Grade = NHG) %>%
  mutate(LN = ifelse(LN > 0, "N+", "N0")) %>%
  mutate(across(c(OS,OSbin,Age,TumSize), as.numeric)) %>%
  mutate(Grade = factor(Grade, levels = c("1","2","3"))) %>%
  mutate(PAM50 = relevel(factor(PAM50), ref = "LumA")) %>%
  mutate(LN = factor(LN, levels = c("N0","N+")))

#######################################################################
# Part 1: IDFS
#######################################################################

# Fig A: IDFS; 2x (KM + uniCox) + 2x FP of mvCox

# outcome measure
OM <- "IDFS"
OMbin <- "IDFSbin"

# sub cohort
sub.cohort <- "RS1"

######################################
# Investigate the EC treatment group #
######################################

# data, surv object, and fit
EC.dat <- svdata.rs1 %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
EC.fit <- survminer::surv_fit(EC.surv~PAM50, data=EC.dat, 
                           conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; CT+ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: CT+ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(EC.dat$PAM50)),"\n\n"))
# add the median OM of censored patients
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for CT+ET censored patients = ",median.EC),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(EC.surv~PAM50, data=EC.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(EC.fit, 
           data = EC.dat,
           pval = TRUE, conf.int = FALSE,
           xlab = paste0(OM," (years)"), 
           break.x.by = 1,
           break.y.by = 0.1,
           ylab = paste0(OM," probability"),
           ylim = c(0,1),
           title = plot.title,
           legend = c(0.9,0.96),
           legend.title = "Subtypes",
           palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(EC.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=EC.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
         data=EC.dat,
         main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

########################################
# Investigate the Endo treatment group
########################################

# data, surv object, and fit
E.dat <- svdata.rs1 %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
E.fit <- survminer::surv_fit(E.surv~PAM50, data=E.dat, 
                           conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(E.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for ET censored patients = ",median.E),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(E.surv~PAM50, data=E.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(E.fit, 
                   data = E.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.96),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(E.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=E.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

#######################################################################
# Part 2: OS
#######################################################################

# - Fig B: OS; 2x (KM + uniCox) + 2x FP of mvCox

OM <- "OS"
OMbin <- "OSbin"

# sub cohort
sub.cohort <- "Rel4"

######################################
# Investigate the EC treatment group #
######################################

# data, surv object, and fit
EC.dat <- svdata.rel4 %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
EC.fit <- survminer::surv_fit(EC.surv~PAM50, data=EC.dat, 
                              conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; CT+ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: CT+ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(EC.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for CT+ET censored patients = ",median.EC),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(EC.surv~PAM50, data=EC.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(EC.fit, 
                   data = EC.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.3),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(EC.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=EC.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

########################################
# Investigate the Endo treatment group
########################################

# data, surv object, and fit
E.dat <- svdata.rel4 %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
E.fit <- survminer::surv_fit(E.surv~PAM50, data=E.dat, 
                             conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(E.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for ET censored patients = ",median.E),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(E.surv~PAM50, data=E.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(E.fit, 
                   data = E.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.96),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(E.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=E.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

#######################################################################
#######################################################################
#######################################################################
# Part 3: DFRI
#######################################################################

# - Fig C: DRFI; 2x (KM + uniCox) + 2x FP of mvCox

OM <- "DRFI"
OMbin <- "DRFIbin"

# sub cohort
sub.cohort <- "Rel4"

######################################
# Investigate the EC treatment group #
######################################

# data, surv object, and fit
EC.dat <- svdata.rel4 %>% 
  filter(Treatment == "CE")
EC.surv <- Surv(EC.dat[[OM]], EC.dat[[OMbin]])
EC.fit <- survminer::surv_fit(EC.surv~PAM50, data=EC.dat, 
                              conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; CT+ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: CT+ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(EC.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.EC <- median(EC.dat[which(EC.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for CT+ET censored patients = ",median.EC),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(EC.surv~PAM50, data=EC.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=EC.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(EC.fit, 
                   data = EC.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.3),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
#exp(confint(main.all)) # the problem is that grade2/3 have infinite CI bounds
main.all <- coxph(EC.surv~PAM50+Age+LN+TumSize, #+Grade
                  data=EC.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=EC.dat,
                 main=plot.title)

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

########################################
# Investigate the Endo treatment group
########################################

# data, surv object, and fit
E.dat <- svdata.rel4 %>% 
  filter(Treatment == "E")
E.surv <- Surv(E.dat[[OM]], E.dat[[OMbin]])
E.fit <- survminer::surv_fit(E.surv~PAM50, data=E.dat, 
                             conf.type="log-log") 
# label output
plot.title <- paste0("ERpHER2n; ",cohort,"; ",sub.cohort,"; ET; ",OM)
txt.out <- append(txt.out,
                  c(plot.title, "\n", 
                    paste0("Analyses with clinical endpoint: ",OM, " in treatment group: ET"),"\n###########################################\n"))
# subtype numbers
txt.out <- append(txt.out,
                  c(capture.output(table(E.dat$PAM50)),"\n###########################################\n"))
# add the median OM of censored patients
median.E <- median(E.dat[which(E.dat[[OMbin]]==0),][[OM]])
txt.out <- append(txt.out,c(paste0("Median ",OM, " for ET censored patients = ",median.E),"\n###########################################\n"))

##########################

# uv cox
main.pam50 <- coxph(E.surv~PAM50, data=E.dat)
res <- summary(main.pam50)
plot <- ggforest(main.pam50,data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

##########################

# KM plot
plot <- ggsurvplot(E.fit, 
                   data = E.dat,
                   pval = TRUE, conf.int = FALSE,
                   xlab = paste0(OM," (years)"), 
                   break.x.by = 1,
                   break.y.by = 0.1,
                   ylab = paste0(OM," probability"),
                   ylim = c(0,1),
                   title = plot.title,
                   legend = c(0.9,0.3),
                   legend.title = "Subtypes",
                   palette = color.palette)[["plot"]]

plot.list <- append(plot.list,list(plot))

##########################

# mv cox
main.all <- coxph(E.surv~PAM50+Age+LN+TumSize+Grade, 
                  data=E.dat) 
res <- summary(main.all)
plot <- ggforest(main.all,
                 data=E.dat,
                 main=plot.title) + theme_bw()

plot.list <- append(plot.list,list(plot))
txt.out <- append(txt.out,c(capture.output(res),"\n###########################################\n"))

#######################################################################
#######################################################################

# save plots
pdf(file = plot.file, onefile = TRUE) 
for (i in 1:length(plot.list)) {
  print(plot.list[[i]])
}
dev.off()

# save text
writeLines(txt.out, txt.file)



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
